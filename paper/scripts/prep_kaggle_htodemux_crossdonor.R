# | Dataset                                 | Description                                                      | Use                                         |
#   | --------------------------------------- | ---------------------------------------------------------------- | ------------------------------------------- |
#   | **PBMC (10x Genomics CITE-seq)**        | Human peripheral blood mononuclear cells with RNA + 30–200 ADTs. | Core benchmark — train on RNA, predict ADT.

# FMLE validation setup (extension)

# | Validation Axis                         | Dataset                                               | Purpose                                                                                          |
# | --------------------------------------- | ----------------------------------------------------- | ------------------------------------------------------------------------------------------------ |
# | **(a) scRNA–ADT (standard)**            | PBMC, CBMC, MALT                                      | Show FMLE reduces bias & improves calibration vs scLinear and TotalVI.  


# https://www.kaggle.com/datasets/alexandervc/citeseq-scrnaseq-proteins-human-pbmcs-2019


#----------------------------------------#
#  Data preprocessing
#----------------------------------------#

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(clusterProfiler)
  library(IRanges)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(FMLE)
})

files <- c("_config.R", "preprocess_helpers_fmle.R", "preprocess_fmle.R")
paths <- file.path(here::here(), "paper", "scripts", files)
stopifnot(all(file.exists(paths)))
invisible(lapply(paths, source))

out_base <- file.path(cfg$out_root, "citeseq_v1")
out_ctp  <- file.path(cfg$out_root, "ctp")
dir.create(out_base, recursive=TRUE, showWarnings=FALSE)
dir.create(out_ctp,  recursive=TRUE, showWarnings=FALSE)

path_rna <- file.path(cfg$data_root, "citeseq_rna_counts.tsv")
path_adt <- file.path(cfg$data_root, "citeseq_adt_counts.tsv")
path_hto <- file.path(cfg$data_root, "citeseq_hto_counts.tsv")


## -------------------------------
## 1) Readers + barcode harmonization
## -------------------------------
read_feature_by_cell <- function(path){
  bcs <- scan(path, what = "", nlines = 1, quiet = TRUE, sep = "\t")
  if (length(bcs) && bcs[1] == "") bcs <- bcs[-1]
  dt <- fread(path, sep = "\t", header = FALSE, skip = 1)
  feats <- dt[[1]]; dt[[1]] <- NULL
  mat <- as.matrix(dt); storage.mode(mat) <- "double"
  rownames(mat) <- make.unique(as.character(feats))
  colnames(mat) <- make.unique(as.character(bcs))
  mat
}
normalize_barcodes <- function(v){
  v <- as.character(v)
  repeat { v2 <- sub("^[A-Za-z0-9]+_", "", v); if (identical(v2,v)) break; v <- v2 }
  v <- sub("[\\.-]1$", "", v); toupper(v)
}

rna_mat <- read_feature_by_cell(path_rna)      # genes × cells
adt_mat <- read_feature_by_cell(path_adt)      # ADTs  × cells
hto_mat <- read_feature_by_cell(path_hto)      # HTOs  × cells
stopifnot(ncol(rna_mat)>0, ncol(adt_mat)>0, ncol(hto_mat)>0)

bc_rna <- normalize_barcodes(colnames(rna_mat))
bc_adt <- normalize_barcodes(colnames(adt_mat))
bc_hto <- normalize_barcodes(colnames(hto_mat))
common <- Reduce(intersect, list(bc_rna, bc_adt, bc_hto))
if (!length(common)) stop("No overlapping cell barcodes after normalization.")

reindex <- function(M, bc_now, common){
  idx <- match(common, bc_now); if (anyNA(idx)) stop("Reindex failure.")
  M[, idx, drop = FALSE]
}
rna_mat <- reindex(rna_mat, bc_rna, common)
adt_mat <- reindex(adt_mat, bc_adt, common)
hto_mat <- reindex(hto_mat, bc_hto, common)
colnames(rna_mat) <- colnames(adt_mat) <- colnames(hto_mat) <- common
cat(sprintf("Aligned cells = %d  | genes=%d  ADTs=%d  HTOs=%d\n",
            length(common), nrow(rna_mat), nrow(adt_mat), nrow(hto_mat)))



seu <- CreateSeuratObject(counts = rna_mat, assay = "RNA", project = "PBMC10k_CITE_Kaggle")

common <- intersect(colnames(seu), colnames(adt_mat))
cat("common cells:", length(common), "\n")
stopifnot(length(common) > 0)
seu <- subset(seu, cells = common)
adt_mat <- adt_mat[, colnames(seu), drop=FALSE]

# Add ADT assay
seu[["ADT"]] <- CreateAssayObject(counts = adt_mat)
stopifnot(identical(colnames(seu[["RNA"]]), colnames(seu[["ADT"]])))


# ============================================================
# 2) QC on RNA + optional ADT-total outlier filter
# ============================================================
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

summary(seu$nFeature_RNA)
summary(seu$nCount_RNA)
summary(seu$percent.mt)

seu <- subset(
  seu,
  subset =
    nFeature_RNA >= 200  & nFeature_RNA <= 3500  &
    nCount_RNA   >= 500  & nCount_RNA   <= 30000 &
    percent.mt   <= 15
)

# Optional: remove extreme ADT-total outliers (winsor at 99.5%)
adt_counts0 <- GetAssayData(seu, assay="ADT", layer="counts")
cs <- Matrix::colSums(adt_counts0)
summary(cs)
thr <- as.numeric(quantile(cs, cfg$adt_total_q))
keep_cells <- names(cs)[cs <= thr]
seu <- subset(seu, cells = keep_cells)

stopifnot(identical(colnames(seu[["RNA"]]), colnames(seu[["ADT"]])))


# ============================================================
# 3) preprocess_fmle (this FINALIZES the cell set)
# ============================================================

seu_final <- preprocess_fmle(
  object = seu,
  annotation_selfCluster  = TRUE,
  remove_doublets = TRUE,
  low_qc_cell_removal = TRUE,
  integrate_data = FALSE,
  remove_empty_droplets = FALSE,
  resolution = 0.8,
  seed = 42,
  return_plots = FALSE,
  print_plots = FALSE,
  species = "Hs"
)

# seu_final <- scLinear::prepare_data(
#   seu,
#   integrate_data          = FALSE,
#   annotation_selfCluster  = TRUE,
#   remove_empty_droplets   = FALSE
# )

stopifnot(identical(colnames(seu_final[["RNA"]]), colnames(seu_final[["ADT"]])))

saveRDS(seu_final, file.path(out_base, "seu_final.rds"))


# ============================================================
# 4) Build X (HVG lognorm), Z (PCA), Y (raw ADT counts) from FINAL object
# ============================================================
seu_final <- NormalizeData(seu_final, assay="RNA", normalization.method="LogNormalize", scale.factor=cfg$scale_factor)
seu_final <- FindVariableFeatures(seu_final, assay="RNA", selection.method="vst", nfeatures=cfg$hvg_n, verbose=FALSE)
hvg <- VariableFeatures(seu_final)

seu_final <- ScaleData(seu_final, assay="RNA", features=hvg, verbose=FALSE)
ElbowPlot(seu_final, ndims = 50)
seu_final <- RunPCA(seu_final, assay="RNA", features=hvg, npcs=cfg$pca_npcs, verbose=FALSE)

Zfull <- Embeddings(seu_final, "pca")
gate <- cfg$gate_pcs
Z <- Zfull[, seq_len(min(gate, ncol(Zfull))), drop=FALSE]

X <- t(as.matrix(GetAssayData(seu_final, assay="RNA", layer="data")[hvg, , drop=FALSE]))
rownames(X) <- colnames(seu_final)

Y_counts <- GetAssayData(seu_final, assay="ADT", layer="counts")  # proteins x cells

# Align again (paranoia)
common <- Reduce(intersect, list(rownames(X), rownames(Z), colnames(Y_counts)))
common <- sort(common)

X <- X[common, , drop=FALSE]
Z <- Z[common, , drop=FALSE]
Y_counts <- Y_counts[, common, drop=FALSE]

stopifnot(
  identical(rownames(X), rownames(Z)),
  identical(rownames(X), colnames(Y_counts))
)

## -------------------------------
## 5) HTO demultiplex → donor_id
## -------------------------------
hto_sp <- as(hto_mat, "dgCMatrix")
stopifnot(!is.null(rownames(hto_sp)), !is.null(colnames(hto_sp)))

# 2) Create a Seurat object with HTO as the ONLY assay (no dummy counts needed)
shto <- CreateSeuratObject(counts = hto_sp, assay = "HTO", project = "HTO_only")

# 3) Normalize and demultiplex on the HTO assay
shto <- NormalizeData(shto, assay = "HTO", normalization.method = "CLR")
shto <- HTODemux(shto, assay = "HTO", positive.quantile = 0.97)  # tweak 0.95–0.98 if needed

table(shto$HTO_maxID)
table(shto$HTO_classification.global)

# 4) Extract classifications and winner hashtag per cell


# 5) Keep only Singlets (recommended for modeling)
keep_idx <- shto$HTO_classification.global == "Singlet"
keep_cells <- colnames(shto)[keep_idx]

groups <- factor(shto$HTO_maxID[keep_idx], exclude = NULL)
names(groups) <- keep_cells

length(groups)                 # how many cells?
length(unique(groups))         # how many groups? (you saw 1)
table(groups, useNA = "ifany") # distribution

# (optional) drop tiny groups that break folds
tab <- table(groups)
small <- names(tab)[tab < 20]  # tweak threshold
if (length(small)) {
  keep_cells <- keep_cells[!(groups %in% small)]
  groups <- droplevels(groups[!(groups %in% small)])
}


# 6) Align & subset ALL matrices by barcode (order by keep_cells!)
#    (Do NOT subset by a logical vector unless it’s in the same order)
keep_cells_aligned <- intersect(keep_cells, rownames(X))
length(keep_cells)          # 19208
length(keep_cells_aligned) 


X <- X[keep_cells_aligned, , drop=FALSE]
Z <- Z[keep_cells_aligned, , drop=FALSE]
Y_counts <- Y_counts[, keep_cells_aligned, drop=FALSE]
groups <- groups[keep_cells_aligned]   # align donor labels
groups <- droplevels(groups) 
donors <- levels(groups)
stopifnot(
  nrow(X) == nrow(Z),
  identical(rownames(X), rownames(Z)),
  identical(colnames(Y_counts), rownames(X)),
  length(groups) == nrow(X)
)

# ============================================================
# 6) Train/test split + save RDS for FMLE/scLinear baselines
# ============================================================
set.seed(42)
test_donors  <- sample(donors, size = 2)
train_donors <- setdiff(donors, test_donors)
train_cells <- names(groups)[groups %in% train_donors]
test_cells  <- names(groups)[groups %in% test_donors]

saveRDS(train_donors, file.path(out_base, "train_donors.rds"))
saveRDS(test_donors,  file.path(out_base, "test_donors.rds"))
saveRDS(train_cells, file.path(out_base, "train_cells.rds"))
saveRDS(test_cells,  file.path(out_base, "test_cells.rds"))
saveRDS(X,           file.path(out_base, "X.rds"))        # cells x genes
saveRDS(Z,           file.path(out_base, "Z.rds"))        # cells x PCs
saveRDS(Y_counts,    file.path(out_base, "adt_mat.rds"))  # proteins x cells
saveRDS(groups,      file.path(out_base, "groups.rds"))   # factor named by cell

# ============================================================
# 7) Export cTPnet CSVs (cells x features)
# ============================================================
rna_train <- X[train_cells, , drop=FALSE]
rna_test  <- X[test_cells,  , drop=FALSE]

adt_train <- t(Y_counts[, train_cells, drop=FALSE])  # cells x proteins
adt_test  <- t(Y_counts[, test_cells,  drop=FALSE])  # cells x proteins

stopifnot(identical(rownames(rna_train), rownames(adt_train)))
stopifnot(identical(rownames(rna_test),  rownames(adt_test)))

write.csv(rna_train, file.path(out_ctp, "rna_train.csv"))
write.csv(rna_test,  file.path(out_ctp, "rna_test.csv"))
write.csv(adt_train, file.path(out_ctp, "adt_train.csv"))
write.csv(adt_test,  file.path(out_ctp, "adt_test.csv"))

cat("DONE\n",
    "Final cells:", length(common), "\n",
    "Proteins:", nrow(Y_counts), "\n",
    "HVGs:", ncol(X), "\n",
    "Saved RDS:", out_base, "\n",
    "Saved cTPnet CSV:", out_ctp, "\n", sep="")



writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_root, "sessionInfo.txt"))
saveRDS(cfg, file.path(cfg$out_root, "config_used.rds"))


