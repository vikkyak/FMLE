

# PBMC 10k CITE-seq (10x Genomics, v3 chemistry)

# https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_protein_v3/pbmc_10k_protein_v3_filtered_feature_bc_matrix.tar.gz


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


source(file.path(here::here(), "paper", "scripts", "_config.R"))

out_base <- file.path(cfg$out_root, "citeseq_v1")
out_ctp  <- file.path(cfg$out_root, "ctp")
dir.create(out_base, recursive=TRUE, showWarnings=FALSE)
dir.create(out_ctp,  recursive=TRUE, showWarnings=FALSE)

data_dir <- file.path(cfg$data_root, "filtered_feature_bc_matrix/")
raw <- Read10X(data.dir = data_dir)
names(raw)

drop_proteins <- function(prot_names) {
  drop_idx <- grepl("^total$", prot_names, ignore.case=TRUE) |
    grepl("isotype|control|IgG", prot_names, ignore.case=TRUE)
  prot_names[!drop_idx]
}

seu <- CreateSeuratObject(
  counts  = raw[["Gene Expression"]],
  assay   = "RNA",
  project = "PBMC10k_CITE"
)

if ("Antibody Capture" %in% names(raw)) {
  adt <- as(raw[["Antibody Capture"]], "dgCMatrix")
} else {
  stopifnot(file.exists(path.expand(adt_csv_gz)))
  adt_df <- read.csv(gzfile(path.expand(adt_csv_gz)), row.names=1, check.names=FALSE)
  adt <- Matrix::Matrix(as.matrix(adt_df), sparse=TRUE)
}

common <- intersect(colnames(seu), colnames(adt))
cat("common cells:", length(common), "\n")
stopifnot(length(common) > 0)
seu <- subset(seu, cells = common)
adt <- adt[, colnames(seu), drop=FALSE]

# Add ADT assay
seu[["ADT"]] <- CreateAssayObject(counts = adt)
stopifnot(identical(colnames(seu[["RNA"]]), colnames(seu[["ADT"]])))

# Filter proteins once, consistently
keep_prot <- drop_proteins(rownames(seu[["ADT"]]))
seu[["ADT"]] <- subset(seu[["ADT"]], features = keep_prot)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")

VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

summary(seu$nFeature_RNA)
summary(seu$nCount_RNA)
summary(seu$percent.mt)

seu <- subset(
  seu,
  subset =
    nFeature_RNA > 500 &      # remove very poor cells
    nFeature_RNA < 4000 &     # cut off multiplets (6k–7k genes)
    nCount_RNA   < 30000 &    # drop extreme UMI outliers (73k)
    percent.mt   < 15         # remove high-mito damaged cells
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
# 3) preprocess
# ============================================================
seu_final <- preprocess(
  object = seu,
  remove_doublets = TRUE,
  low_qc_cell_removal = TRUE,
  remove_empty_droplets = FALSE,
  resolution = 0.8,
  do_clustering = FALSE,
  seed = 42,
  return_plots = FALSE,
  print_plots = FALSE,
  species = "Hs"
)

# Re-apply the same protein filter after preprocess (in case it modified assays)
adt_counts <- GetAssayData(seu_final, assay="ADT", layer="counts")
keep_prot <- drop_proteins(rownames(adt_counts))
adt_counts <- adt_counts[keep_prot, , drop=FALSE]
seu_final[["ADT"]] <- CreateAssayObject(counts = adt_counts)

stopifnot(identical(colnames(seu_final[["RNA"]]), colnames(seu_final[["ADT"]])))

saveRDS(seu_final, file.path(out_base, "seu_final.rds"))


# ============================================================
# 4) Build X (HVG lognorm), Z (PCA), Y (raw ADT counts) from FINAL object
# ============================================================
seu_final <- NormalizeData(seu_final, assay="RNA", normalization.method="LogNormalize", scale.factor=cfg$scale_factor)
seu_final <- FindVariableFeatures(seu_final, assay="RNA", selection.method="vst", nfeatures=cfg$hvg_n, verbose=FALSE)
hvg <- VariableFeatures(seu_final)

# ============================================================
# # Cross trasfer common genes from all three datasets are gene_panel_2000.rds
# bench <- 2
# base <- path.expand(sprintf("~/cfg$data_root/FMLE/benchmarks_%d", bench))
# scl    <- file.path(base, "sclinear") 
# ds      <- "citeseq_v1"
# seu_final <- readRDS(file.path(base, ds, "seu_final.rds"))  
# hvg <- readRDS("~/cfg$data_root/FMLE/transfer_preds/gene_panel_2000.rds")
# stopifnot(all(hvg %in% rownames(GetAssayData(seu_final, assay="RNA", layer="data"))))
# stopifnot(anyDuplicated(hvg) == 0)
# seu_final <- NormalizeData(seu_final, assay="RNA", normalization.method="LogNormalize", scale.factor=1e4)
# ============================================================


seu_final <- ScaleData(seu_final, assay="RNA", features=hvg, verbose=FALSE)
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

# ============================================================
# 5) Train/test split + save RDS for FMLE/scLinear baselines
# ============================================================

all_cells <- rownames(X)
train_cells <- sample(all_cells, size = floor(0.7 * length(all_cells)))
test_cells  <- setdiff(all_cells, train_cells)

saveRDS(train_cells, file.path(out_base, "train_cells.rds"))
saveRDS(test_cells,  file.path(out_base, "test_cells.rds"))
saveRDS(X,           file.path(out_base, "X.rds"))        # cells x genes
saveRDS(Z,           file.path(out_base, "Z.rds"))        # cells x PCs
saveRDS(Y_counts,    file.path(out_base, "adt_mat.rds"))  # proteins x cells

# ============================================================
# 6) Export cTPnet CSVs (cells x features)
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





