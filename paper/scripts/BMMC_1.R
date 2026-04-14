
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
  use_condaenv("any", required = TRUE)
  py_config()
})


# ============================================================
# BMMC CITE-seq benchmark (GSE194122)
# ============================================================
# wget -c https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz
# --2026-02-25 16:32:19--  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE194nnn/GSE194122/suppl/GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz

source(file.path(here::here(), "paper", "scripts", "_config.R"))

out_base <- file.path(cfg$out_root, "citeseq_v1")
out_ctp  <- file.path(cfg$out_root, "ctp")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ctp, recursive = TRUE, showWarnings = FALSE)



anndata <- import("anndata", convert = FALSE)
scipy   <- import("scipy.sparse", convert = FALSE)

path <- file.path(
  cfg$data_root,
  "GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad"
)

adata <- anndata$read_h5ad(path)
ftype <- py_to_r(adata$var$feature_types)
table(ftype)
rna_idx <- which(ftype == "GEX")
adt_idx <- which(ftype == "ADT")

length(rna_idx)
length(adt_idx)
obs_df <- py_to_r(adata$obs)


counts_py <- adata$layers$get("counts")
counts_coo <- scipy$coo_matrix(counts_py)
i <- as.integer(py_to_r(counts_coo$row)) + 1
j <- as.integer(py_to_r(counts_coo$col)) + 1
x <- as.numeric(py_to_r(counts_coo$data))

dims <- dim(py_to_r(counts_py))

counts_clean <- sparseMatrix(
  i = i,
  j = j,
  x = x,
  dims = dims,
  giveCsparse = TRUE
)

counts_clean <- t(counts_clean)

var_names <- as.character(py_to_r(adata$var_names$to_list()))
obs_names <- as.character(py_to_r(adata$obs_names$to_list()))

dimnames(counts_clean) <- list(var_names, obs_names)

rna_counts <- counts_clean[rna_idx, ]
adt_counts <- counts_clean[adt_idx, ]
dim(rna_counts)
dim(adt_counts)
seu <- CreateSeuratObject(counts = rna_counts, assay="RNA")
seu[["ADT"]] <- CreateAssayObject(counts = adt_counts)



seu <- NormalizeData(seu, assay="RNA", normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)

Layers(seu[["RNA"]])

dim(GetAssayData(seu, assay="RNA", layer="counts"))
dim(GetAssayData(seu, assay="RNA", layer="data"))
ncol(seu)

stopifnot(identical(colnames(seu[["RNA"]]), colnames(seu[["ADT"]])))

# =========================
# QC (RNA + ADT-total outliers) on merged object
# =========================

DefaultAssay(seu) <- "RNA"
seu$nCount_RNA   <- Matrix::colSums(GetAssayData(seu, assay="RNA", layer="counts"))
seu$nFeature_RNA <- Matrix::colSums(GetAssayData(seu, assay="RNA", layer="counts") > 0)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")
cat("QC summaries:\n")
print(summary(seu$nFeature_RNA))
print(summary(seu$nCount_RNA))
print(summary(seu$percent.mt))

# choose cutoffs (match your PBMC style; adjust if needed)
f_hi <- as.numeric(quantile(seu$nFeature_RNA, 0.995))
c_hi <- as.numeric(quantile(seu$nCount_RNA,   0.995))

seu <- subset(
  seu,
  subset =
    nFeature_RNA >= 500 &
    nCount_RNA   >= 1500 &
    nFeature_RNA <= f_hi &
    nCount_RNA   <= c_hi &
    percent.mt   <= 15
)

# ADT-total winsor filter (counts)
adt_counts0 <- GetAssayData(seu, assay="ADT", layer="counts")
cs <- Matrix::colSums(adt_counts0)
cat("ADT total summary:\n")
print(summary(cs))
high <- as.numeric(quantile(cs, 0.995))
low <- 300
seu <- subset(seu, cells = names(cs)[cs >= low & cs <= high])

stopifnot(identical(colnames(seu[["RNA"]]), colnames(seu[["ADT"]])))

# =========================
# preprocess 
# =========================
set.seed(42)
seu <- preprocess(
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



# ensure ADT still clean and aligned
stopifnot("ADT" %in% names(seu@assays))
stopifnot(identical(colnames(seu[["RNA"]]), colnames(seu[["ADT"]])))

obs_df$cell_type <- as.character(py_to_r(adata$obs$cell_type))
seu$cell_type <- obs_df[colnames(seu), "cell_type"]

all(colnames(seu) %in% rownames(obs_df))

seu <- NormalizeData(seu, assay="RNA", normalization.method="LogNormalize", scale.factor=1e4)
seu <- FindVariableFeatures(seu, assay="RNA", selection.method="vst", nfeatures=2000, verbose=FALSE)
hvg <- VariableFeatures(seu)

seu <- ScaleData(seu, assay="RNA", features=hvg, verbose=FALSE)
ElbowPlot(seu, ndims = 30)
seu <- RunPCA(seu, assay="RNA", features=hvg, npcs=20, verbose=FALSE)

Zfull <- Embeddings(seu, "pca")
Z <- Zfull[, seq_len(min(20, ncol(Zfull))), drop=FALSE]

X <- t(as.matrix(GetAssayData(seu, assay="RNA", layer="data")[hvg, , drop=FALSE]))
rownames(X) <- colnames(seu)

Y_counts <- GetAssayData(seu, assay="ADT", layer="counts")  # proteins x cells

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
# Train/test split + save RDS for FMLE/scLinear baselines
# ============================================================
set.seed(42)
all_cells <- rownames(X)
train_cells <- sample(all_cells, size = floor(0.7 * length(all_cells)))
test_cells  <- setdiff(all_cells, train_cells)

saveRDS(train_cells, file.path(out_base, "train_cells.rds"))
saveRDS(test_cells,  file.path(out_base, "test_cells.rds"))
saveRDS(X,           file.path(out_base, "X.rds"))        # cells x genes
saveRDS(Z,           file.path(out_base, "Z.rds"))        # cells x PCs
saveRDS(Y_counts,    file.path(out_base, "adt_mat.rds"))  # proteins x cells
saveRDS(seu, file.path(out_base, "seu_final.rds"))
# ============================================================
# Export cTPnet CSVs (cells x features)
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




