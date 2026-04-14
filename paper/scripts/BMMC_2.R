
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
  library(Azimuth)
})


# ============================================================
# Bone-marrow Ab-seq dataset (97 proteins + RNA)
# Bone marrow Ab-seq multimodal dataset
# ============================================================
# https://figshare.com/articles/dataset/Expression_of_97_surface_markers_and_RNA_transcriptome_wide_in_13165_cells_from_a_healthy_young_bone_marrow_donor/13397987

source(file.path(here::here(), "paper", "scripts", "_config.R"))

out_base <- file.path(cfg$out_root, "citeseq_v1")
out_ctp  <- file.path(cfg$out_root, "ctp")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ctp, recursive = TRUE, showWarnings = FALSE)

seu <- readRDS(file.path(cfg$data_root, "WTA_projected.rds"))
rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
ab_counts  <- GetAssayData(seu, assay = "AB",  layer = "counts")
seu_clean <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "BM_AbSeq"
)
seu_clean[["ADT"]] <- CreateAssayObject(counts = ab_counts)
stopifnot(identical(colnames(seu_clean[["RNA"]]), colnames(seu_clean[["ADT"]])))
seu_clean
Assays(seu_clean)

seu_clean <- NormalizeData(seu_clean, assay="RNA", normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)

Layers(seu_clean[["RNA"]])

dim(GetAssayData(seu_clean, assay="RNA", layer="counts"))
dim(GetAssayData(seu_clean, assay="RNA", layer="data"))
ncol(seu_clean)

# =========================
# QC (RNA + ADT-total outliers) on merged object
# =========================

DefaultAssay(seu_clean) <- "RNA"
seu_clean$nCount_RNA   <- Matrix::colSums(GetAssayData(seu_clean, assay="RNA", layer="counts"))
seu_clean$nFeature_RNA <- Matrix::colSums(GetAssayData(seu_clean, assay="RNA", layer="counts") > 0)
seu_clean[["percent.mt"]] <- PercentageFeatureSet(seu_clean, pattern="^MT-")
VlnPlot(seu_clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
cat("QC summaries:\n")
print(summary(seu_clean$nFeature_RNA))
print(summary(seu_clean$nCount_RNA))
print(summary(seu_clean$percent.mt))
n_before <- ncol(seu_clean)
# choose cutoffs (match your PBMC style; adjust if needed)
f_hi <- as.numeric(quantile(seu_clean$nFeature_RNA, 0.995))
c_hi <- as.numeric(quantile(seu_clean$nCount_RNA,   0.995))

seu_clean <- subset(
  seu_clean,
  subset =
    nFeature_RNA >= 500 &
    nCount_RNA   >= 700 &
    nFeature_RNA <= f_hi &
    nCount_RNA   <= c_hi 
)
n_after <- ncol(seu_clean)
c(n_before = n_before, n_after = n_after, removed = n_before - n_after)
# ADT-total winsor filter (counts)
adt_counts0 <- GetAssayData(seu_clean, assay="ADT", layer="counts")
cs <- Matrix::colSums(adt_counts0)
cat("ADT total summary:\n")
print(summary(cs))

high <- as.numeric(quantile(cs, 0.995))
sum(cs > high)
low  <- as.numeric(quantile(cs, 0.001))
sum(cs < low)
quantile(cs, c(0.001, 0.005, 0.01, 0.02, 0.05))
seu_clean <- subset(seu_clean, cells = names(cs)[cs >= low & cs <= high])

# =========================
# preprocess
# =========================
set.seed(42)
seu_clean <- preprocess(
  object = seu_clean,
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
stopifnot("ADT" %in% names(seu_clean@assays))
stopifnot(identical(colnames(seu_clean[["RNA"]]), colnames(seu_clean[["ADT"]])))

seu_clean <- NormalizeData(seu_clean, assay="RNA", normalization.method="LogNormalize", scale.factor=1e4)
seu_clean <- FindVariableFeatures(seu_clean, assay="RNA", selection.method="vst", nfeatures=2000, verbose=FALSE)
hvg <- VariableFeatures(seu_clean)

seu_clean <- ScaleData(seu_clean, assay="RNA", features=hvg, verbose=FALSE)

seu_clean <- RunPCA(seu_clean, assay="RNA", features=hvg, npcs=20, verbose=FALSE)

seu_clean <- RunAzimuth(seu_clean, reference = "bonemarrowref")

Zfull <- Embeddings(seu_clean, "pca")
Z <- Zfull[, seq_len(min(20, ncol(Zfull))), drop=FALSE]

X <- t(as.matrix(GetAssayData(seu_clean, assay="RNA", layer="data")[hvg, , drop=FALSE]))
rownames(X) <- colnames(seu_clean)

Y_counts <- GetAssayData(seu_clean, assay="ADT", layer="counts")  # proteins x cells

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
saveRDS(seu_clean, file.path(out_base, "seu_final.rds"))
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






