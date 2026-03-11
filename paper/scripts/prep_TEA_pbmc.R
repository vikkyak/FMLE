#!/usr/bin/env Rscript

# ============================================================
# TEA-seq PBMC (GSE158013 / GSM4949911): RNA+ADT benchmark export
# Produces:
#   1) seu_final.rds                 (final Seurat object after preprocess_fmle)
#   2) citeseq_v1/{X.rds,Z.rds,adt_mat.rds,train_cells.rds,test_cells.rds}
#   3) ctp/{rna_train.csv,rna_test.csv,adt_train.csv,adt_test.csv}  (cells x features)
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(reticulate)
  library(scLinear)
})

files <- c("_config.R", "preprocess_helpers_fmle.R", "preprocess_fmle.R")
paths <- file.path(here::here(), "paper", "scripts", files)
stopifnot(all(file.exists(paths)))
invisible(lapply(paths, source))
# -----------------------
# Paths
# -----------------------
tea_dir <- file.path(cfg$data_root, "tea_seq")
h5 <- file.path(
  tea_dir,
  "GSM4949911_X061-AP0C1W1_leukopak_perm-cells_tea_fulldepth_cellranger-arc_filtered_feature_bc_matrix.h5"
)

adt_csv_gz <- file.path(
  tea_dir,
  "GSM4949911_adt_counts.csv.gz"
)   # only used if H5 lacks Antibody Capture


out_base <- file.path(cfg$out_root, "citeseq_v1")
out_ctp  <- file.path(cfg$out_root, "ctp")
dir.create(out_base, recursive=TRUE, showWarnings=FALSE)
dir.create(out_ctp,  recursive=TRUE, showWarnings=FALSE)

# -----------------------
# Helpers
# -----------------------
drop_proteins <- function(prot_names) {
  drop_idx <- grepl("^total$", prot_names, ignore.case=TRUE) |
    grepl("isotype|control|IgG", prot_names, ignore.case=TRUE)
  prot_names[!drop_idx]
}

# ============================================================
# 1) Load RNA+ADT (CORRECT, ROBUST)
#    Final goal: ADT = proteins x cells, colnames(ADT) == colnames(seu)
# ============================================================
stopifnot(file.exists(path.expand(h5)))
multi <- Read10X_h5(path.expand(h5))
stopifnot("Gene Expression" %in% names(multi))

rna <- multi[["Gene Expression"]]
seu <- CreateSeuratObject(counts = rna)

# ---- load ADT ----
if ("Antibody Capture" %in% names(multi)) {
  adt <- as(multi[["Antibody Capture"]], "dgCMatrix")
} else {
  stopifnot(file.exists(path.expand(adt_csv_gz)))
  adt_df <- read.csv(gzfile(path.expand(adt_csv_gz)), row.names=1, check.names=FALSE)
  adt <- Matrix::Matrix(as.matrix(adt_df), sparse=TRUE)
}

# ---- ensure proteins x cells (only transpose if proteins are in COLNAMES) ----
# Your debug showed colnames(adt) were "CD10-1", "CD3-1", ... so transpose ONCE.
if (any(grepl("^(total|CD|HLA|Ig|TCR|FceRI|KLRG1)", colnames(adt)))) {
  adt <- t(adt)
}

# ---- harmonize "-1" suffix to match RNA barcodes ----
rna_cells <- colnames(seu)
if (mean(grepl("-1$", rna_cells)) > 0.5 && mean(grepl("-1$", colnames(adt))) < 0.5) {
  colnames(adt) <- paste0(colnames(adt), "-1")
} else if (mean(grepl("-1$", rna_cells)) < 0.5 && mean(grepl("-1$", colnames(adt))) > 0.5) {
  colnames(adt) <- sub("-1$", "", colnames(adt))
}

# ---- overlap + align ----
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

# ============================================================
# 2) QC on RNA + optional ADT-total outlier filter
# ============================================================
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")

qc_ok <- seu$nFeature_RNA >= 300 &
  seu$nFeature_RNA <= 6000 &
  seu$percent.mt <= 20

seu <- subset(seu, cells = colnames(seu)[qc_ok])

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
# ensure reticulate env is set (adjust name if needed)
seu_final <- preprocess_fmle(
  object = seu,
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

# Re-apply the same protein filter after preprocess (in case it modified assays)
adt_counts <- GetAssayData(seu_final, assay="ADT", layer="counts")
keep_prot <- drop_proteins(rownames(adt_counts))
adt_counts <- adt_counts[keep_prot, , drop=FALSE]
seu_final[["ADT"]] <- CreateAssayObject(counts = adt_counts)

stopifnot(identical(colnames(seu_final[["RNA"]]), colnames(seu_final[["ADT"]])))

saveRDS(seu_final, file.path(out_root, "seu_final.rds"))

# ============================================================
# 4) Build X (HVG lognorm), Z (PCA), Y (raw ADT counts) from FINAL object
# ============================================================
seu_final <- NormalizeData(seu_final, assay="RNA", normalization.method="LogNormalize", scale.factor=cfg$scale_factor)
seu_final <- FindVariableFeatures(seu_final, assay="RNA", selection.method="vst", nfeatures=cfg$hvg_n, verbose=FALSE)
hvg <- VariableFeatures(seu_final)

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




writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_root, "sessionInfo.txt"))
saveRDS(cfg, file.path(cfg$out_root, "config_used.rds"))






