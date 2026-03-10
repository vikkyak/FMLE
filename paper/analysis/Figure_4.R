
# lymph_node

# wget -L -O V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5 \
# https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5
# 
# # spatial assets (images + coordinates)
# wget -L -O V1_Human_Lymph_Node_spatial.tar.gz \
# https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz
# 
# tar -xvzf V1_Human_Lymph_Node_spatial.tar.gz


# tonsils 
# https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard

# 1) Feature / barcode matrix (filtered)

# 2) Spatial imaging data


library(Seurat)
library(Matrix)
library(ggplot2)
library(data.table)



# 1) tonsil
mat <- Read10X(data.dir = "/home/vikas/Desktop/FMLE/spatial/tonsil/filtered_feature_bc_matrix")
names(mat)
seu_tonsil <- CreateSeuratObject(counts = mat[["Gene Expression"]], assay = "Spatial")
seu_tonsil[["ADT"]] <- CreateAssayObject(counts = mat[["Antibody Capture"]])

Assays(seu_tonsil)
dim(seu_tonsil[["ADT"]])

img <- Read10X_Image("/home/vikas/Desktop/FMLE/spatial/tonsil/spatial", filter.matrix = TRUE)
img <- img[colnames(seu_tonsil)]   # subset image to EXACT same barcodes
seu_tonsil[["slice1"]] <- img


pos <- fread("/home/vikas/Desktop/FMLE/spatial/tonsil/spatial/tissue_positions.csv")

# keep only spots in object + align order
pos <- pos[barcode %in% colnames(seu_tonsil)]
setkey(pos, barcode)
pos <- pos[colnames(seu_tonsil)]

# store coords (fullres pixel coords)
seu_tonsil$imagerow <- pos$pxl_row_in_fullres
seu_tonsil$imagecol <- pos$pxl_col_in_fullres
seu_tonsil$in_tissue <- pos$in_tissue

head(seu_tonsil@meta.data[, c("imagerow","imagecol","in_tissue")])

# 2) lymph
# seu_ln <- Load10X_Spatial(
#   data.dir = "/home/vikas/Desktop/FMLE/spatial/lymph_node",
#   filename = "V1_Human_Lymph_Node_filtered_feature_bc_matrix.h5"
# )

# 1) normalize Visium same as PBMC scale
# seu_analy <- seu_ln
seu_analy <- seu_tonsil
DefaultAssay(seu_analy) <- "Spatial"

seu_analy <- NormalizeData(seu_analy, assay="Spatial", normalization.method="LogNormalize",
  scale.factor=1e4, verbose=FALSE
)

# seu_analy <- NormalizeData(seu_analy, assay="ADT", normalization.method="CLR", margin=2, verbose=FALSE)
seu_analy[["percent.mt"]] <- PercentageFeatureSet(seu_analy, pattern="^MT-")
summary(seu_analy$percent.mt)
summary(seu_analy$nFeature_Spatial)
table(seu_analy$in_tissue)
seu_analy <- subset(seu_analy, subset = in_tissue == 1)

# Optional: remove extreme ADT-total outliers (winsor at 99.5%)
# adt_counts0 <- GetAssayData(seu_analy, assay="ADT", layer="counts")
# cs <- Matrix::colSums(adt_counts0)
# summary(cs)
# cor(log1p(cs), log1p(seu_analy$nCount_Spatial))
# cor(log1p(cs), seu_analy$nFeature_Spatial)
# thr <- as.numeric(quantile(cs, 0.995))
# keep_cells <- names(cs)[cs <= thr]
bench <- 1  # set 1 or 2 or 3
base1  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_spatial_%d", bench))
ds    <- "citeseq_v1"
dir.create(base1, recursive = TRUE, showWarnings = FALSE)

train_cells <- readRDS(file.path(base1, ds, "train_cells.rds"))
test_cells  <- readRDS(file.path(base1, ds, "test_cells.rds"))
X           <- readRDS(file.path(base1, ds, "X.rds"))
Z           <- readRDS(file.path(base1, ds, "Z.rds"))
adt_mat     <- readRDS(file.path(base1, ds, "adt_mat.rds"))
seu_prep    <- readRDS(file.path(base1, ds, "seu_final.rds"))
if (bench == 1) {
  groups    <- readRDS(file.path(base1, ds, "groups.rds"))
} else {
  groups <- NULL
}

X_train <- X[train_cells, , drop=FALSE]
Z_train <- Z[train_cells, , drop=FALSE]
X_test  <- X[test_cells,  , drop=FALSE]
Z_test  <- Z[test_cells,  , drop=FALSE]

stopifnot(identical(rownames(X_train), rownames(Z_train)))
stopifnot(identical(rownames(X_test),  rownames(Z_test)))

# Use the SAME scaled X for both FMLE and PCA projection
genes_train <- colnames(X_train)  

mu  <- Matrix::colMeans(as(X_train, "dgCMatrix"))
mu2 <- Matrix::colMeans(as(X_train, "dgCMatrix")^2)
sd  <- sqrt(pmax(mu2 - mu^2, 0))
sd[!is.finite(sd) | sd == 0] <- 1

# Build spatial X in the same gene space/order
X_sp_log <- t(as(GetAssayData(seu_analy, assay="Spatial", layer="data"), "dgCMatrix"))
stopifnot(all(genes_train %in% colnames(X_sp_log)))

X_sp <- X_sp_log[, genes_train, drop=FALSE]
stopifnot(identical(colnames(X_sp), genes_train))

# Scale spatial with PBMC train stats
X_sp_scaled <- sweep(X_sp, 2, mu, "-")
X_sp_scaled <- sweep(X_sp_scaled, 2, sd, "/")

pcs <- colnames(Z_train)                 # e.g. PC_1..PC_k
# L <- Loadings(seu_prep, "pca")[genes_train, pcs, drop=FALSE] # genes_in_pca x pcs
L <- readRDS(file.path(base1, "citeseq_v1", "pca_loadings.rds"))
stopifnot(ncol(L) == length(pcs))
stopifnot(
  identical(colnames(X_sp_scaled), names(mu)),
  identical(colnames(X_sp_scaled), rownames(L)),
  identical(rownames(L), genes_train)
)

Z_sp <- as.matrix(X_sp_scaled %*% as.matrix(L))
# check 
# You should see non-zero SDs roughly similar across PCs.
# If a PC SD ≈ 0 → gene mismatch or scaling error.
apply(Z_sp, 2, sd)

colnames(Z_sp) <- pcs


library(glue)

rds_dir   <- "~/Desktop/FMLE/benchmarks_spatial_1/FMLE/citeseq_v1_final"        # where final_*.rds live
sp_outdir <- "~/Desktop/FMLE/benchmarks_spatial_1/FMLE/spatial_preds"
dir.create(sp_outdir, showWarnings = FALSE, recursive = TRUE)

coords <- GetTissueCoordinates(seu_analy)  # x,y + barcode

stopifnot(identical(rownames(X_sp_scaled), rownames(coords)))
# list of proteins you have
rds_files <- list.files(rds_dir, pattern = "^final_.*\\.rds$", full.names = TRUE)

map_pbmc_to_tonsil <- c(
  # exact same antigen (name differs by suffix/case)
  "CD14"        = "CD14.1",
  "CD16"        = "FCGR3A.1",
  "CD19"        = "CD19.1",
  "CD4"         = "CD4.1",
  "CD8a"        = "CD8A.1",
  "CD27"        = "CD27.1",
  
  # PBMC panel name vs tonsil gene-specific/proxy
  "CD197-CCR7"  = "CCR7.1",
  "HLA.DR"      = "HLA-DRA",
  "CD3"         = "CD3E.1"
)

keep_ln     <- c("CD11a","CD11c","CD123","CD127-IL7Ra","CD14","CD16","CD161","CD19",
                 "CD197-CCR7","CD25","CD27","CD278-ICOS","CD28","CD3","CD34","CD38",
                 "CD4","CD45RA","CD45RO","CD56","CD57","CD69","CD79b","CD8a","HLA.DR")

dataset <- "tonsil"   # or "lymph_node"
keep_tonsil <- names(map_pbmc_to_tonsil) 

pred_maps <- list()

for (f in rds_files) {
  obj <- readRDS(f)
  prot <- obj$protein
  if (dataset == "tonsil") {
    if (!(prot %in% keep_tonsil)) next
  } else if (dataset == "lymph_node") {
    if (!(prot %in% keep_ln)) next
  }
  fit  <- obj$fit 
  
  pred <- fmle_predict(fit, X_new = X_sp_scaled, Z_new = Z_sp, return_se = TRUE)
  
  y_hat   <- as.numeric(pred$mean)         # standardized score per spot
  alpha   <- pred$alpha                    # spots × R
  entropy <- -rowSums(alpha * log(alpha + 1e-12))
  entropy <- entropy / log(ncol(alpha))
  
  df <- data.frame(
    spot = rownames(X_sp),
    protein = prot,
    y_hat = y_hat,
    entropy = entropy
  )
  
  # optional: dominant expert per spot
  dom_k <- max.col(alpha, ties.method = "first")
  df$dom_expert <- dom_k
  
  # attach coordinates
  df <- merge(df, cbind(spot = rownames(coords), coords), by="spot", all.x=TRUE, sort=FALSE)
  df <- df[match(rownames(X_sp_scaled), df$spot), ]
  stopifnot(
    nrow(df) == nrow(X_sp_scaled),
    identical(df$spot, rownames(X_sp_scaled)),
    all(!is.na(df$x)) & all(!is.na(df$y))    # tissue coords present
  )
  pred_maps[[prot]] <- list(df = df, alpha = alpha)
  
  saveRDS(pred_maps[[prot]], file.path(sp_outdir, glue("visium_pred_{prot}.rds")))
}

obj <- readRDS(file.path(rds_dir, "final_CD3.rds"))
tf  <- obj$tf_y
fit <- obj$fit

pred <- fmle_predict(fit, X_new = X_sp_scaled, Z_new = Z_sp, return_se = TRUE)
alpha_sp <- pred$alpha
a1 <- alpha_sp[,1]
a2 <- alpha_sp[,2]
margin <- apply(alpha_sp, 1, function(v) sort(v, decreasing=TRUE)[1] - sort(v, decreasing=TRUE)[2])

seu_analy$alpha1 <- a1
seu_analy$alpha2 <- a2
seu_analy$alpha_margin <- margin

R <- ncol(alpha_sp)
entropy <- -rowSums(alpha_sp * log(alpha_sp + 1e-12)) / log(R)

stopifnot(identical(rownames(seu_analy@meta.data), rownames(alpha_sp)))
seu_analy$entropy <- entropy

SpatialFeaturePlot(seu_analy, features=c("alpha1","alpha2","alpha_margin"))

SpatialFeaturePlot(seu_analy, c("entropy","alpha_margin"))

cor(seu_analy$alpha_margin, log1p(seu_analy$nCount_ADT), method="spearman")



plot_pred_vs_adt_tf <- function(seu, adt_name, pred_vec, tf, title_tag="") {
  
  # ---- DIAG 1: ADT feature name ----
  if (!(adt_name %in% rownames(seu[["ADT"]]))) {
    cat("\n[DIAG] ADT name not found:", adt_name, "\n")
    cat("[DIAG] First ADT features:\n")
    print(head(rownames(seu[["ADT"]]), 30))
    stop("Fix adt_name.")
  }
  stopifnot(identical(names(pred_vec), colnames(seu)) || is.null(names(pred_vec)))
  
  # raw ADT (in Seurat spot order)
  y_raw <- as.numeric(GetAssayData(seu, assay="ADT", layer="counts")[adt_name, ])
  
  # predicted (you provide it)
  x <- as.numeric(pred_vec)
  
  # ---- DIAG 2: length/order ----
  cat("\n[DIAG] ", title_tag, "\n", sep="")
  cat("[DIAG] nspots=", ncol(seu),
      " len(pred)=", length(x),
      " len(adt)=", length(y_raw), "\n", sep="")
  stopifnot(length(x) == length(y_raw))
  
  # ---- DIAG 3: tf transform ----
  y_tf <- cap_and_scale_apply(y_raw, tf)
  cat("[DIAG] pred summary:\n"); print(summary(x))
  cat("[DIAG] adt_raw summary:\n"); print(summary(y_raw))
  cat("[DIAG] adt_tf summary:\n"); print(summary(y_tf))
  
  ok <- is.finite(x) & is.finite(y_tf)
  cat("[DIAG] finite pairs:", sum(ok), "/", length(ok), "\n", sep="")
  
  rP <- suppressWarnings(cor(x[ok], y_tf[ok], method="pearson"))
  rS <- suppressWarnings(cor(x[ok], y_tf[ok], method="spearman"))
  cat("[DIAG] Pearson(tf)=", rP, " Spearman(tf)=", rS, "\n", sep="")
  
  ggplot(data.frame(pred=x[ok], adt_tf=y_tf[ok]), aes(pred, adt_tf)) +
    geom_point(size=0.5, alpha=0.35) +
    labs(
      title=paste0(title_tag, "  Pearson(tf)=", sprintf("%.3f", rP),
                   "  Spearman(tf)=", sprintf("%.3f", rS)),
      x="FMLE predicted (tf-scale)",
      y=paste0("Measured ADT (tf-scale): ", adt_name)
    ) +
    theme_classic()
}

obj <- readRDS(file.path(rds_dir, "final_CD3.rds"))
tf  <- obj$tf_y

p_cd3 <- plot_pred_vs_adt_tf(
  seu = seu_analy,
  adt_name = map_pbmc_to_tonsil["CD3"],     # "CD3E.1"
  pred_vec = seu_analy[["FMLE_CD3"]][,1],
  tf = tf,
  title_tag = "CD3(PBMC model) vs CD3E.1(Tonsil)"
)

