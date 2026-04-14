
# This fuction is used for validation of within dataset and cross donor dataset. 
# The base folder is changed as per datatset folder

suppressPackageStartupMessages({
  library(FMLE)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(Matrix)
  library(glmnet)
  library(glue)
  library(reticulate)
  use_condaenv("sclinear", required = TRUE)
  py_config()  
  library(scLinear)
})

#==================================
# 0) Load fixed benchmark artifacts
#==================================
bench <- 1  # set 1, 2, or 3  in case of PBMC and BMMC
ds <- "citeseq_v1"

base <- file.path(cfg$out_root, sprintf("benchmarks_%d", bench))
ds_dir <- file.path(base, ds)
cv_dir <- file.path(base, "FMLE", paste0(ds, "_cv"))

dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cv_dir, recursive = TRUE, showWarnings = FALSE)


train_cells <- readRDS(file.path(ds_dir, "train_cells.rds"))
test_cells  <- readRDS(file.path(ds_dir, "test_cells.rds"))
X      <- readRDS(file.path(ds_dir, "X.rds"))        # cells x genes
Z      <- readRDS(file.path(ds_dir, "Z.rds"))        # cells x PCs
adt_mat<- readRDS(file.path(ds_dir, "adt_mat.rds"))  # proteins x cells
seu_prep <- readRDS(file.path(base, ds, "seu_final.rds"))  
hvg <- readRDS(file.path(base, "hvg_genes.rds"))
stopifnot(all(train_cells %in% colnames(seu_prep)), all(test_cells %in% colnames(seu_prep)))

if (bench == 1) {
  groups <- readRDS(file.path(ds_dir, "groups.rds"))
} else {
  groups <- NULL
}

#===========================
# 1) Safety checks
#===========================
stopifnot(length(intersect(train_cells, test_cells)) == 0)
stopifnot(identical(rownames(X), rownames(Z)))
stopifnot(all(train_cells %in% rownames(X)), all(test_cells %in% rownames(X)))
stopifnot(all(train_cells %in% colnames(adt_mat)), all(test_cells %in% colnames(adt_mat)))

# consistent ordering everywhere
X_train <- X[train_cells, , drop=FALSE]
Z_train <- Z[train_cells, , drop=FALSE]
X_test  <- X[test_cells,  , drop=FALSE]
Z_test  <- Z[test_cells,  , drop=FALSE]

stopifnot(identical(rownames(X_train), rownames(Z_train)))
stopifnot(identical(rownames(X_test),  rownames(Z_test)))

#===========================
# 2) Helper functions
#===========================
calc_metrics <- function(yt, yp) {
  ok <- is.finite(yt) & is.finite(yp)
  yt <- yt[ok]; yp <- yp[ok]
  if (length(yt) < 20) {
    return(data.frame(MSE=NA, RMSE=NA, MAE=NA, Pearson=NA, Spearman=NA, R2=NA))
  }
  mse  <- mean((yt - yp)^2)
  rmse <- sqrt(mse)
  mae  <- mean(abs(yt - yp))
  pear <- suppressWarnings(cor(yt, yp, use="pairwise.complete.obs"))
  spear<- suppressWarnings(cor(yt, yp, method="spearman", use="pairwise.complete.obs"))
  den  <- sum((yt - mean(yt))^2)
  r2   <- if (den > 0) 1 - sum((yt - yp)^2)/den else NA_real_
  data.frame(MSE=mse, RMSE=rmse, MAE=mae, Pearson=pear, Spearman=spear, R2=r2)
}


extract_beta <- function(fit, feature_names = NULL) {
  R <- fit$R
  beta0 <- vapply(fit$experts, `[[`, numeric(1), "beta0")
  beta  <- do.call(cbind, lapply(fit$experts, `[[`, "beta"))  # p × R
  colnames(beta) <- paste0("expert", seq_len(R))
  if (!is.null(feature_names)) rownames(beta) <- feature_names
  list(beta0 = beta0, beta = beta)
}

topk_per_expert <- function(beta_mat, k = 50) {
  # returns a list of character vectors of top |beta|
  lapply(seq_len(ncol(beta_mat)), function(r) {
    b <- beta_mat[, r]
    ord <- order(abs(b), decreasing = TRUE)
    names(b)[ord][seq_len(min(k, length(ord)))]
  })
}

#=================================================================================#
# FMLE validation
#=================================================================================#
# Best hyperparams from CV (protein, R, m, lambda, mse_cv)
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cv_best_all <- readr::read_csv(file.path(cv_dir, "cv_best_all_adts.csv"), show_col_types = FALSE)

results_list <- vector("list", nrow(cv_best_all))
names(results_list) <- cv_best_all$protein

results_ct_list <- vector("list", nrow(cv_best_all))
names(results_ct_list) <- cv_best_all$protein

meta <- seu_prep@meta.data

for (i in seq_len(nrow(cv_best_all))) {
  prot <- cv_best_all$protein[i]
  best <- cv_best_all[i, ]
  
  cat("\n[FMLE FINAL] ", prot, " | R=", best$R, " m=", best$m, " lam=", best$lambda, "\n", sep="")
  
  # ---- y: train-only cap+log1p+z (FMLE internal scale) ----
  y_train_raw <- as.numeric(adt_mat[prot, train_cells])
  tf          <- FMLE:::cap_and_scale_fit(y_train_raw, q = 0.995)
  y_train_tf  <- FMLE:::cap_and_scale_apply(y_train_raw, tf)
  
  y_test_raw  <- as.numeric(adt_mat[prot, test_cells])
  y_test_tf   <- FMLE:::cap_and_scale_apply(y_test_raw, tf)
  
  # ---- fit final model on TRAIN ----
  fit <- fmle_train(
    X = X_train,
    y = y_train_tf,
    Z = Z_train,
    R           = as.integer(best$R),
    m           = as.numeric(best$m),
    lambda_l1   = as.numeric(best$lambda),
    ridge       = 1e-6,
    standardize = TRUE,
    seed        = 1,
    verbose     = FALSE
  )
  
  # ---- predict TEST ----
  pr <- fmle_predict(
    fit,
    X_new = X_test,
    Z_new = Z_test,
    return_se = TRUE
  )
  
  y_hat_test <- as.numeric(pr$mean)
  
  # ---- metrics on FMLE eval scale ----
  met <- calc_metrics(y_test_tf, y_hat_test)
  
  # ---- extract betas for interpretation ----
  bb  <- extract_beta(fit, feature_names = colnames(X_train))
  top <- topk_per_expert(bb$beta, k = 50)
  
  # ---- overall result row ----
  res_row <- data.frame(
    protein = prot,
    R       = as.integer(best$R),
    m       = as.numeric(best$m),
    lambda  = as.numeric(best$lambda),
    mse_cv  = as.numeric(best$mse_cv),
    met
  )
  results_list[[prot]] <- res_row
  
  # =========================================================
  # within-cell-type Pearson on TEST
  # =========================================================
  keep <- intersect(test_cells, rownames(meta))
  
  if (length(keep) > 0) {
    idx <- match(keep, test_cells)
    
    df_ct <- data.frame(
      protein   = prot,
      cell      = keep,
      cell_type = meta[keep, "cell_type", drop = TRUE],
      yt        = y_test_tf[idx],
      yp        = y_hat_test[idx],
      stringsAsFactors = FALSE
    )
    
    ct_metrics <- df_ct %>%
      dplyr::group_by(protein, cell_type) %>%
      dplyr::summarise(
        Pearson_fmle = {
          nn <- dplyr::n()
          if (nn < 2) NA_real_ else suppressWarnings(cor(yt, yp, use = "pairwise.complete.obs"))
        },
        n = dplyr::n(),
        .groups = "drop"
      )
  } else {
    ct_metrics <- NULL
  }
  
  results_ct_list[[prot]] <- ct_metrics
  
  saveRDS(
    list(
      protein = prot,
      best = best,
      tf_y = tf,
      fit = fit,
      y_test = y_test_tf,
      yhat_test = y_hat_test,
      se_test = if (!is.null(pr$se)) as.numeric(pr$se) else NULL,
      alpha_test = if (!is.null(pr$alpha)) as.matrix(pr$alpha) else NULL,
      beta0 = bb$beta0,
      beta  = bb$beta,
      top_genes = top,
      train_cells = train_cells,
      test_cells = test_cells,
      ct_metrics = ct_metrics
    ),
    file.path(out_dir, paste0("final_", prot, ".rds"))
  )
}

# ---- bind and save overall test metrics ----
results_fmle <- dplyr::bind_rows(results_list)
readr::write_csv(results_fmle, file.path(out_dir, "fmle_test_metrics.csv"))

# ---- bind and save within-cell-type Pearson ----
results_fmle_ct <- dplyr::bind_rows(results_ct_list)
readr::write_csv(results_fmle_ct, file.path(out_dir, "fmle_within_celltype_pearson.csv"))

#=================================================================================#
# scLinear validation
#=================================================================================#

scl    <- file.path(base, "sclinear") 
seu_prep_hvg <- seu_prep
seu_prep_hvg[["RNA"]] <- subset(seu_prep[["RNA"]], features = hvg)
Assays(seu_prep_hvg)
nrow(seu_prep_hvg[["RNA"]])
nrow(seu_prep_hvg[["ADT"]])
seu_train <- seu_prep_hvg[, train_cells]
seu_test  <- seu_prep_hvg[, test_cells]

# ---- train scLinear predictor on TRAIN ----
pipe <- scLinear::create_adt_predictor()

gexp_train <- Seurat::GetAssayData(seu_train, assay="RNA", layer="counts")  # genes x train
adt_train  <- Seurat::GetAssayData(seu_train, assay="ADT", layer="counts")  # proteins x train

pipe <- scLinear::fit_predictor(
  pipe          = pipe,
  gexp_train    = as.matrix(gexp_train),
  adt_train     = as.matrix(adt_train),   # raw ADT counts
  normalize_gex = TRUE,
  normalize_adt = FALSE                  
)
# ---- predict TRAIN ----
pred_train <- scLinear::adt_predict(pipe, gexp_train, layer ="counts", normalize=TRUE)  # proteins x test
pred_train <- as.matrix(pred_train@data)

# ---- predict TEST ----
gexp_test <- Seurat::GetAssayData(seu_test, assay="RNA", layer="counts")
pred_test <- scLinear::adt_predict(pipe, gexp_test, layer ="counts", normalize=TRUE)  # proteins x test
pred_test <- as.matrix(pred_test@data)

# ---- align cells & proteins ----
train_cells_eval <- intersect(train_cells, colnames(pred_train))
test_cells_eval  <- intersect(test_cells,  colnames(pred_test))
stopifnot(length(train_cells_eval) > 0, length(test_cells_eval) > 0)

pred_train <- pred_train[, train_cells_eval, drop=FALSE]
pred_test  <- pred_test[,  test_cells_eval,  drop=FALSE]

prot_use <- intersect(rownames(pred_test), rownames(adt_mat))
stopifnot(length(prot_use) > 0)
#----------------------------------------------------------------------------------#
## in case of PBMC10K
#----------------------------------------------------------------------------------#
# adt_controls <- grep("IgG", rownames(adt_mat), value = TRUE)
# adt_use      <- setdiff(rownames(adt_mat), adt_controls)
# prot_use <- adt_use  
#----------------------------------------------------------------------------------#
res_sclinear <- lapply(prot_use, function(p){
  
  y_train_raw <- as.numeric(adt_mat[p, train_cells_eval])
  tf <- FMLE:::cap_and_scale_fit(y_train_raw, q=0.995)
  y_train_tf <- FMLE:::cap_and_scale_apply(y_train_raw, tf)
  
  y_test_raw <- as.numeric(adt_mat[p, test_cells_eval])
  y_true <- FMLE:::cap_and_scale_apply(y_test_raw, tf)
  
  yhat_train_raw <- as.numeric(pred_train[p, train_cells_eval])
  yhat_test_raw  <- as.numeric(pred_test[p,  test_cells_eval])
  
  ok <- is.finite(y_train_tf) & is.finite(yhat_train_raw)
  if (sum(ok) < 20 || sd(yhat_train_raw[ok]) < 1e-8) {
    return(list(
      overall = data.frame(
        protein=p, MSE=NA, RMSE=NA, MAE=NA,
        Pearson=NA, Spearman=NA, R2=NA
      ),
      within_ct = data.frame(
        protein=p, cell_type=NA, Pearson=NA, n=NA
      )
    ))
  }
  
  fit_cal <- lm(y_train_tf[ok] ~ yhat_train_raw[ok])
  y_hat <- as.numeric(coef(fit_cal)[1] + coef(fit_cal)[2] * yhat_test_raw)
  
  ok2 <- is.finite(y_true) & is.finite(y_hat)
  yt <- y_true[ok2]; yp <- y_hat[ok2]
  
  mse <- mean((yt-yp)^2)
  den <- sum((yt-mean(yt))^2)
  
  overall_row <- data.frame(
    protein=p,
    MSE=mse,
    RMSE=sqrt(mse),
    MAE=mean(abs(yt-yp)),
    Pearson=suppressWarnings(cor(yt, yp)),
    Spearman=suppressWarnings(cor(yt, yp, method="spearman")),
    R2=if (den>0) 1 - sum((yt-yp)^2)/den else NA_real_
  )
  ct_all <- meta[test_cells_eval, "cell_type", drop=TRUE]
  ct_ok2 <- ct_all[ok2]
  
  within_ct <- data.frame(
    protein   = p,
    cell_type = as.character(ct_ok2),
    yt = yt,
    yp = yp,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(protein, cell_type) %>%
    dplyr::summarise(
      Pearson = suppressWarnings(cor(yt, yp)),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  # return both; overall stays identical to your current output
  list(overall = overall_row, within_ct = within_ct)
  
})

res_sclinear_overall <- dplyr::bind_rows(lapply(res_sclinear, `[[`, "overall"))
res_sclinear_ct      <- dplyr::bind_rows(lapply(res_sclinear, `[[`, "within_ct"))

#===============================================================================
# cross donor
#===============================================================================
# res_sclinear <- dplyr::bind_rows(res_sclinear)
# res_sclinear_overall <- res_sclinear$overall %>%
#   distinct(protein, .keep_all = TRUE)
# res_sclinear_overall <- res_sclinear_overall %>%
#   filter(
#     is.finite(Pearson),
#     is.finite(Spearman),
#     is.finite(R2)
#   )
# res_sclinear_ct <- res_sclinear$within_ct
# res_sclinear_ct <- res_sclinear_ct %>%
#   filter(is.finite(Pearson))
#===============================================================================


# write files
readr::write_csv(res_sclinear_overall, file.path(scl, "scLinear_test_metrics.csv"))
readr::write_csv(res_sclinear_ct,      file.path(scl, "scLinear_within_celltype_pearson.csv"))


#=================================================================================#
# ctpnet validation
#=================================================================================#

rm(pred_train, pred_test, train_cells_eval, test_cells_eval, prot_use )

ctp    <- file.path(base, "ctp") 

pred_train <- as.matrix(read.csv(file.path(ctp, "kaggle_ctpnet_pred_train.csv"),
                                 row.names=1, check.names=FALSE))
pred_test  <- as.matrix(read.csv(file.path(ctp, "kaggle_ctpnet_pred_test.csv"),
                                 row.names=1, check.names=FALSE))

stopifnot(length(intersect(train_cells, test_cells)) == 0)

train_cells_eval <- intersect(train_cells, colnames(pred_train))
test_cells_eval  <- intersect(test_cells,  colnames(pred_test))

cat("train_cells:", length(train_cells), "\n")
cat("train_cells_eval:", length(train_cells_eval), "\n")
cat("test_cells:", length(test_cells), "\n")
cat("test_cells_eval:", length(test_cells_eval), "\n")


stopifnot(length(train_cells_eval) > 100, length(test_cells_eval) > 100)

if (bench == 2) {
rownames(pred_train) <- gsub("\\.", "-", rownames(pred_train))
rownames(pred_test)  <- gsub("\\.", "-", rownames(pred_test))
}

map <- setNames(rownames(adt_mat), make.names(rownames(adt_mat)))
fix_names <- function(mat) {
  rn <- rownames(mat)
  hit <- rn %in% names(map)
  rn[hit] <- unname(map[rn[hit]])
  rownames(mat) <- rn
  mat
}

pred_train <- fix_names(pred_train)
pred_test  <- fix_names(pred_test)

# Now intersection gives 45
prot_use <- Reduce(intersect, list(rownames(adt_mat), rownames(pred_train), rownames(pred_test)))
length(prot_use) 

stopifnot(length(prot_use) > 0)

res_ctpnet <- lapply(prot_use, function(p){
  
  # 1) TRUE train (fit transform only on train)
  y_train_raw <- as.numeric(adt_mat[p, train_cells_eval])
  tf <- FMLE:::cap_and_scale_fit(y_train_raw, q=0.995)
  y_train_tf <- FMLE:::cap_and_scale_apply(y_train_raw, tf)
  
  # 2) TRUE test on same scale
  y_test_raw <- as.numeric(adt_mat[p, test_cells_eval])
  y_true     <- FMLE:::cap_and_scale_apply(y_test_raw, tf)
  
  # 3) cTPnet raw predictions
  yhat_train_raw <- as.numeric(pred_train[p, train_cells_eval])
  yhat_test_raw  <- as.numeric(pred_test[p,  test_cells_eval])
  
  ok_tr <- is.finite(y_train_tf) & is.finite(yhat_train_raw)
  if (sum(ok_tr) < 20 || sd(yhat_train_raw[ok_tr]) < 1e-8) {
    return(list(
      overall = data.frame(
        protein=p, MSE=NA, RMSE=NA, MAE=NA,
        Pearson=NA, Spearman=NA, R2=NA
      ),
      within_ct = data.frame(
        protein=p, cell_type=NA, Pearson=NA, n=NA
      )
    ))
  }
  
  
  # 4) train-only affine calibration (puts predictions onto FMLE scale)
  fit_cal <- lm(y_train_tf[ok_tr] ~ yhat_train_raw[ok_tr])
  yhat_test_cal <- as.numeric(coef(fit_cal)[1] + coef(fit_cal)[2] * yhat_test_raw)
  
  # 5) metrics
  ok <- is.finite(y_true) & is.finite(yhat_test_cal)
  yt <- y_true[ok]; yp <- yhat_test_cal[ok]
  
  mse  <- mean((yt-yp)^2)
  den  <- sum((yt-mean(yt))^2)
  overall_row <- data.frame(
    protein  = p,
    MSE      = mse,
    RMSE     = sqrt(mse),
    MAE      = mean(abs(yt-yp)),
    Pearson  = suppressWarnings(cor(yt, yp)),
    Spearman = suppressWarnings(cor(yt, yp, method="spearman")),
    R2       = if (den > 0) 1 - sum((yt-yp)^2)/den else NA_real_
  )
  ct_all <- meta[test_cells_eval, "cell_type", drop=TRUE]
  ct_ok  <- ct_all[ok]
  within_ct <- data.frame(
    protein   = p,
    cell_type = as.character(ct_ok),
    yt = yt,
    yp = yp,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::group_by(protein, cell_type) %>%
    dplyr::summarise(
      Pearson = suppressWarnings(cor(yt, yp)),
      n = dplyr::n(),
      .groups="drop"
    ) %>%
    dplyr::filter(cell_type != "Unassigned")
  
  # return both
  list(overall = overall_row, within_ct = within_ct)
  
})

#===============================================================================
# cross donor
# res_ctpnet <- dplyr::bind_rows(res_ctpnet)
# res_ctpnet_overall <- res_ctpnet$overall %>%
#   distinct(protein, .keep_all = TRUE)
# res_ctpnet_overall <- res_ctpnet_overall %>%
#   filter(
#     is.finite(Pearson),
#     is.finite(Spearman),
#     is.finite(R2)
#   )
# res_ctpnet_ct <- res_ctpnet$within_ct
# res_ctpnet_ct <- res_ctpnet_ct %>%
#   filter(is.finite(Pearson))
#===============================================================================
res_ctpnet_overall <- dplyr::bind_rows(lapply(res_ctpnet, `[[`, "overall"))
res_ctpnet_ct      <- dplyr::bind_rows(lapply(res_ctpnet, `[[`, "within_ct"))

write_csv(res_ctpnet_overall,
          file.path(ctp, "ctpnet_test_metrics.csv"))

write_csv(res_ctpnet_ct,
          file.path(ctp, "ctpnet_within_celltype_pearson.csv"))


# ============================================================
#  SAVE example-protein for scatter data plots 
# (truth + preds, matched test cells)
# ============================================================


# --- choose your  proteins here ---
# PBMC 10k Kaggle
ex_prots <- c("CD14",  "CD127-IL7Ra", "CD197-CCR7" ) 

# BMMC CITE-seq dataset
ex_prots <- c("CD40", "CD11c") 

# Young healthy BM
ex_prots <- c("CD11c-AB", "CD14-AB") 

# Cross donor
ex_prots <- c("anti-human-CD58-totalC", "anti-human-CD27-totalC")  


# where to save
fig_dir <- file.path(base, "paper_figures", ds)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
fig_out <- file.path(fig_dir, "fig_examples_same.rds")

# -----------------------------
# 0) Helpers
# -----------------------------
safe_name <- function(x) {
  x <- gsub("/", "_", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x
}
# -----------------------------
# 1) Load FMLE saved objects per protein
# -----------------------------
get_fmle_final_rds <- function(p) {
  f1 <- file.path(out_dir, paste0("final_", safe_name(p), ".rds"))
  if (file.exists(f1)) return(readRDS(f1))
  
  # try dot<->hyphen variants (covers HLA.DR vs HLA-DR)
  p2 <- gsub("\\.", "-", p)
  f2 <- file.path(out_dir, paste0("final_", safe_name(p2), ".rds"))
  if (file.exists(f2)) return(readRDS(f2))
  
  p3 <- gsub("-", ".", p, fixed = TRUE)
  f3 <- file.path(out_dir, paste0("final_", safe_name(p3), ".rds"))
  if (file.exists(f3)) return(readRDS(f3))
  
  stop("Missing FMLE final rds for protein: ", p,
       " tried: ", f1, " | ", f2, " | ", f3)
}

#====================================
# 2) Recompute scLinear predictions
#====================================
message("[FIG SAVE] Recomputing scLinear predictions for saving example scatter data...")

seu_prep2 <- seu_prep_hvg
stopifnot(all(train_cells %in% colnames(seu_prep2)), all(test_cells %in% colnames(seu_prep2)))

seu_train2 <- seu_prep2[, train_cells]
seu_test2  <- seu_prep2[, test_cells]

pipe2 <- scLinear::create_adt_predictor()

gexp_train2 <- Seurat::GetAssayData(seu_train2, assay="RNA", layer="counts")
adt_train2  <- Seurat::GetAssayData(seu_train2, assay="ADT", layer="counts")

pipe2 <- scLinear::fit_predictor(
  pipe          = pipe2,
  gexp_train    = as.matrix(gexp_train2),
  adt_train     = as.matrix(adt_train2),
  normalize_gex = TRUE,
  normalize_adt = FALSE
)

pred_scl_train2 <- scLinear::adt_predict(pipe2, gexp_train2, layer="counts", normalize=TRUE)
pred_scl_train2 <- as.matrix(pred_scl_train2@data)

gexp_test2 <- Seurat::GetAssayData(seu_test2, assay="RNA", layer="counts")
pred_scl_test2 <- scLinear::adt_predict(pipe2, gexp_test2, layer="counts", normalize=TRUE)
pred_scl_test2 <- as.matrix(pred_scl_test2@data)

# -----------------------------
# 3) loaded cTPnet matrices (pred_train/pred_test)
# -----------------------------
stopifnot(exists("pred_train"), exists("pred_test"))

# -----------------------------
# 4) Build per-protein aligned vectors and save
# -----------------------------
fig_examples <- list(
  proteins      = ex_prots,
  truth         = list(),
  pred_fmle     = list(),
  pred_sclinear = list(),
  pred_ctpnet   = list(),
  cells         = list()
)


# ============================================================
# cross donor in case of cTPnet
rownames(pred_train) <- gsub("\\.", "-", rownames(pred_train))
rownames(pred_test)  <- gsub("\\.", "-", rownames(pred_test))
# ============================================================

for (p in ex_prots) {
  
  # must exist everywhere with SAME key now
  if (!p %in% rownames(adt_mat)) {
    message("[FIG SAVE] Skip (not in adt_mat): ", p)
    next
  }
  if (!p %in% rownames(pred_train) || !p %in% rownames(pred_test)) {
    message("[FIG SAVE] Skip (not in pred_train/pred_test cTPnet): ", p)
    next
  }
  if (!p %in% rownames(pred_scl_train2) || !p %in% rownames(pred_scl_test2)) {
    message("[FIG SAVE] Skip (not in scLinear preds): ", p)
    next
  }
  
  # --- FMLE: load final object (has tf + fit + test_cells) ---
  fm <- get_fmle_final_rds(p)
  tf <- fm$tf_y
  
  # --- define common test cell set across truth + all preds ---
  test_cells_fmle <- fm$test_cells
  if (is.null(test_cells_fmle)) test_cells_fmle <- test_cells
  
  scl_cells <- colnames(pred_scl_test2)
  ctp_cells <- colnames(pred_test)
  
  common_cells <- Reduce(intersect, list(
    test_cells_fmle,
    colnames(adt_mat),
    scl_cells,
    ctp_cells
  ))
  common_cells <- test_cells_fmle[test_cells_fmle %in% common_cells]
  
  if (length(common_cells) < 50) {
    message("[FIG SAVE] Skip (too few aligned test cells) for ", p, ": ", length(common_cells))
    next
  }
  
  # --- truth on FMLE scale using stored tf ---
  y_test_raw <- as.numeric(adt_mat[p, common_cells])
  y_true     <- FMLE:::cap_and_scale_apply(y_test_raw, tf)
  
  # --- FMLE pred on common_cells ---
  pr_fm <- FMLE::fmle_predict(
    fm$fit,
    X_new = X[common_cells, , drop=FALSE],
    Z_new = Z[common_cells, , drop=FALSE]
  )
  yhat_fm <- as.numeric(pr_fm$mean)
  
  # --- scLinear: raw preds -> calibrate to FMLE scale using TRAIN only ---
  y_train_raw <- as.numeric(adt_mat[p, train_cells])
  y_train_tf  <- FMLE:::cap_and_scale_apply(y_train_raw, tf)
  
  tr_cells_ok <- intersect(train_cells, colnames(pred_scl_train2))
  if (length(tr_cells_ok) < 50) stop("Too few scLinear train cells aligned for saving: ", p)
  
  y_train_tf_aligned <- y_train_tf[match(tr_cells_ok, train_cells)]
  
  yhat_scl_tr <- as.numeric(pred_scl_train2[p, tr_cells_ok, drop=FALSE])
  yhat_scl_te <- as.numeric(pred_scl_test2[p, common_cells, drop=FALSE])
  
  ok_tr <- is.finite(y_train_tf_aligned) & is.finite(yhat_scl_tr)
  fit_cal_scl <- lm(y_train_tf_aligned[ok_tr] ~ yhat_scl_tr[ok_tr])
  yhat_scl_te_cal <- as.numeric(coef(fit_cal_scl)[1] + coef(fit_cal_scl)[2] * yhat_scl_te)
  
  # --- cTPnet: raw preds -> calibrate to FMLE scale using TRAIN only ---
  tr_cells_ctp <- intersect(train_cells, colnames(pred_train))
  if (length(tr_cells_ctp) < 50) stop("Too few cTPnet train cells aligned for saving: ", p)
  
  y_train_tf_ctp <- y_train_tf[match(tr_cells_ctp, train_cells)]
  
  yhat_ctp_tr <- as.numeric(pred_train[p, tr_cells_ctp, drop=FALSE])
  yhat_ctp_te <- as.numeric(pred_test[p, common_cells, drop=FALSE])
  
  ok_ctp <- is.finite(y_train_tf_ctp) & is.finite(yhat_ctp_tr)
  fit_cal_ctp <- lm(y_train_tf_ctp[ok_ctp] ~ yhat_ctp_tr[ok_ctp])
  yhat_ctp_te_cal <- as.numeric(coef(fit_cal_ctp)[1] + coef(fit_cal_ctp)[2] * yhat_ctp_te)
  
  # --- store ---
  fig_examples$truth[[p]]         <- y_true
  fig_examples$pred_fmle[[p]]     <- yhat_fm
  fig_examples$pred_sclinear[[p]] <- yhat_scl_te_cal
  fig_examples$pred_ctpnet[[p]]   <- yhat_ctp_te_cal
  fig_examples$cells[[p]]         <- common_cells
  
  message("[FIG SAVE] Saved vectors for ", p, " | n=", length(common_cells))
}

saveRDS(fig_examples, fig_out)
message("[FIG SAVE] Wrote: ", fig_out)







