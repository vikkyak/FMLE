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

# -----------------------------
# 0) Load fixed benchmark artifacts
# -----------------------------
bench <- 3  # set 1 or 2 or 3
# base  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_c_doror_%d", bench))
base <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_%d", bench))

ds      <- "citeseq_v1"
cv_dir  <- file.path(base, "FMLE", paste0(ds, "_cv"))

train_cells <- readRDS(file.path(base, ds, "train_cells.rds"))
test_cells  <- readRDS(file.path(base, ds, "test_cells.rds"))
X      <- readRDS(file.path(base, ds, "X.rds"))        # cells x genes
Z      <- readRDS(file.path(base, ds, "Z.rds"))        # cells x PCs
adt_mat<- readRDS(file.path(base, ds, "adt_mat.rds"))  # proteins x cells
if (bench == 1) {
  groups <- readRDS(file.path(base, ds, "groups.rds"))
} else {
  groups <- NULL
}

# Best hyperparams from CV (protein, R, m, lambda, mse_cv)
cv_best_all <- readr::read_csv(file.path(cv_dir, "cv_best_all_adts.csv"), show_col_types = FALSE)

# -----------------------------
# 1) Safety checks
# -----------------------------
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

# -----------------------------
# 2) Metrics helper
# -----------------------------
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

#----------------------------------------------------------------------------------#
# FMLE
#----------------------------------------------------------------------------------#
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

results_list <- vector("list", nrow(cv_best_all))
names(results_list) <- cv_best_all$protein

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
    R         = as.integer(best$R),
    m         = as.numeric(best$m),
    lambda_l1 = as.numeric(best$lambda),
    ridge     = 1e-6,
    standardize = TRUE,
    seed      = 1,
    verbose   = FALSE
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
  # --- extract betas for interpretation ---
  bb <- extract_beta(fit, feature_names = colnames(X_train))
  top <- topk_per_expert(bb$beta, k = 50)
  
  res_row <- data.frame(
    protein = prot,
    R       = as.integer(best$R),
    m       = as.numeric(best$m),
    lambda  = as.numeric(best$lambda),
    mse_cv  = as.numeric(best$mse_cv),
    met
  )
  
  results_list[[prot]] <- res_row
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
      test_cells = test_cells
    ),
    file.path(out_dir, paste0("final_", prot, ".rds"))
  )
}

results_fmle <- dplyr::bind_rows(results_list)
readr::write_csv(results_fmle, file.path(out_dir, "fmle_test_metrics.csv"))

print(results_fmle %>% arrange(MSE))


#----------------------------------------------------------------------------------#
# scLinear validation (train → predict test → metrics)
#----------------------------------------------------------------------------------#
scl    <- file.path(base, "sclinear") 

seu_prep <- readRDS(file.path(base, ds, "seu_final.rds"))  # <--- save this once
stopifnot(all(train_cells %in% colnames(seu_prep)), all(test_cells %in% colnames(seu_prep)))

seu_train <- seu_prep[, train_cells]
seu_test  <- seu_prep[, test_cells]

# ---- train scLinear predictor on TRAIN ----
pipe <- scLinear::create_adt_predictor()

gexp_train <- Seurat::GetAssayData(seu_train, assay="RNA", layer="counts")  # genes x train
adt_train  <- Seurat::GetAssayData(seu_train, assay="ADT", layer="counts")  # proteins x train

pipe <- scLinear::fit_predictor(
  pipe          = pipe,
  gexp_train    = as.matrix(gexp_train),
  adt_train     = as.matrix(adt_train),   # raw ADT counts
  normalize_gex = TRUE,
  normalize_adt = FALSE                  # IMPORTANT (we evaluate with FMLE scaling ourselves)
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
meta <- seu_prep@meta.data
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

# write files
readr::write_csv(res_sclinear_overall, file.path(scl, "scLinear_test_metrics_fmle_fair.csv"))
readr::write_csv(res_sclinear_ct,      file.path(scl, "scLinear_within_celltype_pearson.csv"))


#----------------------------------------------------------------------------------#
# ctpnet
#----------------------------------------------------------------------------------#

rm(pred_train, pred_test, train_cells_eval, test_cells_eval, prot_use )

ctp    <- file.path(base, "ctp") 

pred_train <- as.matrix(read.csv(file.path(ctp, "kaggle_ctpnet_pred_train.csv"),
                                 row.names=1, check.names=FALSE))
pred_test  <- as.matrix(read.csv(file.path(ctp, "kaggle_ctpnet_pred_test.csv"),
                                 row.names=1, check.names=FALSE))

stopifnot(length(intersect(train_cells, test_cells)) == 0)

train_cells_eval <- intersect(train_cells, colnames(pred_train))
test_cells_eval  <- intersect(test_cells,  colnames(pred_test))

# stopifnot(length(train_cells_eval) == length(train_cells))
# stopifnot(length(test_cells_eval)  == length(test_cells))


# stopifnot(all(train_cells %in% colnames(pred_train)))
# stopifnot(all(test_cells  %in% colnames(pred_test)))
# 
# train_cells_eval <- train_cells
# test_cells_eval  <- test_cells

cat("train_cells:", length(train_cells), "\n")
cat("train_cells_eval:", length(train_cells_eval), "\n")
cat("test_cells:", length(test_cells), "\n")
cat("test_cells_eval:", length(test_cells_eval), "\n")


stopifnot(length(train_cells_eval) > 100, length(test_cells_eval) > 100)

if (bench == 2) {
rownames(pred_train) <- gsub("\\.", "-", rownames(pred_train))
rownames(pred_test)  <- gsub("\\.", "-", rownames(pred_test))
}
# 
# prot_use <- Reduce(intersect, list(
#   rownames(adt_mat),
#   rownames(pred_train),
#   rownames(pred_test)
# ))
#  in case of kaggle 
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

res_ctpnet_overall <- dplyr::bind_rows(lapply(res_ctpnet, `[[`, "overall"))
res_ctpnet_ct      <- dplyr::bind_rows(lapply(res_ctpnet, `[[`, "within_ct"))

write_csv(res_ctpnet_overall,
          file.path(ctp, "ctpnet_test_metrics_FMLEscale.csv"))

write_csv(res_ctpnet_ct,
          file.path(ctp, "ctpnet_within_celltype_pearson.csv"))


# ================================
#   FIGURE 3 — SAME-DATASET
# ================================
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)

results_fmle <- read.csv("~/Desktop/FMLE/benchmarks_1/FMLE/citeseq_v1_final/fmle_test_metrics.csv",
                stringsAsFactors = FALSE)
results_sclinear <- read.csv("~/Desktop/FMLE/benchmarks_1/sclinear/scLinear_test_metrics_fmle_fair.csv",
                         stringsAsFactors = FALSE)
results_ctpnet  <- read.csv("~/Desktop/FMLE/benchmarks_1/ctp/ctpnet_test_metrics_FMLEscale.csv",
                         stringsAsFactors = FALSE)

# Panel A — Pearson distribution

df2 <- bind_rows(
  results_sclinear %>% select(protein, Pearson) %>% mutate(method="scLinear"),
  results_ctpnet   %>% select(protein, Pearson) %>% mutate(method="cTPnet"),
  results_fmle     %>% select(protein, Pearson) %>% mutate(method="FMLE")
) %>%
  mutate(method = factor(method, levels=c("scLinear","cTPnet","FMLE")))

df_wide <- df2 %>%
  pivot_wider(id_cols=protein, names_from=method, values_from=Pearson) %>%
  tidyr::drop_na()

p1 <- wilcox.test(df_wide$FMLE, df_wide$scLinear, paired=TRUE, alternative="greater")$p.value
p2 <- wilcox.test(df_wide$FMLE, df_wide$cTPnet,   paired=TRUE, alternative="greater")$p.value

P_cor <- ggplot(df2, aes(method, Pearson, fill=method)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(width=0.1, alpha=0.6) +
  theme_classic(base_size=14) +
  labs(
    title="Pearson across proteins",
    y="Pearson correlation",
    x=""
    # fill = sprintf(
    #   "paired Wilcoxon\nFMLE > scLinear: p=%.2g\nFMLE > cTPnet: p=%.2g",
    #   p1, p2
    # )
  ) +
  theme(plot.title = element_text(hjust = 0.5))



ggsave(
  filename = file.path(out_dir, "Pearson_correlation.pdf"),
  plot = Pearson_correlation,
  device = "pdf",
  width = 8,
  height = 6,
  units = "in"
)


# Panel B — FMLE vs scLinear per protien

# 1) FMLE: build per-(protein, cell_type) Pearson on TEST

meta <- seu_prep@meta.data  # you already have this
files <- list.files(out_dir, pattern="^final_.*\\.rds$", full.names=TRUE)

results_fmle_ct <- lapply(files, function(f){
  obj <- readRDS(f)
  
  prot  <- obj$protein
  cells <- obj$test_cells
  
  keep <- intersect(cells, rownames(meta))
  if (length(keep) == 0) return(NULL)
  
  idx <- match(keep, cells)
  
  df <- data.frame(
    protein   = prot,
    cell      = keep,
    cell_type = meta[keep, "cell_type", drop=TRUE],
    yt        = as.numeric(obj$y_test)[idx],
    yp        = as.numeric(obj$yhat_test)[idx],
    stringsAsFactors = FALSE
  )
  
  df %>%
    group_by(protein, cell_type) %>%
    summarise(
      Pearson_fmle = suppressWarnings(cor(yt, yp, use="pairwise.complete.obs")),
      n = n(),
      .groups="drop"
    )
}) %>% bind_rows()

readr::write_csv(results_fmle_ct, file.path(out_dir, "fmle_within_celltype_pearson.csv"))


# 2) scLinear: build per-(protein, cell_type) Pearson on TEST

stopifnot(all(test_cells_eval %in% rownames(meta)))
ct_test <- meta[test_cells_eval, "cell_type", drop = TRUE]

# proteins you evaluated
prot_use <- intersect(rownames(pred_test), rownames(adt_mat))
stopifnot(length(prot_use) > 0)

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
    # keep your original return unchanged
    return(list(
      overall = data.frame(protein=p, MSE=NA, RMSE=NA, MAE=NA,
                           Pearson=NA, Spearman=NA, R2=NA),
      within_ct = data.frame(protein=p, cell_type=NA, Pearson=NA)
    ))
  }
  
  fit_cal <- lm(y_train_tf[ok] ~ yhat_train_raw[ok])
  y_hat <- as.numeric(coef(fit_cal)[1] + coef(fit_cal)[2] * yhat_test_raw)
  
  ok2 <- is.finite(y_true) & is.finite(y_hat)
  yt <- y_true[ok2]; yp <- y_hat[ok2]
  
  mse <- mean((yt-yp)^2)
  den <- sum((yt-mean(yt))^2)
  
  # ---- your existing overall row (UNCHANGED values) ----
  overall_row <- data.frame(
    protein=p,
    MSE=mse,
    RMSE=sqrt(mse),
    MAE=mean(abs(yt-yp)),
    Pearson=suppressWarnings(cor(yt, yp)),
    Spearman=suppressWarnings(cor(yt, yp, method="spearman")),
    R2=if (den>0) 1 - sum((yt-yp)^2)/den else NA_real_
  )
  
  # ---- ADDED: within cell type Pearson on TEST (panel B-left) ----
  # align the same filtered test points (ok2) to their cell types
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

# keep your original res_sclinear exactly (overall metrics)
res_sclinear_overall <- dplyr::bind_rows(lapply(res_sclinear, `[[`, "overall"))
res_sclinear_ct      <- dplyr::bind_rows(lapply(res_sclinear, `[[`, "within_ct"))

# write files
readr::write_csv(res_sclinear_overall, file.path(scl, "scLinear_test_metrics_fmle_fair.csv"))
readr::write_csv(res_sclinear_ct,      file.path(scl, "scLinear_within_celltype_pearson.csv"))




scl_ct <- readr::read_csv(file.path(scl, "scLinear_within_celltype_pearson.csv")) %>%
  dplyr::rename(Pearson_scl = Pearson)

fmle_ct <- readr::read_csv(file.path(out_dir, "fmle_within_celltype_pearson.csv"))

cmp_fmle_scl <- inner_join(fmle_ct, scl_ct, by=c("protein","cell_type"))

cmp_fmle_scl <- cmp_fmle_scl %>% dplyr::filter(cell_type != "Unassigned")

p_fmle_scl <- ggplot(cmp_fmle_scl, aes(Pearson_scl, Pearson_fmle, color=protein, shape=cell_type)) +
  geom_point(size=3, alpha=0.9) +
  geom_abline(slope=1, intercept=0, color="grey40") +
  coord_equal(xlim=c(-0.5,1), ylim=c(-0.5,1)) +
  theme_light(base_size=14) +
  labs(x="scLinear", y="FMLE", color="Protein", shape="Cell type")

ctpnet_ct <- readr::read_csv(file.path(ctp, "ctpnet_within_celltype_pearson.csv")) %>%
  dplyr::rename(Pearson_ctp = Pearson)


cmp_fmle_ctpnet <- inner_join(fmle_ct,
                     ctpnet_ct,
                     by=c("protein","cell_type"))


p_fmle_ctp <- ggplot(cmp_fmle_ctpnet, aes(Pearson_ctp, Pearson_fmle, color=protein, shape=cell_type)) +
  geom_point(size=3, alpha=0.9) +
  geom_abline(slope=1, intercept=0, color="grey40") +
  coord_equal(xlim=c(-0.5,1), ylim=c(-0.5,1)) +
  theme_light(base_size=14) +
  labs(x="cTPnet", y="FMLE", color="Protein", shape="Cell type")



library(cowplot) 

p1 <- p_fmle_scl + theme(legend.position = "none")
p2 <- p_fmle_ctp + theme(legend.position = "none")

leg <- cowplot::get_legend(
  p_fmle_scl +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.direction = "horizontal"
    ) +
    guides(
      color = guide_legend(ncol = 4, byrow = TRUE),  # increase to fill width
      shape = guide_legend(nrow = 2)
    )
)

library(ggplot2)
library(cowplot)

# 1) Boxplot (no legend)
P_cor2 <- P_cor +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(0.3, 1.02))

# 2) Scatters: SAME limits, NO coord_equal (so height matches)
p_scl2 <- p_fmle_scl +
  coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) +
  theme(legend.position = "none")

p_ctp2 <- p_fmle_ctp +
  coord_cartesian(xlim = c(-0.5, 1), ylim = c(-0.5, 1)) +
  theme(legend.position = "none")

# 3) One shared legend (Protein row + Cell type row)
leg <- cowplot::get_legend(
  p_fmle_scl +
    theme(
      legend.position  = "bottom",
      legend.box       = "vertical",
      legend.direction = "horizontal",
      legend.box.just  = "center"
    ) +
    guides(
      color = guide_legend(ncol = 3, byrow = TRUE, title.position = "top"),
      shape = guide_legend(nrow = 1,  byrow = TRUE, title.position = "top")
    )
)

# 4) Top row aligned
# top_row <- cowplot::plot_grid(
#   P_cor2, p_scl2, p_ctp2,
#   nrow = 1,
#   align = "hv",
#   axis = "tb",
#   rel_widths = c(1, 1, 1)
# )

top_row <- cowplot::plot_grid(
  p_scl2, p_ctp2,
  ncol = 1,
  align = "hv",
  axis = "tb",
  rel_widths = c(1, 1, 1)
)



# 5) Final (legend below)
final_plot <- cowplot::plot_grid(
  top_row, leg,
  ncol = 1,
  rel_heights = c(0.5, 0.22)
)

final_plot



# 
# 
# 
# leg <- cowplot::get_legend(
#   p_scl2 + theme(
#     legend.position = "bottom",
#     legend.box = "horizontal",
#     legend.direction = "horizontal"
#     # legend.title = element_blank()
#   )
# )
# 
# p1 <- p_scl2 + theme(legend.position = "none")
# p2 <- p_ctp2 + theme(legend.position = "none")
# 
# # 3) one column, three rows (legend is the 3rd row)
# final_plot <- cowplot::plot_grid(
#   p1,
#   p2,
#   leg,
#   ncol = 1,
#   align = "v",
#   axis = "lr",
#   rel_heights = c(1, 1, 0.18)  # adjust legend row height
# )
# 
# final_plot
# 
# 
# 
# 
# 
# 
# 
# 
# 
# p_leg_src <- p_scl2 +
#   guides(
#     color = guide_legend(
#       title = "Protein",
#       nrow = 2,        # <- key: use horizontal space
#       byrow = TRUE
#     ),
#     shape = guide_legend(
#       title = "Cell type",
#       nrow = 1,
#       byrow = TRUE
#     )
#   ) +
#   theme(
#     legend.position  = "bottom",
#     legend.box       = "vertical",   # Protein block then Cell type block
#     legend.direction = "horizontal",
#     legend.key.size  = unit(0.35, "cm"),
#     legend.spacing.x = unit(0.25, "cm"),
#     legend.spacing.y = unit(0.10, "cm"),
#     legend.text      = element_text(size = 9),
#     legend.title     = element_text(size = 10, face = "bold"),
#     legend.margin    = margin(t = 0, r = 0, b = 0, l = 0)
#   )
# 
# leg <- cowplot::get_legend(p_leg_src)
# 
# # 2) Remove legends from panels
# p1 <- p_scl2 + theme(legend.position = "none")
# p2 <- p_ctp2 + theme(legend.position = "none")
# 
# # 3) One column, three rows (legend is 3rd row)
# final_plot <- cowplot::plot_grid(
#   p1,
#   p2,
#   leg,
#   ncol = 1,
#   align = "v",
#   axis = "lr",
#   rel_heights = c(1, 1, 0.30)  # increase a bit so legend doesn't get squashed
# )
# 
# final_plot

# 
# # FMLE: one protein + one cell type panel
# 
# plot_c_fmle <- function(prot, ct, method, out_dir, meta, n_max = 5000) {
#   obj <- readRDS(file.path(out_dir, paste0("final_", prot, ".rds")))
#   cells <- obj$test_cells
#   df <- data.frame(
#     cell = cells,
#     truth = as.numeric(obj$y_test),
#     pred  = as.numeric(obj$yhat_test),
#     cell_type = meta[cells, "cell_type", drop=TRUE],
#     stringsAsFactors = FALSE
#   )
#   df <- df[df$cell_type == ct & is.finite(df$truth) & is.finite(df$pred), , drop=FALSE]
#   if (nrow(df) > n_max) df <- df[sample.int(nrow(df), n_max), ]
#   
#   r <- suppressWarnings(cor(df$pred, df$truth))
#   ggplot(df, aes(x = pred, y = truth)) +
#     geom_point(alpha = 0.6, size = 1.6) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#     theme_classic(base_size = 14) +
#     labs(title = paste0(method, " — ",ct, " — ", prot),
#          subtitle = sprintf("cor: %.2f", r),
#          x = "Prediction", y = "Truth")
# }
# 
# # example:
# meta <- seu_prep@meta.data
# # pC <- plot_c_fmle(prot="CD11c", ct="Monocyte", out_dir=out_dir, meta=meta)
# # pC
# 
# # scLinear 
# pred_train <- scLinear::adt_predict(pipe, gexp_train, layer ="counts", normalize=TRUE)  # proteins x test
# pred_train <- as.matrix(pred_train@data)
# 
# # ---- predict TEST ----
# gexp_test <- Seurat::GetAssayData(seu_test, assay="RNA", layer="counts")
# pred_test <- scLinear::adt_predict(pipe, gexp_test, layer ="counts", normalize=TRUE)  # proteins x test
# pred_test <- as.matrix(pred_test@data)
# 
# # ---- align cells & proteins ----
# train_cells_eval <- intersect(train_cells, colnames(pred_train))
# test_cells_eval  <- intersect(test_cells,  colnames(pred_test))
# stopifnot(length(train_cells_eval) > 0, length(test_cells_eval) > 0)
# 
# pred_train <- pred_train[, train_cells_eval, drop=FALSE]
# pred_test  <- pred_test[,  test_cells_eval,  drop=FALSE]
# 
# # cTPnet
# rm(pred_train, pred_test, train_cells_eval, test_cells_eval, prot_use )
# 
# ctp    <- file.path(base, "ctp") 
# 
# pred_train <- as.matrix(read.csv(file.path(ctp, "kaggle_ctpnet_pred_train.csv"),
#                                  row.names=1, check.names=FALSE))
# pred_test  <- as.matrix(read.csv(file.path(ctp, "kaggle_ctpnet_pred_test.csv"),
#                                  row.names=1, check.names=FALSE))
# 
# stopifnot(length(intersect(train_cells, test_cells)) == 0)
# 
# train_cells_eval <- intersect(train_cells, colnames(pred_train))
# test_cells_eval  <- intersect(test_cells,  colnames(pred_test))
# 
# # stopifnot(length(train_cells_eval) == length(train_cells))
# # stopifnot(length(test_cells_eval)  == length(test_cells))
# 
# 
# # stopifnot(all(train_cells %in% colnames(pred_train)))
# # stopifnot(all(test_cells  %in% colnames(pred_test)))
# # 
# # train_cells_eval <- train_cells
# # test_cells_eval  <- test_cells
# 
# cat("train_cells:", length(train_cells), "\n")
# cat("train_cells_eval:", length(train_cells_eval), "\n")
# cat("test_cells:", length(test_cells), "\n")
# cat("test_cells_eval:", length(test_cells_eval), "\n")
# 
# 
# stopifnot(length(train_cells_eval) > 100, length(test_cells_eval) > 100)
# 
# if (bench == 2) {
#   rownames(pred_train) <- gsub("\\.", "-", rownames(pred_train))
#   rownames(pred_test)  <- gsub("\\.", "-", rownames(pred_test))
# }
# # 
# # prot_use <- Reduce(intersect, list(
# #   rownames(adt_mat),
# #   rownames(pred_train),
# #   rownames(pred_test)
# # ))
# #  in case of kaggle 
# map <- setNames(rownames(adt_mat), make.names(rownames(adt_mat)))
# fix_names <- function(mat) {
#   rn <- rownames(mat)
#   hit <- rn %in% names(map)
#   rn[hit] <- unname(map[rn[hit]])
#   rownames(mat) <- rn
#   mat
# }
# 
# 
# pred_train <- fix_names(pred_train)
# pred_test  <- fix_names(pred_test)
# 
# plot_c_baseline <- function(prot, ct, method, pred_train, pred_test, adt_mat,
#                             train_cells_eval, test_cells_eval, meta, n_max = Inf) {
#   
#   y_train_raw <- as.numeric(adt_mat[prot, train_cells_eval])
#   tf <- FMLE:::cap_and_scale_fit(y_train_raw, q=0.995)
#   y_train_tf <- FMLE:::cap_and_scale_apply(y_train_raw, tf)
#   
#   y_test_raw <- as.numeric(adt_mat[prot, test_cells_eval])
#   y_true <- FMLE:::cap_and_scale_apply(y_test_raw, tf)
#   
#   yhat_train_raw <- as.numeric(pred_train[prot, train_cells_eval])
#   yhat_test_raw  <- as.numeric(pred_test[prot,  test_cells_eval])
#   
#   ok_tr <- is.finite(y_train_tf) & is.finite(yhat_train_raw)
#   fit_cal <- lm(y_train_tf[ok_tr] ~ yhat_train_raw[ok_tr])
#   y_hat <- as.numeric(coef(fit_cal)[1] + coef(fit_cal)[2] * yhat_test_raw)
#   
#   ct_all <- meta[test_cells_eval, "cell_type", drop=TRUE]
#   
#   df <- data.frame(cell=test_cells_eval, truth=y_true, pred=y_hat, cell_type=ct_all)
#   df <- df[df$cell_type==ct & is.finite(df$truth) & is.finite(df$pred), , drop=FALSE]
#   if (nrow(df) > n_max) df <- df[sample.int(nrow(df), n_max), ]
#   
#   r <- suppressWarnings(cor(df$pred, df$truth, use="pairwise.complete.obs"))
#   
#   ggplot(df, aes(pred, truth)) +
#     geom_point(alpha=0.6, size=1.6) +
#     geom_abline(slope=1, intercept=0, linetype="dashed") +
#     theme_classic(base_size=14) +
#     labs(
#       title = paste0(method, " : ", ct, " — ", prot),
#       subtitle = sprintf("cor: %.2f", r),
#       x="Prediction", y="Truth"
#     )
# }
# 
# 
# pred_train_sc <- pred_train
# pred_test_sc  <- pred_test
# train_cells_sc <- train_cells_eval  # train_cells_eval
# test_cells_sc  <-test_cells_eval    # test_cells_eval
# 
# library(cowplot)
# p_scl <- plot_c_baseline("CD197-CCR7","T","scLinear",
#                          pred_train_sc, pred_test_sc, adt_mat,
#                          train_cells_sc, test_cells_sc, meta)
# 
# pred_train_ctp   <- pred_train
# pred_test_ctp    <- pred_test  
# train_cells_ctp  <- train_cells_eval  # train_cells_eval
# test_cells_ctp   <-test_cells_eval    # test_cells_eval
# 
# 
# p_ctp <- plot_c_baseline("CD69","T","cTPnet",
#                          pred_train_ctp, pred_test_ctp, adt_mat,
#                          train_cells_ctp, test_cells_ctp, meta)
# 
# pC <- plot_c_fmle(prot="CD197-CCR7", ct="T", "FMLE", out_dir=out_dir, meta=meta)
# 
# cowplot::plot_grid(pC, p_scl, p_ctp, nrow=1)
# 
# 













# cmp_scl <- inner_join(
#   results_fmle %>% select(protein, Pearson_fmle=Pearson),
#   results_sclinear %>% select(protein, Pearson_scl=Pearson),
#   by="protein"
# )
# 
# FMLL_scLinear <- ggplot(cmp_scl, aes(Pearson_scl, Pearson_fmle)) +
#   geom_point(size=3) +
#   geom_abline(slope=1, intercept=0, linetype="dashed") +
#   geom_text_repel(aes(label=protein), size=3) +
#   theme_classic(base_size=14) +
#   labs(x="scLinear Pearson", y="FMLE Pearson")
# 
# ggsave(
#   filename = file.path(out_dir, "FMLL_scLinear.pdf"),
#   plot = FMLL_scLinear,
#   device = "pdf",
#   width = 8,
#   height = 6,
#   units = "in"
# )
# 
# # Panel C — FMLE vs cTPnet
# cmp_ctp <- inner_join(
#   results_fmle %>% select(protein, Pearson_fmle=Pearson),
#   results_ctpnet %>% select(protein, Pearson_ctp=Pearson),
#   by="protein"
# )
# 
# FMLE_cTPnet <- ggplot(cmp_ctp, aes(Pearson_ctp, Pearson_fmle)) +
#   geom_point(size=3) +
#   geom_abline(slope=1, intercept=0, linetype="dashed") +
#   geom_text_repel(aes(label=protein), size=3) +
#   theme_classic(base_size=14) +
#   labs(x="cTPnet Pearson", y="FMLE Pearson")
# 
# 
# 
# ggsave(
#   filename = file.path(out_dir, "FMLE_cTPnet.pdf"),
#   plot = FMLE_cTPnet,
#   device = "pdf",
#   width = 8,
#   height = 6,
#   units = "in"
# )
# 
# P_main <- Pearson_correlation|FMLL_scLinear|FMLE_cTPnet 
# 
# ggsave(
#   filename = file.path(out_dir, "P_main.pdf"),
#   plot = P_main,
#   device = "pdf",
#   width = 15,
#   height = 5,
#   units = "in"
# )
# 
# 
# # Panel D — Predicted vs true (example proteins)
# plot_protein <- function(p){
#   df <- data.frame(
#     truth = y_true_list[[p]],
#     FMLE  = y_pred_fmle[[p]],
#     scL   = y_pred_scl[[p]]
#   )
#   
#   ggplot(df, aes(truth, FMLE)) +
#     geom_point(alpha=0.5) +
#     geom_abline(slope=1, intercept=0, linetype="dashed") +
#     theme_classic() +
#     labs(title=p)
# }
# 
# plot_protein("CD3")
# plot_protein("CD14")


# ============================================================
# FIGURE 3 (Panel D) — SAVE example-protein scatter DATA
# (truth + preds on FMLE scale, matched test cells)
# ============================================================

# ============================================================
# FIGURE 3 (Panel D) — SAVE example-protein scatter DATA
# (truth + preds on FMLE scale, matched test cells)
# ============================================================

# --- choose your 3 proteins here ---
# ex_prots <- c("CD14", "CD127-IL7Ra", "HLA.DR")
# ex_prots <- c("CD197-CCR7", "CD34", "CD127-IL7Ra")
ex_prots <- c("CD14", "CD19", "CD127-IL7Ra", "CD197-CCR7", "CD69" )

# where to save
fig_dir <- file.path(base, "paper_figures", ds)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
fig3_out <- file.path(fig_dir, "fig3_examples_same.rds")

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

# -----------------------------
# 2) Recompute scLinear predictions
# -----------------------------
message("[FIG3 SAVE] Recomputing scLinear predictions for saving example scatter data...")

seu_prep2 <- seu_prep
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
# 3) Use your already-loaded cTPnet matrices (pred_train/pred_test)
# -----------------------------
stopifnot(exists("pred_train"), exists("pred_test"))

# -----------------------------
# 4) Build per-protein aligned vectors and save
# -----------------------------
fig3_examples <- list(
  proteins      = ex_prots,
  truth         = list(),
  pred_fmle     = list(),
  pred_sclinear = list(),
  pred_ctpnet   = list(),
  cells         = list()
)

for (p in ex_prots) {
  
  # must exist everywhere with SAME key now
  if (!p %in% rownames(adt_mat)) {
    message("[FIG3 SAVE] Skip (not in adt_mat): ", p)
    next
  }
  if (!p %in% rownames(pred_train) || !p %in% rownames(pred_test)) {
    message("[FIG3 SAVE] Skip (not in pred_train/pred_test): ", p)
    next
  }
  if (!p %in% rownames(pred_scl_train2) || !p %in% rownames(pred_scl_test2)) {
    message("[FIG3 SAVE] Skip (not in scLinear preds): ", p)
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
    message("[FIG3 SAVE] Skip (too few aligned test cells) for ", p, ": ", length(common_cells))
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
  fig3_examples$truth[[p]]         <- y_true
  fig3_examples$pred_fmle[[p]]     <- yhat_fm
  fig3_examples$pred_sclinear[[p]] <- yhat_scl_te_cal
  fig3_examples$pred_ctpnet[[p]]   <- yhat_ctp_te_cal
  fig3_examples$cells[[p]]         <- common_cells
  
  message("[FIG3 SAVE] Saved vectors for ", p, " | n=", length(common_cells))
}

saveRDS(fig3_examples, fig3_out)
message("[FIG3 SAVE] Wrote: ", fig3_out)


# ex <- readRDS("~/Desktop/FMLE/benchmarks_1/paper_figures/citeseq_v1/fig3_examples_same.rds")
# names(ex$truth)
# 
# p <- "CD14"
# df <- data.frame(
#   truth = ex$truth[[p]],
#   FMLE  = ex$pred_fmle[[p]],
#   scL   = ex$pred_sclinear[[p]],
#   cTP   = ex$pred_ctpnet[[p]]
# )
# 
# ggplot(df, aes(truth, FMLE)) + geom_point(alpha=0.4) + geom_abline(lty=2) + theme_classic()






# ============================================================
# FIGURE (Main paper) — Panel D
# Example-protein scatter (Real vs FMLE-predicted), colored by cell type
# Uses your saved RDS: fig3_examples_same.rds
# Style matched to scLinear: y=x dashed identity line + clean theme
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# -----------------------------
# Inputs you already have
# -----------------------------
# fig3_out  : path to fig3_examples_same.rds (you created earlier)
# seu_prep  : Seurat object (or any object that has celltype per barcode)
# base/ds/out_dir are irrelevant here except for saving

# -----------------------------
# 1) Load saved example vectors
# -----------------------------
fig3 <- readRDS(fig3_out)

# -----------------------------
# 2) Build a named celltype vector (edit the meta column name!)
# -----------------------------
# Change "celltype" to your metadata column if different
stopifnot(exists("seu_prep"))
stopifnot("cell_type" %in% colnames(seu_prep@meta.data))
celltype_vec <- setNames(as.character(seu_prep@meta.data$cell_type), colnames(seu_prep))

# -----------------------------
# 3) Plot function (paper style)
# -----------------------------
# plot_panelD_fmle <- function(prot,
#                              fig3,
#                              celltype_vec,
#                              n_points = 6000,
#                              add_lm = FALSE) {
#   
#   stopifnot(prot %in% names(fig3$truth),
#             prot %in% names(fig3$pred_fmle),
#             prot %in% names(fig3$cells))
#   
#   cells <- fig3$cells[[prot]]
#   
#   df <- data.frame(
#     pred = as.numeric(fig3$pred_fmle[[prot]]),
#     real = as.numeric(fig3$truth[[prot]]),
#     cell = cells,
#     stringsAsFactors = FALSE
#   ) %>%
#     filter(is.finite(pred), is.finite(real))
#   df <- df[sample(nrow(df), 1500), ]
#   # attach cell types (NA if missing)
#   df$celltype <- unname(celltype_vec[df$cell])
#   df$celltype[is.na(df$celltype)] <- "Unassigned"
#   df <- df[df$celltype != "Unassigned", , drop=FALSE]
#   df$celltype <- factor(df$celltype)
#   df$celltype <- droplevels(df$celltype)
#   # subsample for visual clarity
#   if (nrow(df) > n_points) {
#     set.seed(1)
#     df <- df[sample.int(nrow(df), n_points), ]
#   }
#   
#   r <- suppressWarnings(cor(df$pred, df$real, method="pearson"))
#   rho <- suppressWarnings(cor(df$pred, df$real, method="spearman"))
#   
#   lim <- range(c(df$pred, df$real), na.rm = TRUE)
#   
#   p <- ggplot(df, aes(x = pred, y = real, color = celltype)) +
#     geom_point(alpha = 0.75, size = 1.2) +
#     geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.7, color = "grey25") +
#     coord_equal(xlim = lim, ylim = lim) +
#     theme_classic(base_size = 14) +
#     labs(
#       title = sprintf("%s (r=%.2f, ρ=%.2f)", prot, r, rho),
#       x = "FMLE predicted (FMLE scale)",
#       y = "Real protein (FMLE scale)",
#       color = NULL
#     ) +
#     theme(
#       plot.title = element_text(hjust = 0, face = "bold"),
#       legend.position = "right",
#       legend.key.height = unit(0.35, "cm"),
#       legend.text = element_text(size = 10)
#     )
#   
#   if (add_lm) {
#     p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black", inherit.aes = FALSE,
#                          aes(x = pred, y = real))
#   }
#   
#   p
# }

# -----------------------------
# 4) Choose proteins that actually exist in your saved RDS
# -----------------------------
# Use proteins present in fig3 (this avoids your earlier CD3 missing issue)
available <- names(fig3$truth)
print(available)

# Pick 3 for main paper (edit these to your preferred ones FROM `available`)
# Example: if only CD14 exists, you MUST regenerate fig3_out with more proteins.
# ex_prots <- c("CD14", "CD127-IL7Ra", "HLA.DR")
stopifnot(all(ex_prots %in% available))
if (length(ex_prots) == 0) stop("None of requested proteins exist in fig3_out. Available: ", paste(available, collapse=", "))

# -----------------------------
# 5) Build panel (1×N)
# -----------------------------
plots <- lapply(ex_prots, plot_panelD_fmle, fig3 = fig3, celltype_vec = celltype_vec, n_points = 6000)
panelD <- wrap_plots(plots, nrow = 1)

panelD

# -----------------------------
# 6) Save for main paper
# -----------------------------
out_pdf <- file.path(dirname(fig3_out), "Fig_panelD_FMLE_examples.pdf")
ggsave(out_pdf, panelD, width = 5.5 * length(ex_prots), height = 5.5, units = "in")
message("Wrote: ", out_pdf)


plot_panelD_method <- function(prot,
                               method = c("FMLE","scLinear","cTPnet"),
                               fig3,
                               celltype_vec,
                               n_points = 6000,
                               add_lm = FALSE) {
  
  method <- match.arg(method)
  
  pred_list <- switch(
    method,
    "FMLE"     = fig3$pred_fmle,
    "scLinear" = fig3$pred_sclinear,
    "cTPnet"   = fig3$pred_ctpnet
  )
  
  stopifnot(prot %in% names(fig3$truth),
            prot %in% names(pred_list),
            prot %in% names(fig3$cells))
  
  cells <- fig3$cells[[prot]]
  
  df <- data.frame(
    pred = as.numeric(pred_list[[prot]]),
    real = as.numeric(fig3$truth[[prot]]),
    cell = cells,
    stringsAsFactors = FALSE
  )
  
  df <- df[is.finite(df$pred) & is.finite(df$real), , drop=FALSE]
  
  # attach cell types
  df$celltype <- unname(celltype_vec[df$cell])
  df$celltype[is.na(df$celltype)] <- "Unassigned"
  df <- df[df$celltype != "Unassigned", , drop=FALSE]
  df$celltype <- droplevels(factor(df$celltype))
  
  # subsample for visual clarity
  if (nrow(df) > n_points) {
    set.seed(1)
    df <- df[sample.int(nrow(df), n_points), , drop=FALSE]
  }
  
  r   <- suppressWarnings(cor(df$pred, df$real, method="pearson"))
  rho <- suppressWarnings(cor(df$pred, df$real, method="spearman"))
  
  lim <- range(c(df$pred, df$real), na.rm=TRUE)
  
  p <- ggplot(df, aes(x = pred, y = real, color = celltype)) +
    geom_point(alpha = 0.75, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.7, color = "grey25") +
    coord_equal(xlim = lim, ylim = lim) +
    theme_classic(base_size = 14) +
    labs(
      title = sprintf("%s  (r=%.2f, ρ=%.2f)", prot,  r, rho),
      x = sprintf("%s predicted", method),
      y = "Real protein",
      color = NULL
    ) +
    theme(plot.title = element_text(hjust = 0, face = "bold"),
          legend.position = "right")
  
  if (add_lm) {
    p <- p + geom_smooth(method="lm", se=FALSE, linewidth=0.6, color="black")
  }
  
  p
}


p1 <- plot_panelD_method("CD197-CCR7", "FMLE", fig3, celltype_vec)
p2 <- plot_panelD_method("CD197-CCR7", "scLinear", fig3, celltype_vec)
p3 <- plot_panelD_method("CD197-CCR7", "cTPnet", fig3, celltype_vec)

(p1 | p2 | p3)

p_CD197 <- (p1 / p2 / p3)


# ================================
#  multidatasets 
# ================================

library(dplyr)
library(readr)
library(ggplot2)

# -------- SETTINGS --------
fmle_ds_by_bench <- c(
  "citeseq_v1_final",  # bench 1 FMLE folder
  "citeseq_v1_final",  # bench 2 FMLE folder (edit if different)
  "citeseq_v1_final"   # bench 3 FMLE folder (edit if different)
)

# -------- READ FUNCTION --------
read_one_bench <- function(bench, fmle_ds){
  
  base <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_%d", bench))
  
  fmle_path <- file.path(base, "FMLE", fmle_ds, "fmle_test_metrics.csv")
  scl_path  <- file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv")
  ctp_path  <- file.path(base, "ctp", "ctpnet_test_metrics_FMLEscale.csv")
  
  stopifnot(file.exists(fmle_path),
            file.exists(scl_path),
            file.exists(ctp_path))
  
  bind_rows(
    read_csv(fmle_path, show_col_types=FALSE) %>%
      transmute(dataset=paste0("bench_",bench),
                method="FMLE",
                protein,
                Pearson),
    
    read_csv(scl_path, show_col_types=FALSE) %>%
      transmute(dataset=paste0("bench_",bench),
                method="scLinear",
                protein,
                Pearson),
    
    read_csv(ctp_path, show_col_types=FALSE) %>%
      transmute(dataset=paste0("bench_",bench),
                method="cTPnet",
                protein,
                Pearson)
  )
}

# -------- BUILD DATA --------
df_all <- bind_rows(
  read_one_bench(1, fmle_ds_by_bench[1]),
  read_one_bench(2, fmle_ds_by_bench[2]),
  read_one_bench(3, fmle_ds_by_bench[3])
) %>% 
  filter(is.finite(Pearson))%>%
  mutate(method = factor(method, levels=c("scLinear","cTPnet","FMLE")))

# -------- PANEL A PLOT --------
pA <- ggplot(df_all,
             aes(x=method, y=Pearson, color=method)) +
  
  geom_boxplot(outlier.shape=NA, width=0.6) +
  
  geom_jitter(width=0.15, alpha=0.7, size=2) +
  
  facet_wrap(~dataset, ncol=1) +
  
  theme_classic(base_size=14) +
  
  labs(
    x=NULL,
    y="Per-protein Pearson correlation"
  ) +
  
  theme(
    legend.position="none",
    strip.text=element_text(face="bold")
  )


final <- (pA | (p_scl / p_ctp) | p_CD14 | p_CD127 | p_CD197) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 8)
  )

final






library(patchwork)
library(grid)  # unit()

kill_leg <- theme(legend.position = "none")

pA0      <- pA      + kill_leg
p_scl0   <- p_scl2  + kill_leg
p_ctp0   <- p_ctp2  + kill_leg
p_CD140  <- p_CD14  + kill_leg
p_CD1270 <- p_CD127 + kill_leg
p_CD1970 <- p_CD197 + kill_leg
# 1) FORCE WRAPPING (apply to the plots that own the legends)
p_scl0 <- p_scl0 + guides(color = guide_legend(nrow = 2, byrow = TRUE))
p_ctp0 <- p_ctp0 + guides(color = guide_legend(nrow = 2, byrow = TRUE))

# force cell-type (shape) legend to wrap instead of being clipped
p_CD140  <- p_CD140  + guides(shape = guide_legend(nrow = 1, byrow = TRUE))
p_CD1270 <- p_CD1270 + guides(shape = guide_legend(nrow = 1, byrow = TRUE))
p_CD1970 <- p_CD1970 + guides(shape = guide_legend(nrow = 1, byrow = TRUE))

final <- (pA0 | (p_scl0 / p_ctp0) | p_CD140 | p_CD1270 | p_CD1970) +
  plot_layout(guides = "collect") &
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",     # <<< key: separate rows so celltype isn't cut
    legend.box.just   = "left",
    legend.direction  = "horizontal",
    legend.justification = "left",
    legend.key.width  = unit(0.7, "lines"),
    legend.key.height = unit(0.7, "lines"),
    legend.spacing.x  = unit(0.25, "lines"),
    legend.spacing.y  = unit(0.15, "lines"),
    legend.text       = element_text(size = 9),
    legend.title      = element_text(size = 10),
    plot.margin       = margin(0, 0, 0, 0)
  )

final





# ================================
#   FIGURE 4 — MULTI-DATASET
# ================================

all_res %>%
  group_by(dataset, method) %>%
  summarise(median_Pearson=median(Pearson)) %>%
  ggplot(aes(dataset, median_Pearson, fill=method)) +
  geom_col(position="dodge") +
  theme_classic(base_size=14) +
  labs(y="Median Pearson")

ggplot(all_res, aes(method, Pearson, fill=method)) +
  geom_boxplot() +
  facet_wrap(~dataset) +
  theme_classic(base_size=14)



# ================================
#   FIGURE 5 — TRANSFER
# ================================

# Panel A — A→B distribution
transfer_tbl %>%
  filter(train_ds=="A", test_ds=="B") %>%
  ggplot(aes(method, Pearson, fill=method)) +
  geom_boxplot() +
  theme_classic(base_size=14) +
  labs(title="A → B transfer")

# Panel B — B→A

transfer_tbl %>%
  filter(train_ds=="B", test_ds=="A") %>%
  ggplot(aes(method, Pearson, fill=method)) +
  geom_boxplot() +
  theme_classic(base_size=14) +
  labs(title="B → A transfer")


# Panel C — Win-rate

transfer_tbl %>%
  group_by(train_ds,test_ds,protein) %>%
  slice_max(Pearson, n=1, with_ties=FALSE) %>%
  count(method) %>%
  ggplot(aes(method,n,fill=method))+
  geom_col()+
  theme_classic(base_size=14)+
  labs(y="# proteins best")

# Panel D — Transfer scatter

plot_transfer <- function(p){
  df <- data.frame(
    truth = y_true_transfer[[p]],
    pred  = y_pred_fmle_transfer[[p]]
  )
  
  ggplot(df,aes(truth,pred))+
    geom_point(alpha=0.5)+
    geom_abline(slope=1,intercept=0,lty=2)+
    theme_classic()+
    labs(title=p)
}

plot_transfer("CD3")
plot_transfer("CD14")

plot_transfer <- function(p){
  df <- data.frame(
    truth = y_true_transfer[[p]],
    pred  = y_pred_fmle_transfer[[p]]
  )
  
  ggplot(df,aes(truth,pred))+
    geom_point(alpha=0.5)+
    geom_abline(slope=1,intercept=0,lty=2)+
    theme_classic()+
    labs(title=p)
}

plot_transfer("CD3")
plot_transfer("CD14")


















#----------------------------------------------------------------------------------#
# Compute average improvement:
#----------------------------------------------------------------------------------#

# FIGURE 3
#----------------------------------------------------------------------------------#
# 3A) Paired scatter: FMLE vs scLinear (R²) with diagonal
#----------------------------------------------------------------------------------#
make_cmp_one_dataset <- function(results_fmle, res_sclinear, res_ctpnet=NULL) {
  # standardize colnames expected
  stopifnot(all(c("protein","R2","Pearson") %in% colnames(results_fmle)))
  stopifnot(all(c("protein","R2","Pearson") %in% colnames(res_sclinear)))
  
  cmp <- inner_join(
    res_sclinear %>% select(protein, R2_scLinear = R2, Pearson_scLinear = Pearson),
    results_fmle %>% select(protein, R2_FMLE = R2, Pearson_FMLE = Pearson),
    by = "protein"
  )
  
  if (!is.null(res_ctpnet)) {
    cmp <- inner_join(
      cmp,
      res_ctpnet %>% select(protein, R2_cTPnet = R2, Pearson_cTPnet = Pearson),
      by = "protein"
    )
  }
  
  cmp
}

plot_paired_scatter <- function(cmp, xcol, ycol, title, label_top_k=5) {
  stopifnot(all(c("protein", xcol, ycol) %in% colnames(cmp)))
  
  delta <- cmp[[ycol]] - cmp[[xcol]]
  top_idx <- order(delta, decreasing = TRUE)[seq_len(min(label_top_k, nrow(cmp)))]
  
  ggplot(cmp, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(size = 3, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggrepel::geom_text_repel(
      data = cmp[top_idx, ],
      aes(label = protein),
      size = 3.5,
      max.overlaps = Inf
    ) +
    theme_classic(base_size = 14) +
    labs(x = xcol, y = ycol, title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}

# Example usage:
cmp <- make_cmp_one_dataset(results_fmle, res_sclinear, res_ctpnet)
p3A <- plot_paired_scatter(cmp, "R2_scLinear", "R2_FMLE", "Protein-wise R²: FMLE vs sclinear")
print(p3A)

ggsave(
  filename = file.path(out_dir, "scatter_FMLE_sclinear.pdf"),
  plot = p3A,
  device = "pdf",
  width = 8,
  height = 6,
  units = "in"
)


mean(cmp$R2_FMLE - cmp$R2_scLinear)
mean(cmp$R2_FMLE - cmp$R2_cTPnet)
mean(cmp$Pearson_FMLE - cmp$Pearson_scLinear)
mean(cmp$Pearson_FMLE - cmp$Pearson_cTPnet)

#----------------------------------------------------------------------------------#
# Per-protein gains
#----------------------------------------------------------------------------------#

cmp$delta_R2_scLinear <- cmp$R2_FMLE - cmp$R2_scLinear
cmp$delta_R2_cTPnet <- cmp$R2_FMLE - cmp$R2_cTPnet

summary(cmp$delta_R2_scLinear)
summary(cmp$delta_R2_cTPnet)

#----------------------------------------------------------------------------------#
# Reviewer-proof statistics
#----------------------------------------------------------------------------------#
wilcox.test(cmp$R2_FMLE, cmp$R2_scLinear, paired=TRUE)
wilcox.test(cmp$R2_FMLE, cmp$R2_cTPnet, paired=TRUE)

wilcox.test(cmp$Pearson_FMLE, cmp$Pearson_scLinear, paired=TRUE)
wilcox.test(cmp$Pearson_FMLE, cmp$Pearson_cTPnet, paired=TRUE)


#----------------------------------------------------------------------------------#
# 3B) ΔR² bar plot (FMLE − scLinear), sorted, horizontal
#----------------------------------------------------------------------------------#
plot_delta_bar <- function(cmp, x_delta = "delta_R2", title = expression(Delta*R^2)) {
  stopifnot(all(c("protein", x_delta) %in% colnames(cmp)))
  
  cmp2 <- cmp %>%
    arrange(.data[[x_delta]]) %>%
    mutate(protein = factor(protein, levels = protein))
  
  ggplot(cmp2, aes(x = protein, y = .data[[x_delta]])) +
    geom_col() +
    coord_flip() +
    theme_classic(base_size = 14) +
    labs(x = "", y = x_delta, title = title)+
    theme(plot.title = element_text(hjust = 0.5))
}

cmp$delta_R2 <- cmp$R2_FMLE - cmp$R2_cTPnet
p3B <- plot_delta_bar(cmp, "delta_R2", title=expression(Delta*R^2~"(FMLE - cTPnet)"))
print(p3B)

ggsave(
  filename = file.path(out_dir, "bar_FMLE_cTPnet.pdf"),
  plot = p3B,
  device = "pdf",
  width = 8,
  height = 6,
  units = "in"
)


#----------------------------------------------------------------------------------#
# Mean ± SE table
#----------------------------------------------------------------------------------#
summary_tbl <- tibble(
  Method=c("FMLE","scLinear","cTPnet"),
  
  Mean_R2=c(
    mean(cmp$R2_fmle),
    mean(cmp$R2_scl),
    mean(cmp$R2_ctp)
  ),
  
  SE_R2=c(
    sd(cmp$R2_fmle)/sqrt(nrow(cmp)),
    sd(cmp$R2_scl)/sqrt(nrow(cmp)),
    sd(cmp$R2_ctp)/sqrt(nrow(cmp))
  ),
  
  Mean_Pearson=c(
    mean(cmp$Pearson_fmle),
    mean(cmp$Pearson_scl),
    mean(cmp$Pearson_ctp)
  ),
  
  SE_Pearson=c(
    sd(cmp$Pearson_fmle)/sqrt(nrow(cmp)),
    sd(cmp$Pearson_scl)/sqrt(nrow(cmp)),
    sd(cmp$Pearson_ctp)/sqrt(nrow(cmp))
  )
)

summary_tbl

mean(cmp$R2_fmle > cmp$R2_scl)
mean(cmp$R2_fmle > cmp$R2_ctp)

# Effect size
median(cmp$R2_fmle - cmp$R2_scl)
median(cmp$R2_fmle - cmp$R2_ctp)

# win rate
mean(cmp$R2_fmle > cmp$R2_scl)
mean(cmp$R2_fmle > cmp$R2_ctp)

# high-impact metric and This shows fraction of proteins with meaningful gains.
mean(cmp$delta_R2_scl > 0.03)
mean(cmp$delta_R2_ctp > 0.02)

#----------------------------------------------------------------------------------#
# FIGURE 4 — Cross-dataset generalization (3 datasets)
#----------------------------------------------------------------------------------#

# 4A) Mean ± SE table + plot

stack_dataset <- function(dataset_name, fmle, scl, ctp=NULL) {
  fmle2 <- fmle %>% transmute(dataset=dataset_name, method="FMLE", protein, R2, Pearson)
  scl2  <- scl  %>% transmute(dataset=dataset_name, method="scLinear", protein, R2, Pearson)
  out <- bind_rows(fmle2, scl2)
  
  if (!is.null(ctp)) {
    ctp2 <- ctp %>% transmute(dataset=dataset_name, method="cTPnet", protein, R2, Pearson)
    out <- bind_rows(out, ctp2)
  }
  out
}

summarize_mean_se <- function(df_long) {
  df_long %>%
    group_by(dataset, method) %>%
    summarise(
      Mean_R2 = mean(R2, na.rm=TRUE),
      SE_R2   = sd(R2, na.rm=TRUE)/sqrt(sum(is.finite(R2))),
      Mean_Pearson = mean(Pearson, na.rm=TRUE),
      SE_Pearson   = sd(Pearson, na.rm=TRUE)/sqrt(sum(is.finite(Pearson))),
      n_prot = n(),
      .groups="drop"
    )
}

plot_mean_se <- function(sum_tbl, metric = c("R2","Pearson")) {
  metric <- match.arg(metric)
  
  if (metric == "R2") {
    ggplot(sum_tbl, aes(x = method, y = Mean_R2)) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin = Mean_R2 - SE_R2, ymax = Mean_R2 + SE_R2), width=0.15) +
      facet_wrap(~dataset, nrow=1) +
      theme_classic(base_size=14) +
      labs(x="", y="Mean R² ± SE", title="Cross-dataset performance (R²)")
  } else {
    ggplot(sum_tbl, aes(x = method, y = Mean_Pearson)) +
      geom_point(size=3) +
      geom_errorbar(aes(ymin = Mean_Pearson - SE_Pearson, ymax = Mean_Pearson + SE_Pearson), width=0.15) +
      facet_wrap(~dataset, nrow=1) +
      theme_classic(base_size=14) +
      labs(x="", y="Mean Pearson ± SE", title="Cross-dataset performance (Pearson)")
  }
}

# ds1_long <- stack_dataset("PBMC-25", results_fmle_ds1, res_sclinear_ds1, res_ctpnet_ds1)
# ds2_long <- stack_dataset("PBMC-14", results_fmle_ds2, res_sclinear_ds2, res_ctpnet_ds2)
# ds3_long <- stack_dataset("PBMC-45", results_fmle_ds3, res_sclinear_ds3, res_ctpnet_ds3)

# all_long <- bind_rows(ds1_long, ds2_long, ds3_long)
# sum_tbl <- summarize_mean_se(all_long)
# print(sum_tbl)
# print(plot_mean_se(sum_tbl, "R2"))
# print(plot_mean_se(sum_tbl, "Pearson"))

# 4B) Win-rate per dataset

win_rate <- function(results_fmle, res_sclinear, res_ctpnet=NULL) {
  cmp <- make_cmp_one_dataset(results_fmle, res_sclinear, res_ctpnet)
  
  out <- list(
    win_vs_scl = mean(cmp$R2_fmle > cmp$R2_scl, na.rm=TRUE),
    median_delta_vs_scl = median(cmp$R2_fmle - cmp$R2_scl, na.rm=TRUE),
    mean_delta_vs_scl   = mean(cmp$R2_fmle - cmp$R2_scl, na.rm=TRUE)
  )
  
  if (!is.null(res_ctpnet)) {
    out$win_vs_ctp <- mean(cmp$R2_fmle > cmp$R2_ctp, na.rm=TRUE)
    out$median_delta_vs_ctp <- median(cmp$R2_fmle - cmp$R2_ctp, na.rm=TRUE)
    out$mean_delta_vs_ctp   <- mean(cmp$R2_fmle - cmp$R2_ctp, na.rm=TRUE)
  }
  as.data.frame(out)
}

# 4C) Reviewer-proof paired statistics per dataset
paired_tests <- function(results_fmle, res_sclinear, res_ctpnet=NULL) {
  cmp <- make_cmp_one_dataset(results_fmle, res_sclinear, res_ctpnet)
  
  p_r2_scl  <- wilcox.test(cmp$R2_fmle, cmp$R2_scl, paired=TRUE)$p.value
  p_pr_scl  <- wilcox.test(cmp$Pearson_fmle, cmp$Pearson_scl, paired=TRUE)$p.value
  
  out <- tibble(
    test = c("R2_FMLE_vs_scLinear","Pearson_FMLE_vs_scLinear"),
    p_value = c(p_r2_scl, p_pr_scl)
  )
  
  if (!is.null(res_ctpnet)) {
    p_r2_ctp <- wilcox.test(cmp$R2_fmle, cmp$R2_ctp, paired=TRUE)$p.value
    p_pr_ctp <- wilcox.test(cmp$Pearson_fmle, cmp$Pearson_ctp, paired=TRUE)$p.value
    out <- bind_rows(out, tibble(
      test = c("R2_FMLE_vs_cTPnet","Pearson_FMLE_vs_cTPnet"),
      p_value = c(p_r2_ctp, p_pr_ctp)
    ))
  }
  
  out %>% mutate(p_adj = p.adjust(p_value, method="BH"))
}


#----------------------------------------------------------------------------------#
# FIGURE 5 — Biological utility
#----------------------------------------------------------------------------------#

# 5A) Add predicted proteins as a new assay in Seurat

add_pred_assay <- function(seu, pred_mat, assay_name) {
  # pred_mat: proteins x cells
  common_cells <- intersect(colnames(seu), colnames(pred_mat))
  common_prot  <- rownames(pred_mat)
  
  stopifnot(length(common_cells) > 100)
  
  pred_mat2 <- pred_mat[, common_cells, drop=FALSE]
  
  seu[[assay_name]] <- CreateAssayObject(data = pred_mat2)  # put in "data" slot
  seu
}
# seu <- add_pred_assay(seu, pred_fmle_mat, "FMLE_PRED")
# seu <- add_pred_assay(seu, pred_scl_mat,  "SCL_PRED")
# seu <- add_pred_assay(seu, pred_ctp_mat,  "CTP_PRED")

# 4B) UMAP overlays (FeaturePlot) for selected proteins

plot_umap_proteins <- function(seu, proteins, assays = c("ADT","FMLE_PRED","SCL_PRED","CTP_PRED")) {
  plots <- list()
  for (a in assays) {
    if (!a %in% names(seu@assays)) next
    for (p in proteins) {
      if (!p %in% rownames(seu[[a]])) next
      plots[[paste(a,p,sep="__")]] <- FeaturePlot(seu, features = p, reduction="umap", assay=a) +
        ggtitle(paste0(a, ": ", p))
    }
  }
  plots
}

# Example:
# prots_show <- c("CD3", "CD4", "CD8a", "HLA.DR")
# umap_plots <- plot_umap_proteins(seu, prots_show)
# umap_plots[["ADT__CD3"]]
# umap_plots[["FMLE_PRED__CD3"]]


# 4C) Cell type separation improvement (quantitative panel)

# Simple: how well protein values separate annotated cell types (one-way ANOVA R²)
anova_r2 <- function(values, labels) {
  df <- data.frame(v=values, lab=as.factor(labels))
  df <- df[is.finite(df$v) & !is.na(df$lab), , drop=FALSE]
  if (nrow(df) < 50 || nlevels(df$lab) < 2) return(NA_real_)
  fit <- lm(v ~ lab, data=df)
  summary(fit)$r.squared
}

compare_separation <- function(seu, protein, label_col="seurat_clusters",
                               assays=c("ADT","FMLE_PRED","SCL_PRED","CTP_PRED")) {
  labs <- seu[[label_col]][,1]
  out <- lapply(assays, function(a){
    if (!a %in% names(seu@assays)) return(NULL)
    if (!protein %in% rownames(seu[[a]])) return(NULL)
    v <- as.numeric(GetAssayData(seu, assay=a, slot="data")[protein, ])
    data.frame(protein=protein, assay=a, anova_R2=anova_r2(v, labs))
  })
  bind_rows(out)
}

# Example:
# sep_tbl <- bind_rows(lapply(prots_show, function(p) compare_separation(seu, p, "seurat_clusters")))
# print(sep_tbl)





