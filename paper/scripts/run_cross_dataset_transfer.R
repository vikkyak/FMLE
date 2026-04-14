# ====================================================================
# Cross-dataset transfer (A -> B, B -> A) for FMLE / scLinear / cTPnet
# HARD VERIFICATION + FIXED ALIGNMENT + FAIR CALIBRATION
# ====================================================================

suppressPackageStartupMessages({
  library(FMLE)
  library(Seurat) 
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

# ============================
# 0) Utilities: strict checks
# ============================
assert_all <- function(cond, msg) if (!isTRUE(cond)) stop(msg, call. = FALSE)

check_ds <- function(D, tag="DS") {
  assert_all(is.character(D$train_cells) && length(D$train_cells) > 0, paste0(tag, ": train_cells missing"))
  assert_all(is.character(D$test_cells)  && length(D$test_cells)  > 0, paste0(tag, ": test_cells missing"))
  assert_all(is.matrix(D$X) || inherits(D$X, "Matrix"), paste0(tag, ": X not matrix/Matrix"))
  assert_all(is.matrix(D$Z) || inherits(D$Z, "Matrix"), paste0(tag, ": Z not matrix/Matrix"))
  assert_all(is.matrix(D$adt_mat) || inherits(D$adt_mat, "Matrix"), paste0(tag, ": adt_mat not matrix/Matrix"))
  
  # X/Z rownames must be cell IDs
  assert_all(!is.null(rownames(D$X)) && !is.null(rownames(D$Z)), paste0(tag, ": X/Z need rownames=cells"))
  assert_all(identical(rownames(D$X), rownames(D$Z)), paste0(tag, ": rownames(X) != rownames(Z)"))
  
  # train/test must exist in X/Z
  assert_all(all(D$train_cells %in% rownames(D$X)), paste0(tag, ": some train_cells not in rownames(X)"))
  assert_all(all(D$test_cells  %in% rownames(D$X)), paste0(tag, ": some test_cells not in rownames(X)"))
  assert_all(length(intersect(D$train_cells, D$test_cells)) == 0, paste0(tag, ": train/test overlap"))
  
  # adt_mat columns must be cells
  assert_all(!is.null(colnames(D$adt_mat)), paste0(tag, ": adt_mat needs colnames=cells"))
  assert_all(all(D$train_cells %in% colnames(D$adt_mat)), paste0(tag, ": some train_cells not in colnames(adt_mat)"))
  assert_all(all(D$test_cells  %in% colnames(D$adt_mat)), paste0(tag, ": some test_cells not in colnames(adt_mat)"))
  
  invisible(TRUE)
}

# ============================
# 1) Load dataset from folder
# ============================
load_ds <- function(ds_dir) {
  D <- list(
    train_cells = readRDS(file.path(ds_dir, "train_cells.rds")),
    test_cells  = readRDS(file.path(ds_dir, "test_cells.rds")),
    X           = readRDS(file.path(ds_dir, "X.rds")),
    Z           = readRDS(file.path(ds_dir, "Z.rds")),
    adt_mat     = readRDS(file.path(ds_dir, "adt_mat.rds")),
    seu_prep    = readRDS(file.path(ds_dir, "seu_final.rds")),
    groups      = if (file.exists(file.path(ds_dir, "groups.rds"))) readRDS(file.path(ds_dir, "groups.rds")) else NULL
  )
  
  D$X <- as.matrix(D$X)
  D$Z <- as.matrix(D$Z)
  D$adt_mat <- as.matrix(D$adt_mat)
  
  check_ds(D, tag=basename(ds_dir))
  D
}

# =========================================
# 2) ONE global protein canonicalization
# =========================================
canon_prot_global <- function(x){
  x <- as.character(x)
  x <- sub("\\.TotalSeqB$", "", x)
  x <- sub("-TotalSeqB$", "", x)
  x <- gsub("\\.", "-", x)  # CD127.IL7Ra -> CD127-IL7Ra, PD.1 -> PD-1, CCR7 dots -> hyphens
  x <- gsub("^CD127-IL7Ra$|^CD127$|^CD127-IL7Ra$", "CD127-IL7Ra", x)
  x[x == "CD127"] <- "CD127-IL7Ra"
  x
}

# ===============================
# 3) Metrics (same numeric scale)
# ===============================
calc_metrics <- function(yt, yp) {
  ok <- is.finite(yt) & is.finite(yp)
  yt <- yt[ok]; yp <- yp[ok]
  mse <- mean((yt-yp)^2)
  den <- sum((yt-mean(yt))^2)
  tibble(
    MSE      = mse,
    RMSE     = sqrt(mse),
    MAE      = mean(abs(yt-yp)),
    Pearson  = suppressWarnings(cor(yt, yp)),
    Spearman = suppressWarnings(cor(yt, yp, method="spearman")),
    R2       = if (den > 0) 1 - sum((yt-yp)^2)/den else NA_real_
  )
}

# ===================================================
# 4) Shared Z from A-train PCA basis (transfer-valid)
# ===================================================
build_shared_Z_from_Atrain <- function(DS, genes_common, k = 20, seed = 1) {
  stopifnot(all(c("A","B") %in% names(DS)))
  stopifnot(!is.null(DS$A$train_cells), !is.null(DS$B$test_cells))
  
  XA <- DS$A$X[, genes_common, drop = FALSE]
  XB <- DS$B$X[, genes_common, drop = FALSE]
  
  XA_tr <- XA[DS$A$train_cells, , drop = FALSE]
  mu   <- colMeans(XA_tr)
  sdev <- apply(XA_tr, 2, stats::sd)
  sdev[sdev == 0] <- 1
  
  scale_A <- function(M) sweep(sweep(M, 2, mu, "-"), 2, sdev, "/")
  
  XA_tr_s <- scale_A(XA_tr)
  
  suppressPackageStartupMessages(library(irlba))
  set.seed(seed)
  pc <- irlba::prcomp_irlba(XA_tr_s, n = k, center = FALSE, scale. = FALSE)
  
  ZA <- scale_A(XA) %*% pc$rotation
  ZB <- scale_A(XB) %*% pc$rotation
  
  colnames(ZA) <- paste0("PC", seq_len(k))
  colnames(ZB) <- paste0("PC", seq_len(k))
  
  DS$A$Z <- ZA
  DS$B$Z <- ZB
  DS
}

# =====================================================
# 5) FMLE transfer evaluation (with affine calibration)
# =====================================================
eval_transfer_fmle_cal <- function(A, B, best_tbl_A, q = 0.995,
                                   zscore_Z = TRUE, seed = 1) {
  
  prot_common <- intersect(rownames(A$adt_mat), rownames(B$adt_mat))
  best_tbl_A  <- best_tbl_A %>% dplyr::filter(protein %in% prot_common)
  
  Xtr <- A$X[A$train_cells, , drop = FALSE]
  Xte <- B$X[B$test_cells,  , drop = FALSE]
  
  Ztr <- A$Z[A$train_cells, , drop = FALSE]
  Zte <- B$Z[B$test_cells,  , drop = FALSE]
  if (is.null(Ztr) || is.null(Zte)) stop("A$Z or B$Z is NULL. Build shared Z first.")
  
  if (zscore_Z) {
    zmu <- colMeans(Ztr)
    zsd <- apply(Ztr, 2, stats::sd); zsd[zsd == 0] <- 1
    Ztr <- sweep(sweep(Ztr, 2, zmu, "-"), 2, zsd, "/")
    Zte <- sweep(sweep(Zte, 2, zmu, "-"), 2, zsd, "/")
  }
  
  out <- lapply(best_tbl_A$protein, function(p){
    
    cfg <- best_tbl_A[best_tbl_A$protein == p, ][1, ]
    if (nrow(cfg) == 0) return(NULL)
    
    yA_train_raw <- as.numeric(A$adt_mat[p, A$train_cells])
    tf <- FMLE:::cap_and_scale_fit(yA_train_raw, q = q)
    yA_train_scaled <- FMLE:::cap_and_scale_apply(yA_train_raw, tf)
    
    fit <- FMLE::fmle_train(
      X = Xtr, y = yA_train_scaled, Z = Ztr,
      R = cfg$R, m = cfg$m, lambda_l1 = cfg$lambda,
      ridge = 1e-6, standardize = TRUE, seed = seed
    )
    
    pr_tr <- FMLE::fmle_predict(fit, X_new = Xtr, Z_new = Ztr)
    yhatA <- as.numeric(pr_tr$mean)
    
    ok <- is.finite(yA_train_scaled) & is.finite(yhatA)
    
    yB_test_raw <- as.numeric(B$adt_mat[p, B$test_cells])
    yt <- FMLE:::cap_and_scale_apply(yB_test_raw, tf)
    
    pr_te <- FMLE::fmle_predict(fit, X_new = Xte, Z_new = Zte)
    yhatB <- as.numeric(pr_te$mean)
    
    if (sum(ok) < 20 || length(unique(yhatA[ok])) < 5) {
      yp <- yhatB
    } else {
      cal <- lm(yA_train_scaled[ok] ~ yhatA[ok])
      cc <- coef(cal)
      yp <- if (any(!is.finite(cc))) yhatB else cc[1] + cc[2] * yhatB
    }
    
    cbind(data.frame(protein = p, method = "FMLE"), calc_metrics(yt, yp))
  })
  
  dplyr::bind_rows(out)
}

# ==========================================================
# 6) Generic calibrated transfer eval for scLinear / cTPnet
# ==========================================================
eval_transfer_calibrated <- function(A, B, yhatA_train_raw, yhatB_test_raw,
                                     method_name, q=0.995) {
  
  assert_all(all(A$train_cells %in% colnames(yhatA_train_raw)),
             paste0(method_name, ": missing some A train cells in yhatA_train_raw"))
  assert_all(all(B$test_cells %in% colnames(yhatB_test_raw)),
             paste0(method_name, ": missing some B test cells in yhatB_test_raw"))
  
  yhatA_train_raw <- yhatA_train_raw[, A$train_cells, drop=FALSE]
  yhatB_test_raw  <- yhatB_test_raw[,  B$test_cells,  drop=FALSE]
  
  prot_common <- Reduce(intersect, list(
    rownames(A$adt_mat), rownames(B$adt_mat),
    rownames(yhatA_train_raw), rownames(yhatB_test_raw)
  ))
  assert_all(length(prot_common) > 0, paste0(method_name, ": no common proteins after intersect"))
  
  out <- lapply(prot_common, function(p){
    yA_raw <- as.numeric(A$adt_mat[p, A$train_cells])
    tf <- FMLE:::cap_and_scale_fit(yA_raw, q=q)
    yA_scaled <- FMLE:::cap_and_scale_apply(yA_raw, tf)
    
    yhatA <- as.numeric(yhatA_train_raw[p, A$train_cells])
    ok <- is.finite(yA_scaled) & is.finite(yhatA)
    
    yB_raw <- as.numeric(B$adt_mat[p, B$test_cells])
    yB_scaled <- FMLE:::cap_and_scale_apply(yB_raw, tf)
    
    yhatB <- as.numeric(yhatB_test_raw[p, B$test_cells])
    
    if (sum(ok) < 10 || length(unique(yhatA[ok])) < 5) {
      yB_hat <- yhatB
    } else {
      cal <- lm(yA_scaled[ok] ~ yhatA[ok])
      cc <- coef(cal)
      yB_hat <- if (any(!is.finite(cc))) yhatB else cc[1] + cc[2] * yhatB
    }
    
    tibble(protein=p, method=method_name) %>% bind_cols(calc_metrics(yB_scaled, yB_hat))
  })
  
  bind_rows(out)
}

# =================================================
# 7) RUN: load DS, enforce ONE naming + SAME genes
# =================================================
base_1 <- file.path(cfg$out_root, "benchmarks_trasnfer_1")
base_2 <- file.path(cfg$out_root, "benchmarks_trasnfer_2")
base_3 <- file.path(cfg$out_root, "benchmarks_trasnfer_3")
#  in case of PBMC
base <- file.path(cfg$out_root, "transfer_preds_PBMC")
#  in case of BMMC
base <- file.path(cfg$out_root, "transfer_preds_BMMC")

ds_paths <- list(
  A = file.path(base_1, "citeseq_v1"),
  B = file.path(base_2, "citeseq_v1"),
  C = file.path(base_3, "citeseq_v1")
)

dirs <- c(base_1, base_2, base_3, unlist(ds_paths, use.names = FALSE))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

DS <- lapply(ds_paths, load_ds)

# Apply canonical protein naming ONCE (truth)
rownames(DS$A$adt_mat) <- canon_prot_global(rownames(DS$A$adt_mat))
rownames(DS$B$adt_mat) <- canon_prot_global(rownames(DS$B$adt_mat))
rownames(DS$C$adt_mat) <- canon_prot_global(rownames(DS$C$adt_mat))


DS$A$X_full <- t(as.matrix(GetAssayData(DS$A$seu_prep, assay="RNA", layer="data")))
DS$B$X_full <- t(as.matrix(GetAssayData(DS$B$seu_prep, assay="RNA", layer="data")))
DS$C$X_full <- t(as.matrix(GetAssayData(DS$C$seu_prep, assay="RNA", layer="data")))
# # ================================================================================
# # this is used to get common genes in all datatset for transfer
# genes_A <- colnames(DS$A$X_full)
# genes_B <- colnames(DS$B$X_full)
# genes_C <- colnames(DS$C$X_full)
# genes_ABC <- Reduce(intersect, list(genes_A, genes_B, genes_C))
# length(genes_ABC)
# 
# n_pool <- 4000
# hvgA <- VariableFeatures(FindVariableFeatures(DS$A$seu_prep, assay="RNA", nfeatures=n_pool, verbose=FALSE))
# hvgB <- VariableFeatures(FindVariableFeatures(DS$B$seu_prep, assay="RNA", nfeatures=n_pool, verbose=FALSE))
# hvgC <- VariableFeatures(FindVariableFeatures(DS$C$seu_prep, assay="RNA", nfeatures=n_pool, verbose=FALSE))
# 
# rankA <- setNames(seq_along(hvgA), hvgA)
# pool <- Reduce(union, list(hvgA, hvgB, hvgC))
# pool <- intersect(pool, genes_ABC)
# pool <- pool[order(ifelse(pool %in% names(rankA), rankA[pool], 1e9))]
# 
# genes_use2 <- head(pool, 2000)
# stopifnot(length(genes_use2) == 2000)
# 
# 
# cat("A genes:", length(genes_A),
#     "B genes:", length(genes_B),
#     "C genes:", length(genes_C),
#     "genes_use2:", length(genes_use2), "\n")
# 
# DS$A$X <- DS$A$X_full[, genes_use2, drop=FALSE]
# DS$B$X <- DS$B$X_full[, genes_use2, drop=FALSE]
# DS$C$X <- DS$C$X_full[, genes_use2, drop=FALSE]
# 
# stopifnot(
#   identical(colnames(DS$A$X), genes_use2),
#   identical(colnames(DS$B$X), genes_use2),
#   identical(colnames(DS$C$X), genes_use2)
# )
# 
# genes_common <- genes_use2
# 
# saveRDS(genes_common, file.path("~/Desktop/FMLE/transfer_preds", "gene_panel_2000.rds"))
# 
# write.csv(data.frame(gene = genes_common),
#           "~/Desktop/FMLE/transfer_preds/gene_panel_2000.csv",
#           row.names = FALSE)
# # ================================================================================

genes_common <- readRDS(file.path(base, "gene_panel_2000.rds")) 

# ============================================================
# 8) FMLE (needs best config per direction)
# ============================================================
best_cfg <- list(
  A = readr::read_csv(file.path(base_1, "FMLE", "citeseq_v1_cv", "cv_best_all_adts.csv")),
  B = readr::read_csv(file.path(base_2, "FMLE", "citeseq_v1_cv", "cv_best_all_adts.csv")),
  C = readr::read_csv(file.path(base_3, "FMLE", "citeseq_v1_cv", "cv_best_all_adts.csv")))

req_cols <- c("protein","R","m","lambda")
assert_all(all(req_cols %in% colnames(best_cfg$A)), "best_cfg$A missing columns protein/R/m/lambda")
assert_all(all(req_cols %in% colnames(best_cfg$B)), "best_cfg$B missing columns protein/R/m/lambda")
assert_all(all(req_cols %in% colnames(best_cfg$C)), "best_cfg$C missing columns protein/R/m/lambda")


for (nm in c("A","B","C")) {
  best_cfg[[nm]]$protein <- canon_prot_global(best_cfg[[nm]]$protein)
  best_cfg[[nm]] <- best_cfg[[nm]] %>% dplyr::distinct(protein, .keep_all=TRUE)
}

# Gate
k_c <- 20

DS_AB <- build_shared_Z_from_Atrain(list(A=DS$A, B=DS$B), genes_common, k = k_c , seed=1)
DS_AC <- build_shared_Z_from_Atrain(list(A=DS$A, B=DS$C), genes_common, k = k_c, seed=1)

# Build transfer-valid shared Z (B basis for B->A, B basis for B->C)

DS_BA <- build_shared_Z_from_Atrain(list(A=DS$B, B=DS$A), genes_common, k = k_c, seed=1)
DS_BC <- build_shared_Z_from_Atrain(list(A=DS$B, B=DS$C), genes_common, k = k_c, seed=1)

# Build transfer-valid shared Z (C basis for C->A, C basis for C->B)

DS_CA <- build_shared_Z_from_Atrain(list(A=DS$C, B=DS$A), genes_common, k = k_c, seed=1)
DS_CB <- build_shared_Z_from_Atrain(list(A=DS$C, B=DS$B), genes_common, k = k_c, seed=1)

stopifnot(identical(rownames(DS_AB$A$X), rownames(DS$A$X)))
stopifnot(identical(rownames(DS_AB$B$X), rownames(DS$B$X)))
stopifnot(identical(rownames(DS_AC$A$X), rownames(DS$A$X)))
stopifnot(identical(rownames(DS_AC$B$X), rownames(DS$C$X)))
stopifnot(identical(rownames(DS_BC$A$X), rownames(DS$B$X)))
stopifnot(identical(rownames(DS_BC$B$X), rownames(DS$C$X)))

stopifnot(identical(rownames(DS_BA$A$X), rownames(DS$B$X)))
stopifnot(identical(rownames(DS_BA$B$X), rownames(DS$A$X)))
stopifnot(identical(rownames(DS_CA$A$X), rownames(DS$C$X)))
stopifnot(identical(rownames(DS_CA$B$X), rownames(DS$A$X)))
stopifnot(identical(rownames(DS_CB$A$X), rownames(DS$C$X)))
stopifnot(identical(rownames(DS_CB$B$X), rownames(DS$B$X)))

stopifnot(nrow(DS_AB$A$Z) == nrow(DS_AB$A$X))
stopifnot(nrow(DS_AB$B$Z) == nrow(DS_AB$B$X))
stopifnot(length(intersect(rownames(DS_AB$A$X), rownames(DS_AB$B$X))) == 0)

stopifnot(identical(colnames(DS$A$X), genes_common))
stopifnot(identical(colnames(DS$B$X), genes_common))
stopifnot(identical(colnames(DS$C$X), genes_common))


res_fmle_AB <- eval_transfer_fmle_cal(DS_AB$A, DS_AB$B, best_cfg$A, seed=1) %>%
  dplyr::mutate(train_ds="A", test_ds="B")
res_fmle_AB <- dplyr::bind_rows(res_fmle_AB)
readr::write_csv(res_fmle_AB, file.path(ctpAB_dir, "res_fmle_AB.csv"))

res_fmle_AC <- eval_transfer_fmle_cal(DS_AC$A, DS_AC$B, best_cfg$A, seed=1) %>%
  dplyr::mutate(train_ds="A", test_ds="C")
res_fmle_AC <- dplyr::bind_rows(res_fmle_AC)
readr::write_csv(res_fmle_AC, file.path(ctpAC_dir, "res_fmle_AC.csv"))

res_fmle_BA <- eval_transfer_fmle_cal(DS_BA$A, DS_BA$B, best_cfg$B, seed=1) %>%
  dplyr::mutate(train_ds="B", test_ds="A")
res_fmle_BA <- dplyr::bind_rows(res_fmle_BA)
readr::write_csv(res_fmle_BA, file.path(ctpBA_dir, "res_fmle_BA.csv"))

res_fmle_BC <- eval_transfer_fmle_cal(DS_BC$A, DS_BC$B, best_cfg$B, seed=1) %>%
  dplyr::mutate(train_ds="B", test_ds="C")
res_fmle_BC <- dplyr::bind_rows(res_fmle_BC)
readr::write_csv(res_fmle_BC, file.path(ctpBC_dir, "res_fmle_BC.csv"))

res_fmle_CA <- eval_transfer_fmle_cal(DS_CA$A, DS_CA$B, best_cfg$C, seed=1) %>%
  dplyr::mutate(train_ds="C", test_ds="A")
res_fmle_CA <- dplyr::bind_rows(res_fmle_CA)
readr::write_csv(res_fmle_CA, file.path(ctpCA_dir, "res_fmle_CA.csv"))

res_fmle_CB <- eval_transfer_fmle_cal(DS_CB$A, DS_CB$B, best_cfg$C, seed=1) %>%
  dplyr::mutate(train_ds="C", test_ds="B")
res_fmle_CB <- dplyr::bind_rows(res_fmle_CB)
readr::write_csv(res_fmle_CB, file.path(ctpCB_dir, "res_fmle_CB.csv"))

# =======================================================================
# 9) scLinear (train on source train, predict source train + target test)
#    IMPORTANT: subset RNA genes to genes_common for fairness.
# =======================================================================
get_sclinear_transfer_preds <- function(seu_src, seu_tgt, train_cells_src, test_cells_tgt, genes_common) {
  
  assert_all(all(train_cells_src %in% colnames(seu_src)), "scLinear: some train_cells not in source seu")
  assert_all(all(test_cells_tgt  %in% colnames(seu_tgt)), "scLinear: some test_cells not in target seu")
  
  seu_src_train <- seu_src[, train_cells_src, drop=FALSE]
  pipe <- scLinear::create_adt_predictor()
  
  gexp_src_train <- Seurat::GetAssayData(seu_src_train, assay="RNA", layer="counts")
  adt_src_train  <- Seurat::GetAssayData(seu_src_train, assay="ADT", layer="counts")
  
  # fairness: restrict genes
  gexp_src_train <- gexp_src_train[genes_common, , drop=FALSE]
  
  pipe <- scLinear::fit_predictor(
    pipe          = pipe,
    gexp_train    = as.matrix(gexp_src_train),
    adt_train     = as.matrix(adt_src_train),
    normalize_gex = TRUE,
    normalize_adt = FALSE
  )
  
  # predict on src-train
  pred_src_train <- scLinear::adt_predict(pipe=pipe, gexp=gexp_src_train, normalize=TRUE)
  pred_src_train_raw <- as.matrix(pred_src_train@data) # proteins x cells
  
  # predict on tgt-test
  seu_tgt_test <- seu_tgt[, test_cells_tgt, drop=FALSE]
  gexp_tgt_test <- Seurat::GetAssayData(seu_tgt_test, assay="RNA", layer="counts")
  gexp_tgt_test <- gexp_tgt_test[genes_common, , drop=FALSE]
  
  pred_tgt_test <- scLinear::adt_predict(pipe=pipe, gexp=gexp_tgt_test, normalize=TRUE)
  pred_tgt_test_raw <- as.matrix(pred_tgt_test@data)
  
  list(
    pred_src_train_raw = pred_src_train_raw,
    pred_tgt_test_raw  = pred_tgt_test_raw
  )
}

run_scl_transfer <- function(src, tgt, src_tag, tgt_tag){
  
  tmp <- get_sclinear_transfer_preds(
    src$seu_prep,
    tgt$seu_prep,
    src$train_cells,
    tgt$test_cells,
    genes_common)
  
  pred_train <- tmp$pred_src_train_raw
  pred_test  <- tmp$pred_tgt_test_raw
  
  rownames(pred_train) <- canon_prot_global(rownames(pred_train))
  rownames(pred_test)  <- canon_prot_global(rownames(pred_test))
  
  eval_transfer_calibrated(src, tgt,
                           pred_train, pred_test,
                           "scLinear") %>%
    dplyr::mutate(train_ds=src_tag, test_ds=tgt_tag)
}

res_scl_AB <- run_scl_transfer(DS$A, DS$B, "A","B")
res_scl_AB <- dplyr::bind_rows(res_scl_AB)
readr::write_csv(res_scl_AB, file.path(ctpAB_dir, "res_scl_AB.csv"))

res_scl_AC <- run_scl_transfer(DS$A, DS$C, "A","C")
res_scl_AC <- dplyr::bind_rows(res_scl_AC)
readr::write_csv(res_scl_AC, file.path(ctpAC_dir, "res_scl_AC.csv"))

res_scl_BA <- run_scl_transfer(DS$B, DS$A, "B","A")
res_scl_BA <- dplyr::bind_rows(res_scl_BA)
readr::write_csv(res_scl_BA, file.path(ctpBA_dir, "res_scl_BA.csv"))

res_scl_BC <- run_scl_transfer(DS$B, DS$C, "B","C")
res_scl_BC <- dplyr::bind_rows(res_scl_BC)
readr::write_csv(res_scl_BC, file.path(ctpBC_dir, "res_scl_BC.csv"))

res_scl_CA <- run_scl_transfer(DS$C, DS$A, "C","A")
res_scl_CA <- dplyr::bind_rows(res_scl_CA)
readr::write_csv(res_scl_CA, file.path(ctpCA_dir, "res_scl_CA.csv"))

res_scl_CB <- run_scl_transfer(DS$C, DS$B, "C","B")
res_scl_CB <- dplyr::bind_rows(res_scl_CB)
readr::write_csv(res_scl_CB, file.path(ctpCB_dir, "res_scl_CB.csv"))

# ============================================================
#  cTPnet 
# ============================================================

transfer_root <- file.path(cfg$out_root, "transfer_preds")

ctpAB_dir <- file.path(transfer_root, "ctp_A_to_B")
ctpAC_dir <- file.path(transfer_root, "ctp_A_to_C")
ctpBA_dir <- file.path(transfer_root, "ctp_B_to_A")
ctpBC_dir <- file.path(transfer_root, "ctp_B_to_C")
ctpCA_dir <- file.path(transfer_root, "ctp_C_to_A")
ctpCB_dir <- file.path(transfer_root, "ctp_C_to_B")

dirs <- c(ctpAB_dir, ctpAC_dir, ctpBA_dir, ctpBC_dir, ctpCA_dir, ctpCB_dir)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

read_ctp <- function(path){
  as.matrix(read.csv(path, row.names=1, check.names=FALSE)) # proteins x cells
}

run_ctp_transfer <- function(src, tgt, dir_path, src_tag, tgt_tag) {
  
  pred_train <- read_ctp(file.path(dir_path,
                                   paste0("ctpnet_pred_train_source", src_tag, ".csv")))
  pred_test  <- read_ctp(file.path(dir_path,
                                   paste0("ctpnet_pred_test_target", tgt_tag, ".csv")))
  
  # canonicalize
  rownames(pred_train) <- canon_prot_global(rownames(pred_train))
  rownames(pred_test)  <- canon_prot_global(rownames(pred_test))
  
  # strict cell alignment
  pred_train <- pred_train[, src$train_cells, drop=FALSE]
  pred_test  <- pred_test[,  tgt$test_cells,  drop=FALSE]
  
  eval_transfer_calibrated(src, tgt,
                           pred_train, pred_test,
                           "cTPnet") %>%
    dplyr::mutate(train_ds=src_tag, test_ds=tgt_tag)
}

res_ctp_AB <- run_ctp_transfer(DS$A, DS$B, ctpAB_dir, "A","B")
res_ctp_AB <- dplyr::bind_rows(res_ctp_AB)
readr::write_csv(res_ctp_AB, file.path(ctpAB_dir, "res_ctp_AB.csv"))

res_ctp_AC <- run_ctp_transfer(DS$A, DS$C, ctpAC_dir, "A","C")
res_ctp_AC <- dplyr::bind_rows(res_ctp_AC)
readr::write_csv(res_ctp_AC, file.path(ctpAC_dir, "res_ctp_AC.csv"))


res_ctp_BA <- run_ctp_transfer(DS$B, DS$A, ctpBA_dir, "B","A")
res_ctp_BA <- dplyr::bind_rows(res_ctp_BA)
readr::write_csv(res_ctp_BA, file.path(ctpBA_dir, "res_ctp_BA.csv"))

res_ctp_BC <- run_ctp_transfer(DS$B, DS$C, ctpBC_dir, "B","C")
res_ctp_BC <- dplyr::bind_rows(res_ctp_BC)
readr::write_csv(res_ctp_BC, file.path(ctpBC_dir, "res_ctp_BC.csv"))

res_ctp_CA <- run_ctp_transfer(DS$C, DS$A, ctpCA_dir, "C","A")
res_ctp_CA <- dplyr::bind_rows(res_ctp_CA)
readr::write_csv(res_ctp_CA, file.path(ctpCA_dir, "res_ctp_CA.csv"))

res_ctp_CB <- run_ctp_transfer(DS$C, DS$B, ctpCB_dir, "C","B")
res_ctp_CB <- dplyr::bind_rows(res_ctp_CB)
readr::write_csv(res_ctp_CB, file.path(ctpCB_dir, "res_ctp_CB.csv"))

# ============================================================
# 11) Combine + summarize
# ============================================================
transfer_tbl <- bind_rows(
  res_fmle_AB, res_fmle_AC, res_fmle_BA, res_fmle_BC, res_fmle_CA, res_fmle_CB,
  res_scl_AB, res_scl_AC, res_scl_BA, res_scl_BC, res_scl_CA, res_scl_CB,
  res_ctp_AB, res_ctp_AC, res_ctp_BA, res_ctp_BC, res_ctp_CA, res_ctp_CB,
) %>%
  relocate(train_ds, test_ds, protein, method)

transfer_summary <- transfer_tbl %>%
  group_by(method, train_ds, test_ds) %>%
  summarise(
    med_Pearson  = median(Pearson, na.rm=TRUE),
    mean_Pearson = mean(Pearson, na.rm=TRUE),
    med_R2       = median(R2, na.rm=TRUE),
    mean_R2      = mean(R2, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(train_ds, test_ds, desc(med_Pearson))

transfer_tbl
transfer_summary


# paired Wilcoxon signed-rank test

df <- transfer_tbl %>%
  dplyr::select(train_ds, test_ds, protein, method, Pearson)

wilcox_transfer <- function(df, tr, te) {
  
  sub <- df %>%
    dplyr::filter(train_ds == tr,
                  test_ds  == te,
                  method %in% c("FMLE","cTPnet","scLinear")) %>%
    dplyr::select(train_ds, test_ds, protein, method, Pearson) %>%
    tidyr::pivot_wider(names_from = method, values_from = Pearson)
  
  # FMLE vs cTPnet (paired)
  sub_fc <- sub %>% dplyr::filter(is.finite(FMLE), is.finite(cTPnet))
  p_fc <- if (nrow(sub_fc) >= 3) {
    wilcox.test(sub_fc$FMLE, sub_fc$cTPnet, paired = TRUE, exact = FALSE)$p.value
  } else NA_real_
  
  # FMLE vs scLinear (paired)
  sub_fs <- sub %>% dplyr::filter(is.finite(FMLE), is.finite(scLinear))
  p_fs <- if (nrow(sub_fs) >= 3) {
    wilcox.test(sub_fs$FMLE, sub_fs$scLinear, paired = TRUE, exact = FALSE)$p.value
  } else NA_real_
  
  data.frame(
    train_ds = tr,
    test_ds  = te,
    n_prot_FMLE_cTP = nrow(sub_fc),
    n_prot_FMLE_scL = nrow(sub_fs),
    median_FMLE = median(sub$FMLE, na.rm = TRUE),
    median_cTP  = median(sub$cTPnet, na.rm = TRUE),
    median_scL  = median(sub$scLinear, na.rm = TRUE),
    p_FMLE_vs_cTP = p_fc,
    p_FMLE_vs_scL = p_fs
  )
}



pairs <- unique(df_all[, c("train_ds","test_ds")])

wilcox_results <- purrr::map2_df(
  pairs$train_ds,
  pairs$test_ds,
  ~ wilcox_transfer(df_all, .x, .y)
)

wilcox_results

df_wide <- df_all %>%
  filter(method %in% c("FMLE","cTPnet","scLinear"  )) %>%
  select(train_ds, test_ds, protein, method, Pearson) %>%
  pivot_wider(names_from = method, values_from = Pearson) %>%
  filter(is.finite(FMLE), is.finite(cTPnet))   # paired complete cases

wilcox.test(df_wide$FMLE, df_wide$cTPnet, paired = TRUE, exact = FALSE)
wilcox.test(df_wide$FMLE, df_wide$scLinear, paired = TRUE, exact = FALSE)


writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_root, "sessionInfo.txt"))
saveRDS(cfg, file.path(cfg$out_root, "config_used.rds"))


