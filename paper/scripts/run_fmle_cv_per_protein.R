
suppressPackageStartupMessages({
  library(data.table)   # fread
  library(Matrix)       # sparse matrices
  library(Seurat)       # Seurat v5 objects + PCA/UMAP
  library(FMLE)         # fmle_train / fmle_cv_parallel / helpers
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)        # write_csv
  library(tibble)
  library(ggplot2)      # plots (ElbowPlot etc. is Seurat, but ggplot2 ok)
  library(glue)         # glue("final_{prot}.rds")
  library(glmnet)       # if FMLE uses it internally or you call it directly
  library(future)
  library(future.apply)
})


bench <- 1  # set 1 or 2 or 3
# base  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_c_doror_%d", bench))
# base  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_%d", bench))
base  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_trasnfer_%d", bench))
# base  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_spatial_%d", bench))
ds    <- "citeseq_v1"
dir.create(base, recursive = TRUE, showWarnings = FALSE)


train_cells <- readRDS(file.path(base, ds, "train_cells.rds"))
test_cells  <- readRDS(file.path(base, ds, "test_cells.rds"))
X           <- readRDS(file.path(base, ds, "X.rds"))
Z           <- readRDS(file.path(base, ds, "Z.rds"))
adt_mat     <- readRDS(file.path(base, ds, "adt_mat.rds"))
if (bench == 1) {
  groups    <- readRDS(file.path(base, ds, "groups.rds"))
} else {
  groups <- NULL
}
 

stopifnot(length(intersect(train_cells, test_cells)) == 0)

stopifnot(all(train_cells %in% rownames(X)))
stopifnot(all(test_cells  %in% rownames(X)))

stopifnot(identical(rownames(X), rownames(Z)))

stopifnot(all(train_cells %in% colnames(adt_mat)))
stopifnot(all(test_cells  %in% colnames(adt_mat)))

## if group is present heer in case of kaggle uncommentt the group

if (bench == 1) {
  stopifnot(exists("groups"))
  stopifnot(!is.null(names(groups)))
  stopifnot(all(train_cells %in% names(groups)))
  stopifnot(all(test_cells  %in% names(groups)))
}

out_dir <- file.path(base, "FMLE", "citeseq_v1_cv")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plan(multicore, workers = 16)   # 12–14 is good on your 16-core (no SMT) system
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
options(future.globals.maxSize = 16 * 1024^3)  # bump if X/Z are big
set.seed(1)

R_grid      <- c(2, 3, 4, 5, 6)
m_grid      <- c(1.5, 1.6, 1.7, 1.8)
lambda_grid <- c(0, 1e-3, 1e-2, 3e-2, 5e-2)
X_cv  <- X[train_cells, , drop = FALSE]
Z_cv  <- Z[train_cells, , drop = FALSE]

if (bench == 1) {
  groups_cv <- groups[train_cells]
} else {
  groups_cv <- NULL
}

best_rows <- vector("list", length(rownames(adt_mat)))
names(best_rows) <- rownames(adt_mat)

for (prot in rownames(adt_mat)) {
  message("=== ", prot, " ===")
  y_raw <- as.numeric(adt_mat[prot, train_cells])
  
  res <- fmle_cv_parallel(
    X = X_cv, y = y_raw, Z = Z_cv,
    R_grid = R_grid, m_grid = m_grid, lambda_grid = lambda_grid,
    folds = 5, groups = groups_cv,                                # groups = groups_cv, # if group is present
    ridge = 1e-6, standardize = TRUE, 
    seed = 1, q = 0.995,
    exec = "future", verbose = TRUE
  )
  
  cv_tbl <- as.data.frame(res$table)
  cv_tbl$protein <- prot
  readr::write_csv(cv_tbl, file.path(out_dir, paste0("cv_", prot, ".csv")))
  
  
  best_rows[[prot]] <- data.frame(
    protein = prot,
    R       = res$best$R,
    m       = res$best$m,
    lambda  = res$best$lambda,
    mse_cv  = res$best$mse,
    folds   = 5,
    seed    = 1,
    q       = 0.995,
    standardize = TRUE,
    ridge   = 1e-6
  )
  readr::write_csv(
    dplyr::bind_rows(best_rows),
    file.path(out_dir, "cv_best_all_adts_running.csv")
  )
}

cv_best_all <- dplyr::bind_rows(best_rows)
readr::write_csv(cv_best_all, file.path(out_dir, "cv_best_all_adts.csv"))





writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_root, "sessionInfo.txt"))
saveRDS(cfg, file.path(cfg$out_root, "config_used.rds"))



















