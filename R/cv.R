# -----------------------------
# Internal helper: reuse centers (single task)
# -----------------------------
#' @keywords internal
fmle_train_from_centers <- function(
    X, y, Z, centers, m,
    lambda_l1   = 0,
    ridge       = 1e-6,
    standardize = TRUE,
    sx = NULL, sz = NULL
){
  X <- as.matrix(X)
  y <- as.numeric(y)
  Z <- as.matrix(Z)
  
  if (standardize) {
    stopifnot(!is.null(sx), !is.null(sz))  # must pass same scalers used to compute centers
    Xs <- sweep(sweep(X, 2, sx$mu, "-"), 2, sx$sd, "/")
    Zs <- sweep(sweep(Z, 2, sz$mu, "-"), 2, sz$sd, "/")
  } else {
    Xs <- X
    Zs <- Z
  }
  
  R   <- nrow(centers)
  eps <- 1e-12
  
  # memberships from cached centers (same as in fmle_predict)
  CZ    <- tcrossprod(centers, Zs)        # R x N
  nZ    <- rowSums(Zs^2)                  # N
  nC    <- rowSums(centers^2)             # R
  dist2 <- sweep(sweep(-2 * CZ, 2, nZ, "+"), 1, nC, "+")
  dist2 <- pmax(dist2, 0)
  Dpow  <- (dist2 + eps)^(-2/(m - 1))
  U     <- t(Dpow / matrix(colSums(Dpow) + 1e-18,
                           nrow = R, ncol = nrow(Zs), byrow = TRUE))
  w_all <- U^m
  
  experts <- lapply(seq_len(R), function(r)
    .fit_expert(Xs, y, w_all[, r],
                lambda_l1 = lambda_l1,
                ridge     = ridge)
  )
  
  structure(list(R = R, m = m, centers = centers, experts = experts,
                 ridge = ridge, lambda_l1 = lambda_l1,
                 standardize = standardize, sx = sx, sz = sz, U = U),
            class = "fmle")
}

# -----------------------------
# Single-task CV (parallel-capable)
# -----------------------------

#' Cross-validation for single-task FMLE
#'
#' Performs K-fold cross-validation over grids of R (number of experts),
#' m (fuzzifier), and lambda_l1 (L1 penalty), optionally using grouped
#' folds and future-based parallelization.
#'
#' @inheritParams fmle_train
#' @param X matrix, n x p features.
#' @param y numeric response vector of length n.
#' @param Z matrix, n x d gating features (e.g., PCs).
#' @param R_grid grid of R values.
#' @param m_grid grid of m values.
#' @param lambda_grid grid of lambda_l1 values.
#' @param folds number of CV folds.
#' @param groups optional factor of length n for grouped CV; NULL gives random folds.
#' @param exec character, "sequential" or "future" (for future_lapply).
#' @param verbose logical; print progress messages.
#' @param q Quantile used to cap the response before log1p/scale preprocessing in cross-validation.
#'
#' @return A list with components:
#' \describe{
#'   \item{best}{List with the best (R, m, lambda, mse).}
#'   \item{table}{Data frame of all configs and their CV MSE.}
#' }
#' @export
fmle_cv_parallel <- function(
    X, y, Z,
    R_grid      = c(4, 6, 8),
    m_grid      = c(1.6, 1.8, 2.0),
    lambda_grid = c(0, 1e-3, 1e-2, 3e-2),
    folds       = 5,
    groups      = NULL,
    ridge       = 1e-6,
    standardize = TRUE,
    seed        = 1,
    exec        = c("sequential", "future"),
    verbose     = TRUE,
    q           = 0.995 
){
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  y <- as.numeric(y)
  y_raw <- y
  stopifnot(nrow(X) == nrow(Z), length(y) == nrow(X))
  exec <- match.arg(exec)
  
  set.seed(seed)
  N <- nrow(X)
  
  # --- build validation indices per fold ---
  if (!is.null(groups)) {
    stopifnot(length(groups) == N)
    val_idx_list <- .make_grouped_folds(groups, K = folds, seed = seed)
  } else {
    fold_ids     <- sample(rep(seq_len(folds), length.out = N))
    val_idx_list <- lapply(seq_len(folds), function(f) which(fold_ids == f))
  }
  
  # --- grid over (R, m) only, to cache FCM once per (R,m,fold) ---
  RM_grid <- expand.grid(R = R_grid, m = m_grid, KEEP.OUT.ATTRS = FALSE)
  n_RM    <- nrow(RM_grid)
  fcm_cache <- vector("list", n_RM)
  
  if (verbose) message("Precomputing FCM (gates) per (R, m) and fold ...")
  
  for (k in seq_len(n_RM)) {
    Rv <- RM_grid$R[k]
    mv <- RM_grid$m[k]
    
    fcm_cache[[k]] <- vector("list", length(val_idx_list))
    
    for (f in seq_along(val_idx_list)) {
      va <- val_idx_list[[f]]
      tr <- base::setdiff(seq_len(N), va)
      
      Xtr <- X[tr, , drop = FALSE]
      Ztr <- Z[tr, , drop = FALSE]
      tf  <- cap_and_scale_fit(y_raw[tr], q = q)
      ytr <- cap_and_scale_apply(y_raw[tr], tf)
      
      if (standardize) {
        sx <- .scaler_fit(Xtr)
        sz <- .scaler_fit(Ztr)
        Zs <- .scaler_apply(Ztr, sz)
      } else {
        sx <- NULL
        sz <- NULL
        Zs <- Ztr
      }
      
      # gates only
      fc <- fcm_fit(Zs, R = Rv, m = mv, max_iter = 200, tol = 1e-5, seed = seed, verbose = FALSE)
      
      fcm_cache[[k]][[f]] <- list(
        centers     = fc$centers,
        sx          = sx,
        sz          = sz,
        standardize = standardize,
        m           = mv,
        R           = Rv,
        tf_y        = tf
      )
    }
  }
  
  # --- full grid over (R, m, lambda) for evaluation ---
  grid  <- expand.grid(R = R_grid, m = m_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE)
  n_cfg <- nrow(grid)
  
  key_RM  <- paste(RM_grid$R, RM_grid$m, sep = ":")
  key_cfg <- paste(grid$R,   grid$m,   sep = ":")
  idx_RM_for_cfg <- match(key_cfg, key_RM)
  
  collect_row <- function(Rv, mv, lam, mse) {
    data.frame(
      R      = as.integer(Rv),
      m      = as.numeric(mv),
      lambda = as.numeric(lam),
      MSE    = as.numeric(mse)
    )
  }
  
  eval_one <- function(ii) {
    Rv  <- grid$R[ii]
    mv  <- grid$m[ii]
    lam <- grid$lambda[ii]
    k   <- idx_RM_for_cfg[ii]
    
    mses <- numeric(length(val_idx_list))
    
    for (f in seq_along(val_idx_list)) {
      va <- val_idx_list[[f]]
      tr <- base::setdiff(seq_len(N), va)
      
      Xtr <- X[tr, , drop = FALSE]
      Ztr <- Z[tr, , drop = FALSE]
      # ytr <- y[tr]
      
      Xva <- X[va, , drop = FALSE]
      Zva <- Z[va, , drop = FALSE]
      # yva <- y[va]
      
      cache_f <- fcm_cache[[k]][[f]]
      
      tf <- cache_f$tf_y
      
      ytr <- cap_and_scale_apply(y_raw[tr], tf)
      yva <- cap_and_scale_apply(y_raw[va], tf)
      
      fit <- fmle_train_from_centers(
        Xtr, ytr, Ztr,
        centers     = cache_f$centers,
        m           = cache_f$m,
        lambda_l1   = lam,
        ridge       = ridge,
        standardize = cache_f$standardize,
        sx          = cache_f$sx,
        sz          = cache_f$sz
      )
      
      pr <- fmle_predict(fit, X_new = Xva, Z_new = Zva, return_se = FALSE)
      mses[f] <- mean((pr$mean - yva)^2)
    }
    
    mse <- mean(mses)
    if (verbose) {
      message(sprintf("[CV] R=%d m=%.2f lambda=%g -> MSE=%.6f", Rv, mv, lam, mse))
    }
    collect_row(Rv, mv, lam, mse)
  }
  
  if (exec == "sequential") {
    res_list <- lapply(seq_len(n_cfg), eval_one)
  } else {
    res_list <- future.apply::future_lapply(seq_len(n_cfg), eval_one, future.seed = TRUE)
  }
  
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  
  best_row <- which.min(res$MSE)
  best <- list(
    R      = res$R[best_row],
    m      = res$m[best_row],
    lambda = res$lambda[best_row],
    mse    = res$MSE[best_row]
  )
  
  list(best = best, table = res)
}

# Optional thin wrapper fmle_cv() can just call fmle_cv_parallel() with exec = "sequential".

# -----------------------------
# Multi-task helpers and CV
# -----------------------------

#' @keywords internal
fmle_train_mt_from_centers <- function(
    X, Y, Z, centers, m,
    lambda_l1   = 0,
    ridge       = 1e-6,
    standardize = TRUE,
    sx = NULL, sz = NULL
){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- as.matrix(Z)
  stopifnot(nrow(X) == nrow(Y), nrow(Z) == nrow(X))
  N  <- nrow(X)
  Tt <- ncol(Y)
  
  if (standardize) {
    stopifnot(!is.null(sx), !is.null(sz))
    Xs <- sweep(sweep(X, 2, sx$mu, "-"), 2, sx$sd, "/")
    Zs <- sweep(sweep(Z, 2, sz$mu, "-"), 2, sz$sd, "/")
  } else {
    Xs <- X
    Zs <- Z
  }
  
  R   <- nrow(centers)
  eps <- 1e-12
  
  CZ    <- tcrossprod(centers, Zs)       # R x N
  nZ    <- rowSums(Zs^2)                 # N
  nC    <- rowSums(centers^2)            # R
  dist2 <- sweep(sweep(-2 * CZ, 2, nZ, "+"), 1, nC, "+")
  dist2 <- pmax(dist2, 0)
  Dpow <- (dist2 + eps)^(-2 / (m - 1))
  U    <- t(Dpow / matrix(colSums(Dpow) + 1e-18,
                          nrow = R, ncol = N, byrow = TRUE))
  w_all <- U^m
  
  # normalize lambda_l1 into RxT matrix (same as fmle_train_mt)
  Tt <- ncol(Y)
  lam_mat <- NULL
  if (length(lambda_l1) == 1L) {
    lam_mat <- matrix(lambda_l1, nrow = R, ncol = Tt)
  } else if (is.vector(lambda_l1) && length(lambda_l1) == R) {
    lam_mat <- matrix(as.numeric(lambda_l1), nrow = R, ncol = Tt)
  } else if (is.vector(lambda_l1) && length(lambda_l1) == Tt) {
    lam_mat <- matrix(rep(as.numeric(lambda_l1), each = R), nrow = R, ncol = Tt)
  } else if (is.matrix(lambda_l1)) {
    stopifnot(nrow(lambda_l1) == R, ncol(lambda_l1) == Tt)
    lam_mat <- lambda_l1
  } else if (is.list(lambda_l1) && length(lambda_l1) == Tt) {
    lam_mat <- matrix(0, nrow = R, ncol = Tt)
    for (t in seq_len(Tt)) {
      stopifnot(length(lambda_l1[[t]]) == R)
      lam_mat[, t] <- as.numeric(lambda_l1[[t]])
    }
  } else {
    stop("lambda_l1 must be scalar, length-R vector, length-T vector, RxT matrix, or list of length T with length-R vectors.")
  }
  
  experts_mt <- vector("list", Tt)
  for (t in seq_len(Tt)) {
    y_t <- Y[, t]
    experts_mt[[t]] <- lapply(seq_len(R), function(r) {
      .fit_expert(Xs, y_t, w_all[, r],
                  lambda_l1 = lam_mat[r, t],
                  ridge     = ridge)
    })
  }
  
  model <- list(
    R          = R,
    m          = m,
    centers    = centers,
    experts_mt = experts_mt,
    U          = U,
    ridge      = ridge,
    lambda_l1  = lam_mat,
    T          = Tt,
    standardize = standardize,
    sx         = sx,
    sz         = sz
  )
  class(model) <- "fmle_mt"
  model
}

#' Multi-task cross-validation for FMLE
#'
#' Performs K-fold multi-task cross-validation over grids of R, m, and
#' lambda_l1 using \code{fmle_train_mt} / \code{fmle_predict_mt}, with
#' optional task weighting and grouped folds.
#'
#' @param X matrix, n x p features.
#' @param Y matrix, n x T responses (e.g., multiple proteins).
#' @param Z matrix, n x d gating features.
#' @param R_grid grid of R values.
#' @param m_grid grid of m values.
#' @param lambda_grid grid of lambda_l1 values.
#' @param folds number of folds.
#' @param groups optional factor of length n for grouped CV.
#' @param ridge ridge penalty.
#' @param standardize logical.
#' @param seed random seed.
#' @param exec "sequential" or "future".
#' @param verbose logical.
#' @param task_weights optional length-T vector of task weights; NULL = equal.
#'
#' @return A list with components:
#' \describe{
#'   \item{best}{List with best (R, m, lambda, mse).}
#'   \item{table}{Data frame of configs and their weighted mean MSE.}
#' }
#' @export
fmle_cv_mt_parallel <- function(
    X, Y, Z,
    R_grid      = c(4, 6, 8),
    m_grid      = c(1.6, 1.8, 2.0),
    lambda_grid = c(0, 1e-3, 1e-2, 3e-2),
    folds       = 5,
    groups      = NULL,
    ridge       = 1e-6,
    standardize = TRUE,
    seed        = 1,
    exec        = c("sequential", "future"),
    verbose     = TRUE,
    task_weights = NULL
){
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  Z <- as.matrix(Z)
  stopifnot(nrow(X) == nrow(Y), nrow(Z) == nrow(X))
  
  exec <- match.arg(exec)
  set.seed(seed)
  N  <- nrow(X)
  Tt <- ncol(Y)
  
  if (is.null(task_weights)) {
    task_weights <- rep(1 / Tt, Tt)
  } else {
    stopifnot(length(task_weights) == Tt)
    task_weights <- task_weights / sum(task_weights)
  }
  
  # --- folds ---
  if (!is.null(groups)) {
    stopifnot(length(groups) == N)
    val_idx_list <- .make_grouped_folds(groups, K = folds, seed = seed)
  } else {
    fold_ids     <- sample(rep(seq_len(folds), length.out = N))
    val_idx_list <- lapply(seq_len(folds), function(f) which(fold_ids == f))
  }
  
  RM_grid <- expand.grid(R = R_grid, m = m_grid, KEEP.OUT.ATTRS = FALSE)
  n_RM    <- nrow(RM_grid)
  fcm_cache <- vector("list", n_RM)
  
  if (verbose) message("Precomputing FCM (gates) per (R, m) and fold for multi-task CV ...")
  
  for (k in seq_len(n_RM)) {
    Rv <- RM_grid$R[k]
    mv <- RM_grid$m[k]
    fcm_cache[[k]] <- vector("list", length(val_idx_list))
    
    for (f in seq_along(val_idx_list)) {
      va <- val_idx_list[[f]]
      tr <- base::setdiff(seq_len(N), va)
      
      Xtr <- X[tr, , drop = FALSE]
      Ztr <- Z[tr, , drop = FALSE]
      Ytr <- Y[tr, , drop = FALSE]
      
      fit_gates <- fmle_train_mt(
        X = Xtr, Y = Ytr, Z = Ztr,
        R = Rv, m = mv,
        lambda_l1   = 0,
        ridge       = ridge,
        standardize = standardize,
        fcm_max_iter = 200, fcm_tol = 1e-5,
        seed        = seed,
        verbose     = FALSE
      )
      
      fcm_cache[[k]][[f]] <- list(
        centers     = fit_gates$centers,
        sx          = fit_gates$sx,
        sz          = fit_gates$sz,
        standardize = fit_gates$standardize,
        m           = fit_gates$m,
        R           = fit_gates$R
      )
    }
  }
  
  grid  <- expand.grid(R = R_grid, m = m_grid, lambda = lambda_grid, KEEP.OUT.ATTRS = FALSE)
  n_cfg <- nrow(grid)
  
  key_RM  <- paste(RM_grid$R, RM_grid$m, sep = ":")
  key_cfg <- paste(grid$R,   grid$m,   sep = ":")
  idx_RM_for_cfg <- match(key_cfg, key_RM)
  
  collect_row <- function(Rv, mv, lam, mse) {
    data.frame(
      R      = as.integer(Rv),
      m      = as.numeric(mv),
      lambda = as.numeric(lam),
      MSE    = as.numeric(mse)
    )
  }
  
  eval_one <- function(ii) {
    Rv  <- grid$R[ii]
    mv  <- grid$m[ii]
    lam <- grid$lambda[ii]
    k   <- idx_RM_for_cfg[ii]
    
    fold_losses <- numeric(length(val_idx_list))
    
    for (f in seq_along(val_idx_list)) {
      va <- val_idx_list[[f]]
      tr <- base::setdiff(seq_len(N), va)
      
      Xtr <- X[tr, , drop = FALSE]
      Ztr <- Z[tr, , drop = FALSE]
      Ytr <- Y[tr, , drop = FALSE]
      
      Xva <- X[va, , drop = FALSE]
      Zva <- Z[va, , drop = FALSE]
      Yva <- Y[va, , drop = FALSE]
      
      cache_f <- fcm_cache[[k]][[f]]
      
      fit <- fmle_train_mt_from_centers(
        X = Xtr, Y = Ytr, Z = Ztr,
        centers     = cache_f$centers,
        m           = cache_f$m,
        lambda_l1   = lam,
        ridge       = ridge,
        standardize = cache_f$standardize,
        sx          = cache_f$sx,
        sz          = cache_f$sz
      )
      
      pr     <- fmle_predict_mt(fit, X_new = Xva, Z_new = Zva, return_se = FALSE)
      Y_hat  <- pr$mean    # N_va x T
      mse_t  <- colMeans((Yva - Y_hat)^2)   # length T
      fold_losses[f] <- sum(task_weights * mse_t)
    }
    
    mse <- mean(fold_losses)
    if (verbose) {
      message(sprintf("[MT-CV] R=%d m=%.2f lambda=%g -> weighted mean MSE=%.6f",
                      Rv, mv, lam, mse))
    }
    collect_row(Rv, mv, lam, mse)
  }
  
  if (exec == "sequential") {
    res_list <- lapply(seq_len(n_cfg), eval_one)
  } else {
    res_list <- future.apply::future_lapply(seq_len(n_cfg), eval_one, future.seed = TRUE)
  }
  
  res <- do.call(rbind, res_list)
  rownames(res) <- NULL
  
  best_row <- which.min(res$MSE)
  best <- list(
    R      = res$R[best_row],
    m      = res$m[best_row],
    lambda = res$lambda[best_row],
    mse    = res$MSE[best_row]
  )
  
  list(best = best, table = res)
}
