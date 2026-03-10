#' Train FMLE model for a single task
#'
#' @param X matrix, n x p features (e.g., HVG expression)
#' @param y numeric vector length n (response)
#' @param Z matrix, n x d gating features (e.g., PCs). If NULL, uses X.
#' @param R number of experts
#' @param m fuzzifier (>1)
#' @param lambda_l1 L1 penalty (scalar or length-R vector)
#' @param ridge small ridge penalty for stability
#' @param standardize logical; standardize X and Z
#' @param fcm_max_iter max FCM iterations
#' @param fcm_tol tolerance for FCM convergence
#' @param seed random seed
#' @param verbose logical
#'
#' @return object of class \code{"fmle"}
#' @export
fmle_train <- function(X, y, Z = NULL,
                       R = 6, m = 1.8,
                       lambda_l1 = 0, ridge = 1e-6,
                       standardize = TRUE,
                       fcm_max_iter = 200, fcm_tol = 1e-5,
                       seed = 1, verbose = FALSE) {
  X <- as.matrix(X); y <- as.numeric(y); if (is.null(Z)) Z <- X; Z <- as.matrix(Z)
  stopifnot(nrow(X) == nrow(Z), nrow(X) == length(y))
  
  if (standardize) { sx <- .scaler_fit(X); Xs <- .scaler_apply(X, sx)
  sz <- .scaler_fit(Z); Zs <- .scaler_apply(Z, sz) } else {
    sx <- NULL; sz <- NULL; Xs <- X; Zs <- Z }
  
  fc <- fcm_fit(Zs, R = R, m = m, max_iter = fcm_max_iter, tol = fcm_tol, seed = seed, verbose = verbose)
  U <- fc$U; C <- fc$centers; w_all <- U^m
  # print(summary(w_all))
  ## NEW: normalize lambda vector
  lam_vec <- if (length(lambda_l1) == 1L) rep(lambda_l1, R) else {
    stopifnot(length(lambda_l1) == R); as.numeric(lambda_l1)
  }
  
  experts <- lapply(seq_len(R), function(r) .fit_expert(Xs, y, w_all[, r],
                                                        lambda_l1 = lam_vec[r],
                                                        ridge = ridge))
  
  structure(list(R=R, m=m, centers=C, experts=experts,
                 ridge=ridge, lambda_l1=lam_vec,   # store the vector
                 standardize=standardize, sx=sx, sz=sz, U=U),
            class="fmle")
}


#' Predict from an FMLE model on single task
#'
#' @param model object of class \code{"fmle"}
#' @param X_new new design matrix (cells x p)
#' @param Z_new new gating features (cells x d); if NULL, uses X_new
#' @param return_se logical; if TRUE, also return uncertainty estimates
#'
#' @return list with components \code{mean}, \code{alpha}, \code{mu_r},
#'   and optionally \code{var}, \code{se}
#' @export
fmle_predict <- function(model, X_new, Z_new=NULL, return_se=TRUE){
  stopifnot(inherits(model,"fmle"))
  X_new <- as.matrix(X_new); if (is.null(Z_new)) Z_new <- X_new; Z_new <- as.matrix(Z_new)
  
  if (isTRUE(model$standardize)){ Xs <- .scaler_apply(X_new, model$sx); Zs <- .scaler_apply(Z_new, model$sz) }
  else { Xs <- X_new; Zs <- Z_new }
  
  R <- model$R; m <- model$m; C <- model$centers
  eps <- 1e-12; N <- nrow(Xs)
  
  CZ <- tcrossprod(C, Zs); nZ <- rowSums(Zs^2); nC <- rowSums(C^2)
  dist2 <- sweep(sweep(-2*CZ, 2, nZ, "+"), 1, nC, "+")
  dist2 <- pmax(dist2, 0)
  Dpow  <- (dist2 + eps)^(-2/(m-1))
  Unew  <- t(Dpow / matrix(colSums(Dpow) + 1e-18, nrow=R, ncol=N, byrow=TRUE))
  # print(summary(Unew))
  A     <- Unew^m
  alpha <- A / (rowSums(A) + 1e-18)
  
  mu_r <- matrix(0, N, R)
  for (r in 1:R){
    ex <- model$experts[[r]]
    mu_r[,r] <- ex$beta0 + as.vector(Xs %*% ex$beta)
  }
  mu <- rowSums(alpha * mu_r)
  if (!return_se) return(list(mean=mu, alpha=alpha, mu_r=mu_r))
  
  # sigma2_r <- sapply(model$experts, `[[`, "sigma2"); sigma2_r[is.na(sigma2_r)] <- 0
  # var_ale  <- as.vector((alpha^2) %*% sigma2_r)
  
  sigma2_r <- vapply(model$experts, `[[`, numeric(1), "sigma2")
  sigma2_r[!is.finite(sigma2_r)] <- 0
  
  ## correct aleatoric variance: sum_r alpha_{ir} * sigma_r^2
  var_ale <- as.vector(alpha %*% sigma2_r)
  
  Xtil <- .mat_colbind1(Xs)
  var_par <- numeric(N)
  for (r in 1:R){
    covr <- model$experts[[r]]$cov
    if (!is.matrix(covr) || anyNA(covr)) next
    CX <- Xtil %*% covr
    q  <- rowSums(CX * Xtil)  # diag(Xtil Cov Xtil^T)
    var_par <- var_par + (alpha[,r]^2) * q
  }
  v <- pmax(0, var_ale + var_par)
  list(mean=mu, var=v, se=sqrt(v), alpha=alpha, mu_r=mu_r)
}

#' Train FMLE model for a multi-task
#'
#' @param X matrix, n x p features
#' @param Y matrix, n x T responses (T tasks, e.g., multiple proteins)
#' @param Z matrix, n x d gating features
#' @param R number of experts
#' @param m fuzzifier
#' @param lambda_l1 penalty (scalar, length-R, length-T, RxT matrix, or list)
#' @param ridge ridge penalty
#' @param standardize logical
#' @param fcm_max_iter max FCM iterations
#' @param fcm_tol tolerance
#' @param seed random seed
#' @param verbose logical
#'
#' @return object of class \code{"fmle_mt"}
#' @export
fmle_train_mt <- function(X, Y, Z = NULL,
                          R = 6, m = 1.8,
                          lambda_l1 = 0, ridge = 1e-6,
                          standardize = TRUE,
                          fcm_max_iter = 200, fcm_tol = 1e-5,
                          seed = 1, verbose = FALSE) {
  X <- as.matrix(X); Y <- as.matrix(Y); if (is.null(Z)) Z <- X; Z <- as.matrix(Z)
  stopifnot(nrow(X) == nrow(Y), nrow(Z) == nrow(X))
  Tt <- ncol(Y)
  
  # --- standardize (optional) ---
  if (standardize) {
    sx <- .scaler_fit(X); Xs <- .scaler_apply(X, sx)
    sz <- .scaler_fit(Z); Zs <- .scaler_apply(Z, sz)
  } else {
    sx <- NULL; sz <- NULL; Xs <- X; Zs <- Z
  }
  
  # --- shared FCM (gates) ---
  fc <- fcm_fit(Zs, R = R, m = m, max_iter = fcm_max_iter, tol = fcm_tol, seed = seed, verbose = verbose)
  U <- fc$U; w_all <- U^m
  
  # --- normalize lambda_l1 to an R x T matrix (lam_mat) ---
  lam_mat <- NULL
  if (length(lambda_l1) == 1L) {
    lam_mat <- matrix(lambda_l1, nrow = R, ncol = Tt)
  } else if (is.vector(lambda_l1) && length(lambda_l1) == R) {
    lam_mat <- matrix(as.numeric(lambda_l1), nrow = R, ncol = Tt, byrow = FALSE)
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
  
  # --- fit per-task, per-expert ---
  experts_mt <- vector("list", Tt)
  for (t in 1:Tt) {
    y <- Y[, t]
    experts_mt[[t]] <- lapply(seq_len(R), function(r) {
      .fit_expert(Xs, y, w_all[, r],
                  lambda_l1 = lam_mat[r, t],
                  ridge = ridge)
    })
  }
  
  model <- list(
    R = R, m = m, centers = fc$centers,
    experts_mt = experts_mt, U = U,
    ridge = ridge, lambda_l1 = lam_mat,  # store the RxT matrix
    T = Tt,
    standardize = standardize, sx = sx, sz = sz,
    fcm_objective = fc$J, fcm_iters = fc$iters
  )
  class(model) <- "fmle_mt"
  model
}


#' Predict from a multi-task FMLE model
#'
#' @param model object of class \code{"fmle_mt"}
#' @param X_new new feature matrix
#' @param Z_new new gating feature matrix
#' @param return_se logical
#'
#' @return list with \code{mean} (n x T), \code{alpha}, and optional \code{var}, \code{se}
#' @export
fmle_predict_mt <- function(model, X_new, Z_new=NULL, return_se=TRUE){
  stopifnot(inherits(model,"fmle_mt"))
  X_new <- as.matrix(X_new); if (is.null(Z_new)) Z_new <- X_new; Z_new <- as.matrix(Z_new)
  
  if (isTRUE(model$standardize)){ Xs <- .scaler_apply(X_new, model$sx); Zs <- .scaler_apply(Z_new, model$sz) }
  else { Xs <- X_new; Zs <- Z_new }
  
  R <- model$R; m <- model$m; C <- model$centers
  eps <- 1e-12; N <- nrow(Xs); Tt <- model$T
  
  CZ <- tcrossprod(C, Zs); nZ <- rowSums(Zs^2); nC <- rowSums(C^2)
  dist2 <- sweep(sweep(-2*CZ, 2, nZ, "+"), 1, nC, "+")
  dist2 <- pmax(dist2, 0)
  Dpow  <- (dist2 + eps)^(-2/(m-1))
  Unew  <- t(Dpow / matrix(colSums(Dpow) + 1e-18, nrow=R, ncol=N, byrow=TRUE))
  A     <- Unew^m
  alpha <- A / (rowSums(A) + 1e-18)
  
  means <- matrix(0, N, Tt)
  vars  <- if (return_se) matrix(0, N, Tt) else NULL
  Xtil  <- .mat_colbind1(Xs)
  
  for (t in 1:Tt){
    mu_r <- matrix(0, N, R)
    sigma2_r <- numeric(R)
    for (r in 1:R){
      ex <- model$experts_mt[[t]][[r]]
      mu_r[,r] <- ex$beta0 + as.vector(Xs %*% ex$beta)
      sigma2_r[r] <- ifelse(is.na(ex$sigma2), 0, ex$sigma2)
    }
    means[,t] <- rowSums(alpha * mu_r)
    
    if (return_se){
      # var_ale <- as.vector((alpha^2) %*% sigma2_r)
      var_ale <- as.vector(alpha %*% sigma2_r)
      var_par <- numeric(N)
      for (r in 1:R){
        covr <- model$experts_mt[[t]][[r]]$cov
        if (!is.matrix(covr) || anyNA(covr)) next
        CX <- Xtil %*% covr
        q  <- rowSums(CX * Xtil)
        var_par <- var_par + (alpha[,r]^2) * q
      }
      vars[,t] <- pmax(0, var_ale + var_par)
    }
  }
  if (return_se) list(mean=means, var=vars, se=sqrt(vars), alpha=alpha) else list(mean=means, alpha=alpha)
}


#' @export
print.fmle <- function(x, ...) {
  lam <- x$lambda_l1
  lam_str <- if (length(lam) == 1L) {
    sprintf("%.3g", lam)
  } else {
    sprintf("vector[%d]{min=%.3g, med=%.3g, max=%.3g}",
            length(lam), min(lam), stats::median(lam), max(lam))
  }
  cat(sprintf("FMLE: R=%d, m=%.3f, lambda_l1=%s, ridge=%.1e\n",
              x$R, x$m, lam_str, x$ridge))
  s2 <- vapply(x$experts, function(e) e$sigma2, numeric(1))
  cat(sprintf("Experts residual variance: min=%.4g, median=%.4g, max=%.4g\n",
              min(s2, na.rm=TRUE), stats::median(s2, na.rm=TRUE), max(s2, na.rm=TRUE)))
  invisible(x)
}








