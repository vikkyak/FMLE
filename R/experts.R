#' @keywords internal
.fit_expert <- function(Xs, y, w, lambda_l1=0, ridge=1e-6){
  Xs <- as.matrix(Xs); y <- as.numeric(y); W <- as.numeric(w)
  N <- nrow(Xs); p <- ncol(Xs)
  if (sum(W) <= 0){
    return(list(beta0=0, beta=rep(0,p), sigma2=NA_real_,
                cov=diag(NA_real_, p+1), df=1, method="none"))
  }
  Xtil <- .mat_colbind1(Xs)
  
  if (lambda_l1 <= 0){
    WX   <- Xtil * W
    XtWX <- t(Xtil) %*% WX
    XtWy <- crossprod(Xtil, W * y)
    inv  <- .solve_spd(XtWX, ridge, pen_intercept=FALSE)
    theta <- inv %*% XtWy
    mu <- as.vector(Xtil %*% theta)
    Neff <- .effective_n(W)
    df <- min(sum(W>0), p+1)
    sigma2 <- if (Neff > df) sum(W*(y-mu)^2) / (Neff - df) else NA_real_
    cov <- if (is.na(sigma2)) matrix(NA_real_, p+1, p+1) else sigma2 * inv
    return(list(beta0=as.numeric(theta[1]), beta=as.numeric(theta[-1]),
                sigma2=sigma2, cov=cov, df=df, method="wls"))
  }
  
  if (!requireNamespace("glmnet", quietly=TRUE)){
    stop("glmnet is required when lambda_l1 > 0. Install with install.packages('glmnet').")
  }
  
  if (isTRUE(getOption("fmle.cv_lasso", FALSE))){
    cv <- glmnet::cv.glmnet(Xs, y, alpha=1, weights=W, intercept=TRUE, standardize=FALSE)
    lambda_use <- cv$lambda.min
  } else lambda_use <- lambda_l1
  
  fit <- glmnet::glmnet(Xs, y, alpha=1, lambda=lambda_use,
                        intercept=TRUE, standardize=FALSE,
                        weights=W, thresh=1e-7)
  beta_lasso <- as.vector(fit$beta)
  active <- which(beta_lasso != 0L)
  if (length(active) == 0L) active <- seq_len(p)  # fallback: ridge on all
  
  Xs_sel  <- Xs[, active, drop=FALSE]
  Xtil_sel <- cbind(1, Xs_sel)
  WXs   <- Xtil_sel * W
  XtWXs <- t(Xtil_sel) %*% WXs
  XtWys <- crossprod(Xtil_sel, W * y)
  invs  <- .solve_spd(XtWXs, ridge, pen_intercept=FALSE)
  theta_sel <- invs %*% XtWys
  mu <- as.vector(Xtil_sel %*% theta_sel)
  
  Neff <- .effective_n(W)
  df <- ncol(Xtil_sel)
  sigma2 <- if (Neff > df) sum(W*(y-mu)^2) / (Neff - df) else NA_real_
  
  cov_full <- matrix(0, p+1, p+1)
  idx <- c(1, active + 1)
  cov_full[idx, idx] <- if (is.na(sigma2)) NA_real_ else sigma2 * invs
  
  beta_full <- numeric(p); beta_full[active] <- as.numeric(theta_sel[-1])
  beta0 <- as.numeric(theta_sel[1])
  
  list(beta0=beta0, beta=beta_full, sigma2=sigma2, cov=cov_full, df=df, method="lasso+postridge")
}
