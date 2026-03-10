#' @keywords internal
.scaler_fit   <- function(X){ mu <- colMeans(X); sd <- apply(X,2,sd); sd[sd==0] <- 1; list(mu=mu, sd=sd) }

#' @keywords internal
.scaler_apply <- function(X,sc) sweep(sweep(X,2,sc$mu,"-"),2,sc$sd,"/")

#' @keywords internal
.effective_n  <- function(w){ s1 <- sum(w); s2 <- sum(w^2); if (s2==0) 0 else (s1^2)/s2 }

#' @keywords internal
.mat_colbind1 <- function(X) cbind(Intercept=1, X)

#' @keywords internal
.solve_spd <- function(G, ridge=1e-6, pen_intercept=FALSE){
  P <- nrow(G); Pen <- diag(P); if (!pen_intercept) Pen[1,1] <- 0
  A <- G + ridge * Pen
  R <- tryCatch(chol(A), error=function(e) NULL)
  if (!is.null(R)) chol2inv(R) else qr.solve(A)
}

#' @keywords internal
.n_unique_rows <- function(M, digits = 6) {
  M <- as.matrix(M)
  if (is.numeric(M)) M <- round(M, digits)  # avoid floating-point near-duplicates
  nrow(unique(M))
}

#' @keywords internal
.make_grouped_folds <- function(groups, K = 5, seed = 1) {
  set.seed(seed)
  g    <- factor(groups)
  glev <- levels(g)
  glev <- sample(glev, length(glev))                 # shuffle groups
  bins <- split(glev, rep(seq_len(K), length.out = length(glev)))
  lapply(bins, function(b) which(g %in% b))
}

# public helper if you want
#' Cap and log-scale a response vector
#'
#' @param y numeric response (e.g., ADT counts)
#' @param q upper quantile for capping (default 0.995)
#' @param eps Small positive constant used for numerical stability.
#'
#' @return scaled numeric vector
#' @export
cap_and_scale_fit <- function(y, q = 0.995, eps = 1e-8) {
  y <- as.numeric(y)
  cap <- as.numeric(stats::quantile(y, q, na.rm = TRUE, names = FALSE))
  y_cap <- pmin(y, cap)
  y_log <- log1p(y_cap)
  mu  <- mean(y_log, na.rm = TRUE)
  sdv <- stats::sd(y_log, na.rm = TRUE)
  if (!is.finite(sdv) || sdv < eps) sdv <- 1.0
  list(cap = cap, mu = mu, sd = sdv, q = q)
}

cap_and_scale_apply <- function(y, tf) {
  y <- as.numeric(y)
  y_cap <- pmin(y, tf$cap)
  y_log <- log1p(y_cap)
  (y_log - tf$mu) / tf$sd
}
