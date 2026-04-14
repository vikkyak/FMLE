#' Fuzzy c-means clustering for gating
#'
#' @param Z numeric matrix (cells × latent dimensions, e.g. PCs)
#' @param R number of clusters / experts
#' @param m fuzzifier (>1)
#' @param max_iter maximum iterations
#' @param tol convergence tolerance
#' @param seed random seed
#' @param verbose logical; print progress
#'
#' @return list with membership matrix U, centers, objective J, etc.
#' @export
fcm_fit <- function(Z, R, m = 1.8, max_iter = 200, tol = 1e-5, seed = 1, verbose = FALSE) {
  stopifnot(m > 1, R >= 2)
  set.seed(seed)
  Z <- as.matrix(Z)
  N <- nrow(Z); eps <- 1e-12; pow <- -1 / (m - 1)
  
  km <- stats::kmeans(Z, centers = R, nstart = 5, iter.max = 50)
  C <- km$centers
  U_prev <- NULL; J_prev <- Inf
  
  ## jitter any duplicated centers (rare but deadly for gating)
  dup <- which(duplicated(round(C, 6)))
  if (length(dup)) {
    set.seed(seed + 123)
    C[dup, ] <- C[dup, ] + matrix(stats::rnorm(length(dup) * ncol(C), sd = 1e-3),
                                  nrow = length(dup))
  }
  
  
  dist2_CN <- function(C, Z) {
    CZ <- tcrossprod(C, Z)          # R x N
    nZ <- rowSums(Z^2)              # N
    nC <- rowSums(C^2)              # R
    sweep(sweep(-2 * CZ, 2, nZ, "+"), 1, nC, "+")  # R x N
  }
  
  dist2 <- dist2_CN(C, Z)
  dist2 <- pmax(dist2, 0)
  for (iter in 1:max_iter) {
    ## soft membership first
    Dpow <- (dist2 + eps) ^ pow                    # R x N
    den  <- colSums(Dpow) + 1e-18                  # N (guard)
    U    <- t(Dpow / matrix(den, nrow = nrow(Dpow), ncol = ncol(Dpow), byrow = TRUE))  # N x R
    
    ## overwrite ONLY rows with exact zero distance (one-hot, not global)
    hits <- dist2 < 1e-14                          # R x N logical
    if (any(hits)) {
      zero_cols <- which(colSums(hits) > 0)        # cells with a zero distance
      for (j in zero_cols) {
        rstar <- which.min(dist2[, j])
        U[j, ] <- 0
        U[j, rstar] <- 1
      }
    }
    
    ## centroid update
    Um   <- U ^ m
    num  <- t(Um) %*% Z                            # R x d
    denC <- colSums(Um)                            # R
    C    <- num / denC
    
    ## recompute distances + objective
    dist2 <- dist2_CN(C, Z)
    dist2 <- pmax(dist2, 0)
    J     <- sum((U ^ m) * t(dist2))
    if (verbose && iter %% 10 == 0)
      message(sprintf("[FCM] iter=%d J=%.6g", iter, J))
    
    ## stopping
    if (!is.null(U_prev)) {
      relU <- sqrt(mean((U - U_prev)^2)) / (sqrt(mean(U_prev^2)) + eps)
      relJ <- abs(J_prev - J) / (abs(J_prev) + eps)
      if (relU < tol || relJ < tol) break
    }
    U_prev <- U; J_prev <- J
  }
  
  list(U = U, centers = C, m = m, R = R, J = J, iters = iter)
}
