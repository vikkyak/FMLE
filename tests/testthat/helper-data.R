make_fmle_demo_data <- function(seed = 1, n_train = 120, n_test = 40,
                                p = 30, d = 8, T = 3, R = 3) {
  set.seed(seed)
  
  n_total <- n_train + n_test
  z_latent <- matrix(rnorm(n_total * d), nrow = n_total, ncol = d)
  centers <- matrix(rnorm(R * d, sd = 1.5), nrow = R, ncol = d)
  
  dist2 <- sapply(seq_len(R), function(r) {
    rowSums((z_latent - matrix(centers[r, ], n_total, d, byrow = TRUE))^2)
  })
  logits <- -dist2
  logits <- logits - apply(logits, 1, max)
  gates <- exp(logits)
  gates <- gates / rowSums(gates)
  
  loadings <- matrix(rnorm(d * p, sd = 0.6), nrow = d, ncol = p)
  X <- z_latent %*% loadings + matrix(rnorm(n_total * p, sd = 0.4),
                                      nrow = n_total, ncol = p)
  
  beta <- array(rnorm(R * p * T, sd = 0.15), dim = c(R, p, T))
  intercepts <- matrix(rnorm(R * T, mean = 1.5, sd = 0.3), nrow = R, ncol = T)
  
  Y <- matrix(0, nrow = n_total, ncol = T)
  for (t in seq_len(T)) {
    mu_r <- sapply(seq_len(R), function(r) intercepts[r, t] + X %*% beta[r, , t])
    eta <- rowSums(gates * mu_r) + rnorm(n_total, sd = 0.25)
    
    ## strictly positive response, compatible with cap/log1p in fmle_cv_parallel()
    Y[, t] <- exp(eta / 3)
  }
  
  colnames(X) <- paste0("Gene", seq_len(p))
  colnames(z_latent) <- paste0("PC", seq_len(d))
  colnames(Y) <- paste0("Protein", seq_len(T))
  rownames(X) <- paste0("Cell", seq_len(n_total))
  rownames(z_latent) <- rownames(X)
  rownames(Y) <- rownames(X)
  
  idx_tr <- seq_len(n_train)
  idx_te <- seq.int(n_train + 1, n_total)
  
  list(
    X_train = X[idx_tr, , drop = FALSE],
    Y_train = Y[idx_tr, , drop = FALSE],
    Z_train = z_latent[idx_tr, , drop = FALSE],
    X_test  = X[idx_te, , drop = FALSE],
    Y_test  = Y[idx_te, , drop = FALSE],
    Z_test  = z_latent[idx_te, , drop = FALSE]
  )
}