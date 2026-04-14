entropy_shift_stats <- function(H_tr, H_te, B = 1000, seed = 1) {
  set.seed(seed)
  H_tr <- H_tr[is.finite(H_tr)]
  H_te <- H_te[is.finite(H_te)]
  
  med_tr <- median(H_tr); med_te <- median(H_te)
  dmed   <- med_te - med_tr
  iqr_tr <- quantile(H_tr, c(.25,.75))
  iqr_te <- quantile(H_te, c(.25,.75))
  
  # bootstrap CI for delta-median
  ntr <- length(H_tr); nte <- length(H_te)
  dboot <- replicate(B, {
    median(sample(H_te, nte, replace=TRUE)) - median(sample(H_tr, ntr, replace=TRUE))
  })
  ci <- quantile(dboot, c(.025,.975))
  
  tibble::tibble(
    n_tr = ntr, n_te = nte,
    med_tr = med_tr, q25_tr = iqr_tr[1], q75_tr = iqr_tr[2],
    med_te = med_te, q25_te = iqr_te[1], q75_te = iqr_te[2],
    delta_median = dmed,
    ci_low = ci[1], ci_high = ci[2],
    p_wilcox = suppressWarnings(wilcox.test(H_tr, H_te)$p.value)
  )
}