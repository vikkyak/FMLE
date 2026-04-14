entropy_error_stats <- function(H, y, yhat, hi = 0.9, lo = 0.1) {
  H <- H[is.finite(H)]
  resid <- abs(y - yhat)
  ok <- is.finite(resid) & is.finite(H)
  H <- H[ok]; resid <- resid[ok]
  
  q_hi <- stats::quantile(resid, hi)
  q_lo <- stats::quantile(resid, lo)
  
  grp <- ifelse(resid >= q_hi, "high_error",
                ifelse(resid <= q_lo, "low_error", NA))
  keep <- !is.na(grp)
  df <- data.frame(H=H[keep], grp=grp[keep], resid=resid[keep])
  
  med_hi <- median(df$H[df$grp=="high_error"])
  med_lo <- median(df$H[df$grp=="low_error"])
  dmed   <- med_hi - med_lo
  
  tibble::tibble(
    n_high = sum(df$grp=="high_error"),
    n_low  = sum(df$grp=="low_error"),
    med_high = med_hi,
    med_low  = med_lo,
    delta_median = dmed,
    spearman_H_resid = suppressWarnings(cor(df$H, df$resid, method="spearman")),
    p_wilcox = suppressWarnings(wilcox.test(H ~ grp, data=df)$p.value)
  )
}