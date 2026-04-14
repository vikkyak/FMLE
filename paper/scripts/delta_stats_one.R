# For each protein: fit1 vs fit2
# Save: delta_R2, delta_RSS_frac, p_anova, p_adj

delta_stats_one <- function(df) {
  df <- df[is.finite(df$rna) & is.finite(df$y), ]
  df$hard <- droplevels(factor(df$hard))
  df$celltype <- droplevels(factor(df$celltype))
  if (nlevels(df$hard) < 2 || nlevels(df$celltype) < 2) return(NULL)
  
  fit1 <- lm(y ~ rna*celltype, data=df)
  fit2 <- lm(y ~ rna*celltype + hard + rna:hard, data=df)
  
  a <- anova(fit1, fit2)
  p <- a$`Pr(>F)`[2]
  
  RSS1 <- sum(resid(fit1)^2)
  RSS2 <- sum(resid(fit2)^2)
  
  tibble(
    delta_R2 = 1 - RSS2/RSS1,
    delta_RSS_frac = (RSS1 - RSS2)/RSS1,
    p_anova = p
  )
}

# after looping across proteins:
# res$p_adj <- p.adjust(res$p_anova, method="BH")