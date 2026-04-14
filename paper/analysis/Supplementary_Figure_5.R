suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
  library(purrr)
  library(tidyr)
})

#--------------------------------------------#
# Cross data sets BMMC
#--------------------------------------------#

base <- file.path(cfg$out_root, "transfer_preds_BMMC")
dir.create(base, recursive=TRUE, showWarnings=FALSE)

bench_name <- c(
  `A` = "Young healthy BM",
  `B` = "BMMC CITE-seq"
  
) 

read_one_condition <- function(train_ds, test_ds) {
  cond_dir <- file.path(base, sprintf("ctp_%s_to_%s", train_ds, test_ds))
  tag <- paste0(train_ds, test_ds)  # "AB", "AC", ...
  
  fmle_path <- file.path(cond_dir, sprintf("res_fmle_%s.csv", tag))
  scl_path  <- file.path(cond_dir, sprintf("res_scl_%s.csv",  tag))
  ctp_path  <- file.path(cond_dir, sprintf("res_ctp_%s.csv",  tag))
  
  stopifnot(file.exists(fmle_path), file.exists(scl_path), file.exists(ctp_path))
  
  bind_rows(
    read_csv(fmle_path, show_col_types = FALSE) %>%
      transmute(train_ds, test_ds, method = "FMLE",    protein, Pearson),
    read_csv(scl_path,  show_col_types = FALSE) %>%
      transmute(train_ds, test_ds, method = "scLinear", protein, Pearson),
    read_csv(ctp_path,  show_col_types = FALSE) %>%
      transmute(train_ds, test_ds, method = "cTPnet",  protein, Pearson)
  ) %>%
    mutate(
      condition = paste0(bench_name[train_ds], " \u2192 ", bench_name[test_ds]),
      train_label = bench_name[train_ds],
      test_label  = bench_name[test_ds]
    )
}


conds <- tribble(
  ~train_ds, ~test_ds,
  "A","B",
  "B","A"
)


df_all <- pmap_dfr(conds, read_one_condition) %>%
  filter(is.finite(Pearson)) %>%
  mutate(
    method = factor(method, levels = c("scLinear","cTPnet","FMLE")),
    condition = factor(condition, levels = unique(condition))
  )


p_transfer <- ggplot(df_all, aes(x = method, y = Pearson, color = method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  facet_wrap(~ condition, nrow = 2) +
  theme_classic(base_size = 10) +
  labs(x = NULL, y = "Per-protein Pearson correlation") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )

df_pair2 <- df_all %>%
  filter(method %in% c("FMLE","scLinear","cTPnet")) %>%
  select(train_ds, test_ds, condition, protein, method, Pearson) %>%
  pivot_wider(names_from = method, values_from = Pearson) %>%
  mutate(
    gain_scLinear = FMLE - scLinear,
    gain_cTPnet   = FMLE - cTPnet
  ) %>%
  pivot_longer(
    cols = c(scLinear, cTPnet),
    names_to = "baseline_method",
    values_to = "baseline"
  ) %>%
  mutate(
    gain = FMLE - baseline,
    baseline_method = factor(baseline_method, levels = c("scLinear","cTPnet"))
  ) %>%
  filter(is.finite(FMLE), is.finite(baseline))

df_pair2 <- df_pair2 %>%
  mutate(
    better = ifelse(gain > 0, "FMLE", as.character(baseline_method)),
    better = factor(better, levels = c("FMLE", "scLinear", "cTPnet"))
  )


make_scatter_one <- function(dat, baseline_name, cond_lab, xlab, ylab,
                             show_legend = TRUE, title_text = NULL) {
  sub <- dat %>%
    filter(baseline_method == baseline_name, condition == cond_lab)
  sub$better <- factor(sub$better, levels = c("FMLE", "scLinear", "cTPnet"))
  
  lims <- range(c(sub$baseline, sub$FMLE), na.rm = TRUE)
  
  ggplot(sub, aes(x = baseline, y = FMLE)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = 2, linewidth = 0.5, color = "black"
    ) +
    geom_point(
      aes(color = better),
      shape = 16, size = 2.6, alpha = 0.9
    ) +
    scale_color_manual(
      name = "Better method",
      values = c(
        "FMLE"    = "#38b48b",
        "scLinear"= "#e67e22",
        "cTPnet"  = "#e67e22"
      ),
      breaks = c("FMLE", "scLinear", "cTPnet")
    ) +
    coord_cartesian(xlim = lims, ylim = lims) +
    labs(
      x = xlab,
      y = ylab,
      title = ifelse(is.null(title_text), cond_lab, title_text)
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(size = 10, face = "bold", hjust = 0),
      legend.position = if (show_legend) "right" else "none",
      legend.title    = element_text(size = 10),
      legend.text     = element_text(size = 10)
    )
}
cond_levels <- levels(df_all$condition)


p_ab_scl <- make_scatter_one(df_pair2, "scLinear", cond_levels[1],
                             "scLinear", "FMLE",
                             show_legend = TRUE, title_text = cond_levels[1])

p_ba_scl <- make_scatter_one(df_pair2, "scLinear", cond_levels[2],
                             "scLinear", "FMLE",
                             show_legend = TRUE, title_text = cond_levels[2])

p_ab_ctp <- make_scatter_one(df_pair2, "cTPnet", cond_levels[1],
                             "cTPnet", "FMLE",
                             show_legend = TRUE, title_text = cond_levels[1])

p_ba_ctp <- make_scatter_one(df_pair2, "cTPnet", cond_levels[2],
                             "cTPnet", "FMLE",
                             show_legend = TRUE, title_text = cond_levels[2])


p_transfer_ab <- ggplot(
  df_all %>% filter(condition == unique(as.character(condition))[1]),
  aes(x = method, y = Pearson, color = method)
) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Per-protein Pearson",
    title = unique(as.character(df_all$condition))[1]
  ) +
  theme(
    plot.title      = element_text(size = 10, face = "bold", hjust = 0),
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 10)
  )

p_transfer_ba <- ggplot(
  df_all %>% filter(condition == unique(as.character(condition))[2]),
  aes(x = method, y = Pearson, color = method)
) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Per-protein Pearson",
    title = unique(as.character(df_all$condition))[2]
  ) +
  theme(
    plot.title      = element_text(size = 10, face = "bold", hjust = 0),
    legend.title    = element_text(size = 10),
    legend.text     = element_text(size = 10)
  )


p_ab_scl <- p_ab_scl + theme(legend.position = "right")
p_ab_ctp <- p_ab_ctp + theme(legend.position = "none")
p_ba_scl <- p_ba_scl + theme(legend.position = "none")
p_ba_ctp <- p_ba_ctp + theme(legend.position = "none")

extended_fig <- ((p_transfer_ab | p_ab_scl | p_ab_ctp) /
                   (p_transfer_ba | p_ba_scl | p_ba_ctp)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag = element_text(size = 10, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )


out_dir <- file.path(base)

ggsave(file.path(out_dir, "Supplementary_Figure_5.pdf"),
       plot = extended_fig ,
       device = cairo_pdf,
       width = 11, height = 6, units = "in",
       dpi = 300,
       limitsize = FALSE)





cmp_all <- df_all %>%
  select(condition, protein, method, Pearson) %>%
  tidyr::pivot_wider(names_from = method, values_from = Pearson) %>%
  filter(is.finite(FMLE), is.finite(scLinear), is.finite(cTPnet)) %>%
  mutate(
    delta_scLinear = FMLE - scLinear,
    delta_cTPnet   = FMLE - cTPnet
  )

# 2) mean gain per condition (what you did with mean(cmp$...))
mean_gain_by_cond <- cmp_all %>%
  group_by(condition) %>%
  summarise(
    n_proteins = n(),
    mean_gain_scLinear = mean(delta_scLinear, na.rm = TRUE),
    mean_gain_cTPnet   = mean(delta_cTPnet,   na.rm = TRUE),
    .groups = "drop"
  )


mean_gain_by_cond


summary_gain_by_cond <- cmp_all %>%
  group_by(condition) %>%
  summarise(
    n_proteins = n(),
    
    scL_min = min(delta_scLinear, na.rm=TRUE),
    scL_q1  = quantile(delta_scLinear, 0.25, na.rm=TRUE),
    scL_med = median(delta_scLinear, na.rm=TRUE),
    scL_mean= mean(delta_scLinear, na.rm=TRUE),
    scL_q3  = quantile(delta_scLinear, 0.75, na.rm=TRUE),
    scL_max = max(delta_scLinear, na.rm=TRUE),
    
    cTP_min = min(delta_cTPnet, na.rm=TRUE),
    cTP_q1  = quantile(delta_cTPnet, 0.25, na.rm=TRUE),
    cTP_med = median(delta_cTPnet, na.rm=TRUE),
    cTP_mean= mean(delta_cTPnet, na.rm=TRUE),
    cTP_q3  = quantile(delta_cTPnet, 0.75, na.rm=TRUE),
    cTP_max = max(delta_cTPnet, na.rm=TRUE),
    
    .groups = "drop"
  )

summary_gain_by_cond


wilcox_by_cond <- cmp_all %>%
  group_by(condition) %>%
  summarise(
    p_FMLE_vs_scLinear = wilcox.test(FMLE, scLinear, paired=TRUE)$p.value,
    p_FMLE_vs_cTPnet   = wilcox.test(FMLE, cTPnet,   paired=TRUE)$p.value,
    .groups = "drop"
  )

wilcox_by_cond





