suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})


#================================================================
# 1)  Within hashtag oligonucleotide (HTO) transfer PBMC Kaggle
#================================================================

bench <- 1  
ds      <- "citeseq_v1"
base <- file.path(cfg$out_root, sprintf("benchmarks_hto_%d", bench))
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
ctp    <- file.path(base, "ctp") 
scl    <- file.path(base, "sclinear") 
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmle_ds_by_bench <- c("citeseq_v1_final")

bench_name <- c(`1` = "Kaggle PBMC")

bench_pick <- 1

read_one_bench <- function(bench, fmle_ds){

  fmle_path <- file.path(base, "FMLE", fmle_ds, "fmle_test_metrics.csv")
  scl_path  <- file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv")
  ctp_path  <- file.path(base, "ctp", "ctpnet_test_metrics_FMLEscale.csv")
  
  stopifnot(file.exists(fmle_path),
            file.exists(scl_path),
            file.exists(ctp_path))
  cat("\n[bench]", bench, "\n",
      "FMLE:", fmle_path, "\n",
      "SCL :", scl_path,  "\n",
      "CTP :", ctp_path,  "\n", sep="")
  ds <- bench_name[as.character(bench)]
  bind_rows(
    read_csv(fmle_path, show_col_types=FALSE) %>%
      transmute(dataset=ds,
                method="FMLE",
                protein,
                Pearson),
    
    read_csv(scl_path, show_col_types=FALSE) %>%
      transmute(dataset=ds,
                method="scLinear",
                protein,
                Pearson),
    
    read_csv(ctp_path, show_col_types=FALSE) %>%
      transmute(dataset=ds,
                method="cTPnet",
                protein,
                Pearson)
  )
}

df_kaggle <- read_one_bench(bench_pick, fmle_ds_by_bench[bench_pick]) %>%
  filter(is.finite(Pearson)) %>%
  mutate(
    dataset = factor(dataset, levels = c("Kaggle PBMC")),
    method  = factor(method,  levels = c("scLinear","cTPnet","FMLE"))
  )

pA_kaggle <- ggplot(df_kaggle, aes(x=method, y=Pearson, color=method)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  theme_classic(base_size = 15) +
  labs(
    x = NULL,
    y = "Per-protein Pearson correlation"
  ) +
  theme(legend.position = "none") +
  ggtitle("Kaggle PBMC") +
  theme(plot.title = element_text(face = "bold"))



results_fmle <- read.csv(
  file.path(base, "FMLE", "citeseq_v1_final", "fmle_test_metrics.csv"),
  stringsAsFactors = FALSE
)

results_sclinear <- read.csv(
  file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv"),
  stringsAsFactors = FALSE
)

results_ctpnet <- read.csv(
  file.path(base, "ctp", "ctpnet_test_metrics_FMLEscale.csv"),
  stringsAsFactors = FALSE
)

make_cmp_one_dataset <- function(results_fmle, res_sclinear, res_ctpnet=NULL) {
  stopifnot(all(c("protein","R2","Pearson") %in% colnames(results_fmle)))
  stopifnot(all(c("protein","R2","Pearson") %in% colnames(res_sclinear)))
  
  cmp <- inner_join(
    res_sclinear %>% select(protein, R2_scLinear = R2, Pearson_scLinear = Pearson),
    results_fmle %>% select(protein, R2_FMLE = R2, Pearson_FMLE = Pearson),
    by = "protein"
  )
  
  if (!is.null(res_ctpnet)) {
    cmp <- inner_join(
      cmp,
      res_ctpnet %>% select(protein, R2_cTPnet = R2, Pearson_cTPnet = Pearson),
      by = "protein"
    )
  }
  
  cmp
}


cmp <- make_cmp_one_dataset(results_fmle, results_sclinear, results_ctpnet)

mean(cmp$Pearson_FMLE - cmp$Pearson_scLinear)
mean(cmp$Pearson_FMLE - cmp$Pearson_cTPnet)

#----------------------------------------------------------------------------------#
# Per-protein gains
#----------------------------------------------------------------------------------#

cmp$delta_Pearson_scLinear <- cmp$Pearson_FMLE - cmp$Pearson_scLinear
cmp$delta_Pearson_cTPnet <- cmp$Pearson_FMLE - cmp$Pearson_cTPnet

summary(cmp$delta_Pearson_scLinear)
summary(cmp$delta_Pearson_cTPnet)

#----------------------------------------------------------------------------------#
# Statistical test 
#----------------------------------------------------------------------------------#
wilcox.test(cmp$Pearson_FMLE, cmp$Pearson_scLinear, paired=TRUE)
wilcox.test(cmp$Pearson_FMLE, cmp$Pearson_cTPnet, paired=TRUE)

#----------------------------------------------------------------------------------#
scl_ct <- readr::read_csv(file.path(scl, "scLinear_within_celltype_pearson.csv")) %>%
  dplyr::rename(Pearson_scl = Pearson)

fmle_ct <- readr::read_csv(file.path(out_dir, "fmle_within_celltype_pearson.csv"))

cmp_fmle_scl <- inner_join(fmle_ct, scl_ct, by=c("protein","cell_type"))

cmp_fmle_scl <- cmp_fmle_scl %>% dplyr::filter(cell_type != "Unassigned")

ctpnet_ct <- readr::read_csv(file.path(ctp, "ctpnet_within_celltype_pearson.csv")) %>%
  dplyr::rename(Pearson_ctp = Pearson)


cmp_fmle_ctpnet <- inner_join(fmle_ct,
                              ctpnet_ct,
                              by=c("protein","cell_type"))

lim <- range(
  c(cmp_fmle_scl$Pearson_scl,
    cmp_fmle_scl$Pearson_fmle,
    cmp_fmle_ctpnet$Pearson_ctp,
    cmp_fmle_ctpnet$Pearson_fmle),
  finite = TRUE
)
padx <- 0.1 * diff(lim)
pady <- 0.1 * diff(lim)   # more vertical space

xlim2 <- lim + c(-padx, padx)
ylim2 <- lim + c(-pady, pady)

p_fmle_scl <- ggplot(cmp_fmle_scl, aes(Pearson_scl, Pearson_fmle, color=protein, shape=cell_type)) +
  geom_point(size=3, alpha=0.9) +
  geom_abline(slope=1, intercept=0, color="grey40") +
  coord_cartesian(xlim = xlim2, ylim = ylim2, expand = FALSE)+
  theme_light(base_size=15) +
  labs(x="scLinear", y="FMLE", color="Protein", shape="Cell type")


p_fmle_ctp <- ggplot(cmp_fmle_ctpnet, aes(Pearson_ctp, Pearson_fmle, color=protein, shape=cell_type)) +
  geom_point(size=3, alpha=0.9) +
  geom_abline(slope=1, intercept=0, color="grey40") +
  coord_cartesian(xlim = xlim2, ylim = ylim2, expand = FALSE)  +
  theme_light(base_size=15) +
  labs(x="cTPnet", y="FMLE", color="Protein", shape="Cell type")

p_fmle_scl <- p_fmle_scl +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))
p_fmle_ctp <- p_fmle_ctp +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))



#==============================================
# 2) Cross-dataset transfer PBMC
#==============================================

base <- file.path(cfg$out_root, "transfer_preds_PBMC")
dir.create(base, recursive=TRUE, showWarnings=FALSE)

bench_name <- c(
  A = "Kaggle PBMC",
  B = "10x PBMC 10k",
  C = "TEA-seq PBMC"
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
  "A","C",
  "B","A",
  "B","C",
  "C","A",
  "C","B"
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
  facet_wrap(~ condition, ncol = 3) +
  theme_classic(base_size = 14) +
  labs(x = NULL, y = "Per-protein Pearson correlation") +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )

cond_levels <- unique(df_all$condition)

make_transfer_panel <- function(cond_lab) {
  ggplot(df_all %>% filter(condition == cond_lab),
         aes(x = method, y = Pearson, color = method)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, alpha = 0.7, size = 1.6) +
    theme_classic(base_size = 13) +
    labs(x = NULL, y = NULL, title = cond_lab) +
    theme(
      legend.position = "none",
      plot.title = element_text(face="bold", size=12),
      axis.text.x = element_text(angle=30, hjust=1)
    )
}

pT_list <- lapply(cond_levels, make_transfer_panel)


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
  mutate(better = factor(ifelse(gain > 0, "FMLE", "Baseline"),
                         levels = c("FMLE","Baseline")))


pC_scl <- ggplot(df_pair2 %>% filter(baseline_method=="scLinear"),
                 aes(x = baseline, y = FMLE)) +
  geom_point(aes(color = better),
             shape = 17, alpha = 0.85, size = 2.2) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  scale_color_manual(
    name   = "Better method",
    values = c("FMLE"="#1b9e77","Baseline"="#d95f02")
  ) +
  guides(color = guide_legend(ncol=1, byrow=TRUE)) +
  theme_classic(base_size=13) +
  labs(x=NULL, y=NULL, title="FMLE vs scLinear")  +
  theme(plot.title = element_text(face = "bold"))


pC_ctp <- ggplot(df_pair2 %>% filter(baseline_method=="cTPnet"),
                 aes(x = baseline, y = FMLE)) +
  geom_point(aes(color = better),
             shape = 17, alpha = 0.85, size = 2.2) +
  geom_abline(slope=1, intercept=0, linetype=2) +
  scale_color_manual(
    name   = "Better method",
    values = c("FMLE"="#1b9e77","Baseline"="#d95f02")
  ) +
  guides(color = guide_legend(ncol=1, byrow=TRUE)) +
  theme_classic(base_size=13) +
  labs(x=NULL, y=NULL, title="FMLE vs cTPnet") +
  theme(plot.title = element_text(face = "bold"))




heat_df <- df_all %>%
  group_by(method, train_ds, test_ds) %>%
  summarise(
    med = median(Pearson, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    transfer = paste0(train_ds, "\u2192", test_ds),
    transfer = factor(
      transfer,
      levels = c("A→B","A→C","B→A","B→C","C→A","C→B")
    ),
    method = factor(method, levels = c("FMLE","cTPnet","scLinear"))
  )



pB <- ggplot(heat_df, aes(x = transfer, y = method, fill = med)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", med)), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#d73027", mid = "#fee08b", high = "#1a9850",
    midpoint = 0.6, name = "Median\nPearson"
  ) +
  theme_classic(base_size = 12) +
  labs(x = "Train \u2192 Test", y = NULL) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )


pB_panel <- pB + theme(legend.position="none") + labs(title="Transfer heatmap") +
  theme(plot.title = element_text(face="bold", size=12))

pA_panel <- pA_kaggle + theme(legend.position="none") + labs(title="Kaggle PBMC") +
  theme(plot.title = element_text(face="bold", size=12))

pS_panel <- p_fmle_scl + theme(legend.position="none") + labs(title="FMLE vs scLinear (Kaggle)") +
  theme(plot.title = element_text(face="bold", size=12))

pC_panel <- p_fmle_ctp + theme(legend.position="none") + labs(title="FMLE vs cTPnet (Kaggle)") +
  theme(plot.title = element_text(face="bold", size=12))


dataset_key <- "A = Kaggle PBMC\nB = 10x PBMC 10k\nC = TEA-seq PBMC"

pA_legend_source <- pA_kaggle +
  guides(color = guide_legend(title="Method", ncol=1, byrow=TRUE)) +
  theme(legend.position = "right", legend.box = "vertical")

pA_legend_source <- pA_legend_source +
  geom_point(
    data = data.frame(key = dataset_key),
    mapping = aes(x = Inf, y = Inf, shape = key),
    inherit.aes = FALSE,
    alpha = 0,
    show.legend = TRUE
  ) +
  scale_shape_manual(
    name = NULL,
    values = setNames(NA, dataset_key),
    guide = guide_legend(override.aes = list(alpha = 0))
  )

pA_legend_source <- pA_legend_source +
  scale_shape_manual(
    name = NULL,
    values = setNames(NA, dataset_key),
    guide = guide_legend(
      override.aes = list(alpha = 0),
      label.position = "left",   # puts text left of key box
      label.hjust = 0            # left-justify the text
    )
  ) +
  theme(
    legend.text.align = 0,                 # left align legend text
    legend.box.margin = margin(0, 0, 0, -6), # shift whole legend block left (tune -6)
    legend.spacing.x = unit(0, "pt"),
    legend.key.width = unit(0, "pt")
  )

scale_method <- scale_color_manual(values = c(
  "FMLE" = "#1b9e77",
  "scLinear" = "#d95f02",
  "cTPnet" = "#7570b3"
))


pT_list <- lapply(pT_list, \(p)
                  p + scale_method + guides(color = "none", fill = "none") + theme(legend.position = "none")
)

pS_panel <- pS_panel +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.box = "vertical")

pC_panel <- pC_panel +
  guides(color = guide_legend(ncol = 1, byrow = TRUE)) +
  theme(legend.box = "vertical")

grid12 <- wrap_plots(
  c(
    pT_list[1], pT_list[2], pT_list[3], pA_legend_source,
    pT_list[4], pT_list[5], pT_list[6], pS_panel,
    pB_panel,   pC_scl,     pC_ctp,     pC_panel
  ),
  ncol = 4
) +
  plot_layout(guides = "collect", heights = rep(1,3), widths = c(1,1,1,1)) +
  plot_annotation(
    tag_levels = list(c("a","b","c","j","d","e","f","k","g","h","i","l")),
    tag_prefix = "(",
    tag_suffix = ")"
  ) &
  theme(
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text  = element_text(size = 10),
    plot.tag.position = c(0.05, 0.95),
    plot.tag = element_text(face = "bold", size = 12)
  )
grid12 <- grid12 &
  theme(plot.title = element_text(margin = margin(t = 8)))

gap <- 0.0335   

# left box: 
x_left   <- (0.721 - gap) / 2
w_left   <-  0.721 - gap

# right box: 
gap1 <- 0.055
x_right  <- 0.652 + (0.24 + gap1)/2
w_right  <- 0.275 - gap1

box_h <- 0.97            
y_box <- 0.50           

# compute top edge as a grid unit
y_top_u <- unit(y_box + box_h/2, "npc")

final_fig <- wrap_elements(full =
                             grobTree(
                               patchworkGrob(grid12),
                               
                               rectGrob(
                                 x = unit(x_left, "npc"), y = unit(y_box, "npc"),
                                 width = unit(w_left, "npc"), height = unit(box_h, "npc"),
                                 gp = gpar(col="black", fill=NA, lty=2, lwd=1.2)
                               ),
                               
                               rectGrob(
                                 x = unit(x_right, "npc"), y = unit(y_box, "npc"),
                                 width = unit(w_right, "npc"), height = unit(box_h, "npc"),
                                 gp = gpar(col="black", fill=NA, lty=2, lwd=1.2)
                               ),
                               
                               # labels (no clipping)
                               grobTree(
                                 textGrob(
                                   "Cross-dataset transfer",
                                   x = unit(x_left, "npc"),
                                   y = y_top_u + unit(1.0, "mm"),   # <- above border, but in mm not npc
                                   just = c("center","bottom"),
                                   gp = gpar(fontface="bold", fontsize=16)
                                 ),
                                 textGrob(
                                   "Cross-donor transfer",
                                   x = unit(x_right, "npc"),
                                   y = y_top_u + unit(1.0, "mm"),
                                   just = c("center","bottom"),
                                   gp = gpar(fontface="bold", fontsize=16)
                                 ),
                                 vp = viewport(clip = "off")        # <- CRITICAL
                               )
                             )
)



final


ggsave(file.path(out_dir, "Figure_3.pdf"),
       plot = final,
       device = cairo_pdf,
       width = 19, height = 13.2, units = "in",
       dpi = 300,
       limitsize = FALSE)

# =========================================================
# Comparision with baseline method in cross data transfer
# =========================================================

cmp_all <- df_all %>%
  select(condition, protein, method, Pearson) %>%
  tidyr::pivot_wider(names_from = method, values_from = Pearson) %>%
  filter(is.finite(FMLE), is.finite(scLinear), is.finite(cTPnet)) %>%
  mutate(
    delta_scLinear = FMLE - scLinear,
    delta_cTPnet   = FMLE - cTPnet
  )


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





