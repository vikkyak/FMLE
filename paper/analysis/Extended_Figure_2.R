suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(grid)
})

# ---- Dataset 1: BMMC CITE-seq
bench1 <- 1
ds <- "citeseq_v1"

base1 <- file.path(cfg$out_root, sprintf("benchmarks_%d", bench1))
out_dir1 <- file.path(base1, "FMLE", paste0(ds, "_final"))
ctp1 <- file.path(base1, "ctp")
scl1 <- file.path(base1, "sclinear")
fig_dir1 <- file.path(base1, "paper_figures", ds)

dir.create(out_dir1, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir1, recursive = TRUE, showWarnings = FALSE)

fig2_out1 <- file.path(fig_dir1, "fig_examples_same.rds") # from validation function
fig2_1 <- readRDS(fig2_out1)

seu_prep1 <- readRDS(file.path(base1, ds, "seu_final.rds"))
stopifnot("cell_type" %in% colnames(seu_prep1@meta.data))
celltype_vec1 <- setNames(as.character(seu_prep1@meta.data$cell_type), colnames(seu_prep1))

# ---- Dataset 2: Young healthy BM
bench2 <- 2
base2 <- file.path(cfg$out_root, sprintf("benchmarks_%d", bench2))
out_dir2 <- file.path(base2, "FMLE", paste0(ds, "_final"))
ctp2 <- file.path(base2, "ctp")
scl2 <- file.path(base2, "sclinear")
fig_dir2 <- file.path(base2, "paper_figures", ds)

dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir2, recursive = TRUE, showWarnings = FALSE)

fig2_out2 <- file.path(fig_dir2, "fig_examples_same.rds")   # from validation function
fig2_2 <- readRDS(fig2_out2)

seu_prep2 <- readRDS(file.path(base2, ds, "seu_final.rds"))
stopifnot("predicted.celltype.l1" %in% colnames(seu_prep2@meta.data))
celltype_vec2 <- setNames(as.character(seu_prep2@meta.data$predicted.celltype.l1), colnames(seu_prep2))



# ============================================================
# Overall benchmark summary
# ============================================================
fmle_ds_by_bench <- c(
  "citeseq_v1_final",
  "citeseq_v1_final",
  "citeseq_v1_final"
)

bench_name <- c(
  `1` = "BMMC CITE-seq",
  `2` = "Young healthy BM"
)

read_one_bench <- function(bench, fmle_ds){
  
  base <- file.path(out_root, sprintf("benchmarks_%d", bench))
  fmle_path <- file.path(base, "FMLE", fmle_ds, "fmle_test_metrics.csv")
  scl_path  <- file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv")
  ctp_path  <- file.path(base, "ctp", "ctpnet_test_metrics_FMLEscale.csv")
  
  stopifnot(file.exists(fmle_path), file.exists(scl_path), file.exists(ctp_path))
  
  ds_name <- bench_name[as.character(bench)]
  
  bind_rows(
    read_csv(fmle_path, show_col_types = FALSE) %>%
      transmute(dataset = ds_name, method = "FMLE", protein, Pearson),
    read_csv(scl_path, show_col_types = FALSE) %>%
      transmute(dataset = ds_name, method = "scLinear", protein, Pearson),
    read_csv(ctp_path, show_col_types = FALSE) %>%
      transmute(dataset = ds_name, method = "cTPnet", protein, Pearson)
  )
}

# -------- BUILD DATA --------
df_all <- bind_rows(
  read_one_bench(1, fmle_ds_by_bench[1]),
  read_one_bench(2, fmle_ds_by_bench[2])
) %>%
  filter(is.finite(Pearson)) %>%
  mutate(
    dataset = factor(dataset, levels = c("BMMC CITE-seq", "Young healthy BM")),
    method  = factor(method, levels = c("scLinear", "cTPnet", "FMLE"))
  )


# ===============================================================
# UPPER PANEL A "BMMC CITE-seq", UPPER PANEL D "Young healthy BM" 
# ===============================================================
pA1 <- ggplot(
  df_all %>% filter(dataset == "BMMC CITE-seq"),
  aes(x = method, y = Pearson, color = method)
) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.5) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.5) +
  facet_wrap(~dataset, ncol = 1) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Per-protein Pearson correlation"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.6),
    axis.text.x = element_text(size = 11),
    plot.margin = margin(1, 1, 1, 1)
  )

pA2 <- ggplot(
  df_all %>% filter(dataset == "Young healthy BM"),
  aes(x = method, y = Pearson, color = method)
) +
  geom_boxplot(outlier.shape = NA, width = 0.58, linewidth = 0.5) +
  geom_jitter(width = 0.10, alpha = 0.55, size = 1.5) +
  facet_wrap(~dataset, ncol=1) +
  theme_classic(base_size = 13) +
  labs(
    x=NULL,
    y="Per-protein Pearson correlation"
  )+
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.6),
    axis.text.x = element_text(size = 11),
    plot.margin = margin(1, 1, 1, 1)
  )

# ===============================================================
# Panel A and D — FMLE vs scLinear and FMLE vs cTPnet per protien
# ===============================================================

# ---- Dataset 1------  "BMMC CITE-seq" ----  Panel A

scl_ct1_tbl <- read_csv(file.path(scl1, "scLinear_within_celltype_pearson.csv"),
                        show_col_types = FALSE) %>%
  rename(Pearson_scl = Pearson)

fmle_ct1_tbl <- read_csv(file.path(out_dir1, "fmle_within_celltype_pearson.csv"),
                         show_col_types = FALSE)

ctpnet_ct1_tbl <- read_csv(file.path(ctp1, "ctpnet_within_celltype_pearson.csv"),
                           show_col_types = FALSE) %>%
  rename(Pearson_ctp = Pearson)

cmp_fmle_scl_1 <- inner_join(fmle_ct1_tbl, scl_ct1_tbl, by = c("protein", "cell_type")) %>%
  filter(
    cell_type != "Unassigned",
    is.finite(Pearson_fmle),
    is.finite(Pearson_scl)
  )

cmp_fmle_ctp_1 <- inner_join(fmle_ct1_tbl, ctpnet_ct1_tbl, by = c("protein", "cell_type")) %>%
  filter(
    cell_type != "Unassigned",
    is.finite(Pearson_fmle),
    is.finite(Pearson_ctp)
  )

# ---- Dataset 2------     "Young healthy BM"------Panel D

scl_ct2_tbl <- read_csv(file.path(scl2, "scLinear_within_celltype_pearson.csv"),
                        show_col_types = FALSE) %>%
  rename(Pearson_scl = Pearson)

fmle_ct2_tbl <- read_csv(file.path(out_dir2, "fmle_within_celltype_pearson.csv"),
                         show_col_types = FALSE)

ctpnet_ct2_tbl <- read_csv(file.path(ctp2, "ctpnet_within_celltype_pearson.csv"),
                           show_col_types = FALSE) %>%
  rename(Pearson_ctp = Pearson)

cmp_fmle_scl_2 <- inner_join(fmle_ct2_tbl, scl_ct2_tbl, by = c("protein", "cell_type")) %>%
  filter(
    cell_type != "Unassigned",
    is.finite(Pearson_fmle),
    is.finite(Pearson_scl)
  )

cmp_fmle_ctp_2 <- inner_join(fmle_ct2_tbl, ctpnet_ct2_tbl, by = c("protein", "cell_type")) %>%
  filter(
    cell_type != "Unassigned",
    is.finite(Pearson_fmle),
    is.finite(Pearson_ctp)
  )


# ==========================================================
# panel B, C and E, F
# ==========================================================

get_shared_lim <- function(prot, fig2) {
  vals <- c(
    as.numeric(fig2$truth[[prot]]),
    as.numeric(fig2$pred_fmle[[prot]]),
    as.numeric(fig2$pred_sclinear[[prot]]),
    as.numeric(fig2$pred_ctpnet[[prot]])
  )
  vals <- vals[is.finite(vals)]
  rng <- range(vals)
  pad <- 0.04 * diff(rng)
  rng + c(-pad, pad)
}

plot_panel_clean <- function(prot,
                             method = c("FMLE", "scLinear", "cTPnet"),
                             fig2, celltype_vec,
                             n_points = 4000) {
  method <- match.arg(method)
  
  pred_list <- switch(
    method,
    "FMLE"     = fig2$pred_fmle,
    "scLinear" = fig2$pred_sclinear,
    "cTPnet"   = fig2$pred_ctpnet
  )
  
  stopifnot(
    prot %in% names(fig2$truth),
    prot %in% names(pred_list),
    prot %in% names(fig2$cells)
  )
  
  df <- data.frame(
    pred = as.numeric(pred_list[[prot]]),
    real = as.numeric(fig2$truth[[prot]]),
    cell = fig2$cells[[prot]],
    stringsAsFactors = FALSE
  )
  
  df <- df[is.finite(df$pred) & is.finite(df$real), , drop = FALSE]
  df$celltype <- unname(celltype_vec[df$cell])
  df$celltype[is.na(df$celltype)] <- "Unassigned"
  df <- df[df$celltype != "Unassigned", , drop = FALSE]
  df$celltype <- droplevels(factor(df$celltype))
  
  if (nrow(df) > n_points) {
    set.seed(1)
    df <- df[sample.int(nrow(df), n_points), , drop = FALSE]
  }
  
  r   <- suppressWarnings(cor(df$pred, df$real, method = "pearson"))
  rho <- suppressWarnings(cor(df$pred, df$real, method = "spearman"))
  
  lims <- get_shared_lim(prot, fig2)
  
  title_txt <- sprintf("%s (r = %.2f, ρ = %.2f)", prot, r, rho)
  
  ggplot(df, aes(x = pred, y = real, color = celltype)) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "grey35", linewidth = 0.6
    ) +
    geom_point(alpha = 0.55, size = 1.0) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE) +
    theme_classic(base_size = 12) +
    labs(
      title = title_txt,
      x = paste0(method, " predicted"),
      y = "True protein",
      color = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      aspect.ratio = 1,
      plot.margin = margin(1, 1, 1, 1),
      legend.position = "none"
    )
}

make_protein_block_clean <- function(prot, fig2, celltype_vec) {
  p1 <- plot_panel_clean(prot, "FMLE", fig2, celltype_vec)
  p2 <- plot_panel_clean(prot, "scLinear", fig2, celltype_vec)
  p3 <- plot_panel_clean(prot, "cTPnet", fig2, celltype_vec)
  p1 / p2 / p3
}

# ---- Dataset 1 representative proteins

p_CD11c_1 <- make_protein_block_clean("CD11c", fig2_1, celltype_vec1)

p_CD40_1  <- make_protein_block_clean("CD40",  fig2_1, celltype_vec1)

# ---- Dataset 2 representative proteins

p_CD11c_2 <- make_protein_block_clean("CD11c-AB", fig2_2, celltype_vec2)
p_CD14_2  <- make_protein_block_clean("CD14-AB",  fig2_2, celltype_vec2)

# ============================================================
# Curated summary subset for dataset 1
# ============================================================
proteins_show_1 <- c(
  "CD71", "CD94", "CD36", "CD19", "CD56",
  "CD16", "CD14", "IgM", "CD57", "CD8"
)

celltypes_show_1 <- c(
  "Naive CD20+ B IGKC+",
  "CD14+ Mono",
  "CD16+ Mono",
  "CD4+ T naive",
  "CD8+ T naive",
  "NK",
  "pDC",
  "Erythroblast"
)

sub_scl_1 <- cmp_fmle_scl_1 %>%
  filter(protein %in% proteins_show_1, cell_type %in% celltypes_show_1) %>%
  filter(is.finite(Pearson_scl), is.finite(Pearson_fmle)) %>%
  mutate(
    protein   = factor(protein, levels = proteins_show_1),
    cell_type = factor(cell_type, levels = celltypes_show_1)
  )

sub_ctp_1 <- cmp_fmle_ctp_1 %>%
  filter(protein %in% proteins_show_1, cell_type %in% celltypes_show_1) %>%
  filter(is.finite(Pearson_ctp), is.finite(Pearson_fmle)) %>%
  mutate(
    protein   = factor(protein, levels = proteins_show_1),
    cell_type = factor(cell_type, levels = celltypes_show_1)
  )

stopifnot(nrow(sub_scl_1) > 0, nrow(sub_ctp_1) > 0)

lim_mid_1 <- range(
  c(sub_scl_1$Pearson_scl, sub_scl_1$Pearson_fmle,
    sub_ctp_1$Pearson_ctp, sub_ctp_1$Pearson_fmle),
  finite = TRUE
)
pad_mid_1 <- 0.06 * diff(lim_mid_1)
xylim_mid_1 <- lim_mid_1 + c(-pad_mid_1, pad_mid_1)

protein_cols_1 <- c(
  "CD71" = "#D73027",
  "CD94" = "#7B3294",
  "CD36" = "#E08214",
  "CD19" = "#1B9E77",
  "CD56" = "#66A61E",
  "CD16" = "#4DBBD5",
  "CD14" = "#E6AB02",
  "IgM"  = "#A6761D",
  "CD57" = "#7570B3",
  "CD8"  = "#1F78B4"
)

celltype_shapes_1 <- c(
  "Naive CD20+ B IGKC+" = 15,
  "CD14+ Mono"         = 17,
  "CD16+ Mono"         = 8,
  "CD4+ T naive"       = 3,
  "CD8+ T naive"       = 16,
  "NK"                 = 18,
  "pDC"                = 4,
  "Erythroblast"       = 7
)

# ============================================================
# Curated summary subset for dataset 2
# ============================================================
proteins_show_2 <- c(
  "CD1c-AB", "CD56-AB", "CD14-AB", "CD9-AB", "CD32-AB",
  "CD226-AB", "CD155-AB", "CD163-AB", "IgD-AB", "CD39-AB"
)

celltypes_show_2 <- c(
  "B", "CD4 T", "CD8 T", "Mono",
  "NK", "DC", "HSPC", "other T"
)


sub_scl_2 <- cmp_fmle_scl_2 %>%
  filter(protein %in% proteins_show_2, cell_type %in% celltypes_show_2) %>%
  filter(is.finite(Pearson_scl), is.finite(Pearson_fmle)) %>%
  mutate(
    protein   = factor(protein, levels = proteins_show_2),
    cell_type = factor(cell_type, levels = celltypes_show_2)
  )

sub_ctp_2 <- cmp_fmle_ctp_2 %>%
  filter(protein %in% proteins_show_2, cell_type %in% celltypes_show_2) %>%
  filter(is.finite(Pearson_ctp), is.finite(Pearson_fmle)) %>%
  mutate(
    protein   = factor(protein, levels = proteins_show_2),
    cell_type = factor(cell_type, levels = celltypes_show_2)
  )

stopifnot(nrow(sub_scl_2) > 0, nrow(sub_ctp_2) > 0)

lim_mid_2 <- range(
  c(sub_scl_2$Pearson_scl, sub_scl_2$Pearson_fmle,
    sub_ctp_2$Pearson_ctp, sub_ctp_2$Pearson_fmle),
  finite = TRUE
)
pad_mid_2 <- 0.06 * diff(lim_mid_2)
xylim_mid_2 <- lim_mid_2 + c(-pad_mid_2, pad_mid_2)

protein_cols_2 <- c(
  "CD1c-AB"   = "#D73027",
  "CD56-AB"   = "#1B9E77",
  "CD14-AB"   = "#E6AB02",
  "CD9-AB"    = "#7B3294",
  "CD32-AB"   = "#4DBBD5",
  "CD226-AB"  = "#E08214",
  "CD155-AB"  = "#A6761D",
  "CD163-AB"  = "#7570B3",
  "IgD-AB"    = "#66A61E",
  "CD39-AB"   = "#1F78B4"
)

celltype_shapes_2 <- c(
  "B"       = 15,
  "CD4 T"   = 3,
  "CD8 T"   = 16,
  "Mono"    = 17,
  "NK"      = 18,
  "DC"      = 4,
  "HSPC"    = 8,
  "other T" = 7
)
# ============================================================
# Generic summary panel function
# ============================================================
make_summary_panel <- function(df, xvar, xlab, protein_cols, celltype_shapes, xylim_use) {
  ggplot(
    df,
    aes(x = .data[[xvar]], y = Pearson_fmle, color = protein, shape = cell_type)
  ) +
    geom_abline(slope = 1, intercept = 0, color = "grey65", linewidth = 0.7) +
    geom_point(size = 2.8, alpha = 0.88) +
    scale_color_manual(values = protein_cols, drop = FALSE) +
    scale_shape_manual(values = celltype_shapes, drop = FALSE) +
    coord_cartesian(xlim = xylim_use, ylim = xylim_use, expand = FALSE) +
    theme_light(base_size = 13) +
    labs(
      x = xlab,
      y = "FMLE",
      color = "Protein",
      shape = "Cell type"
    ) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.40, "lines"),
      plot.margin = margin(1, 1, 1, 1)
    ) +
    guides(
      color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
      shape = guide_legend(order = 2, nrow = 2, byrow = TRUE)
    )
}
# ---- Dataset 1 middle block
p_fmle_scl_ext_1 <- make_summary_panel(
  sub_scl_1, "Pearson_scl", "scLinear",
  protein_cols = protein_cols_1,
  celltype_shapes = celltype_shapes_1,
  xylim_use = xylim_mid_1
)

p_fmle_ctp_ext_1 <- make_summary_panel(
  sub_ctp_1, "Pearson_ctp", "cTPnet",
  protein_cols = protein_cols_1,
  celltype_shapes = celltype_shapes_1,
  xylim_use = xylim_mid_1
)

# ---- Dataset 2 middle block
p_fmle_scl_ext_2 <- make_summary_panel(
  sub_scl_2, "Pearson_scl", "scLinear",
  protein_cols = protein_cols_2,
  celltype_shapes = celltype_shapes_2,
  xylim_use = xylim_mid_2
)

p_fmle_ctp_ext_2 <- make_summary_panel(
  sub_ctp_2, "Pearson_ctp", "cTPnet",
  protein_cols = protein_cols_2,
  celltype_shapes = celltype_shapes_2,
  xylim_use = xylim_mid_2
)

legend_plot_1 <- p_fmle_scl_ext_1
legend_plot_2 <- p_fmle_scl_ext_2

leg1 <- cowplot::get_legend(legend_plot_1)
leg2 <- cowplot::get_legend(legend_plot_2)

p_fmle_scl_ext_1 <- p_fmle_scl_ext_1 + theme(legend.position = "none")
p_fmle_ctp_ext_1 <- p_fmle_ctp_ext_1 + theme(legend.position = "none")

p_fmle_scl_ext_2 <- p_fmle_scl_ext_2 + theme(legend.position = "none")
p_fmle_ctp_ext_2 <- p_fmle_ctp_ext_2 + theme(legend.position = "none")

mid_block_1 <- p_fmle_scl_ext_1 / p_fmle_ctp_ext_1
mid_block_2 <- p_fmle_scl_ext_2 / p_fmle_ctp_ext_2

left_block_1  <- pA1
left_block_2  <- pA2

left_stack_1 <- (left_block_1 / mid_block_1) +
  plot_layout(heights = c(0.55, 1.65))

left_stack_2 <- (left_block_2 / mid_block_2) +
  plot_layout(heights = c(0.55, 1.65))

right_block_1 <- p_CD11c_1 | p_CD40_1
right_block_2 <- p_CD11c_2 | p_CD14_2

top_block_1 <- (
  wrap_elements(full = left_stack_1) |
    wrap_elements(full = right_block_1)
) + plot_layout(widths = c(1.6, 2.2))

top_block_2 <- (
  wrap_elements(full = left_stack_2) |
    wrap_elements(full = right_block_2)
) + plot_layout(widths = c(1.6, 2.2))



final_clean_1 <- (
  wrap_elements(full = top_block_1) /
    wrap_elements(full = leg1)
) + plot_layout(heights = c(1, 0.06))

final_clean_2 <- (
  wrap_elements(full = top_block_2) /
    wrap_elements(full = leg2)
) + plot_layout(heights = c(1, 0.06))

final <- (
  wrap_elements(full = final_clean_1) |
    wrap_elements(full = final_clean_2)
) + plot_layout(heights = c(1, 1))

final 


ggsave(file.path(out_dir, "Figure_2_extended.pdf"),
       plot = final,
       device = cairo_pdf,
       width = 21, height = 10, units = "in",
       dpi = 300,
       limitsize = FALSE)















