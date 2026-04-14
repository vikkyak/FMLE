suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

ds <- "citeseq_v1"

# ----------------------------
# donor-specific base folders
# ----------------------------
# Donor NS1
base1 <- file.path(cfg$out_root, ds, "LODO_1_test_NS1")
dir.create(base, recursive=TRUE, showWarnings=FALSE)
# Donor TP7
base2 <- file.path(cfg$out_root, ds, "LODO_2_test_TP7")
dir.create(base, recursive=TRUE, showWarnings=FALSE)
# Donor TS5
base3 <- file.path(cfg$out_root, ds, "LODO_3_test_TS5")
dir.create(base, recursive=TRUE, showWarnings=FALSE)


bench_base <- c(
  `1` = base1,
  `2` = base2,
  `3` = base3
)

bench_name <- c(
  `1` = "NS1",
  `2` = "TP7",
  `3` = "TS5"
)

fmle_ds_by_bench <- c(
  `1` = "citeseq_v1_final",
  `2` = "citeseq_v1_final",
  `3` = "citeseq_v1_final"
)

# ----------------------------
# method folders
# ----------------------------
# cTPnet
# Donor NS1
ctp1 <- file.path(base1, "ctp")
dir.create(ctp1, recursive=TRUE, showWarnings=FALSE)
# Donor TP7
ctp2 <- file.path(base2, "ctp")
dir.create(ctp2, recursive=TRUE, showWarnings=FALSE)
# Donor TS5
ctp3 <- file.path(base3, "ctp")
dir.create(ctp3, recursive=TRUE, showWarnings=FALSE)

ctp_base <- c(
  `1` = ctp1,
  `2` = ctp2,
  `3` = ctp3
)

# scLinear
scl1 <- file.path(base1, "sclinear")
dir.create(scl1, recursive=TRUE, showWarnings=FALSE)
scl2 <- file.path(base2, "sclinear")
dir.create(scl2, recursive=TRUE, showWarnings=FALSE)
scl3 <- file.path(base3, "sclinear")
dir.create(scl3, recursive=TRUE, showWarnings=FALSE)

# FMLE

fmle1 <- file.path(base1, "FMLE", "citeseq_v1_final")
fmle2 <- file.path(base2, "FMLE", "citeseq_v1_final")
fmle3 <- file.path(base3, "FMLE", "citeseq_v1_final")

# ----------------------------
# figure example files
# ----------------------------
fig_dir1 <- file.path(base1, "paper_figures", "citeseq_v1")
fig_dir2 <- file.path(base2, "paper_figures", "citeseq_v1")
fig_dir3 <- file.path(base3, "paper_figures", "citeseq_v1")

fig2_out1 <- file.path(fig_dir1, "fig_examples_same.rds")   # from validation function
fig2_out2 <- file.path(fig_dir2, "fig_examples_same.rds")   # from validation function
fig2_out3 <- file.path(fig_dir3, "fig_examples_same.rds")   # from validation function

fig21 <- readRDS(fig2_out1)
fig22 <- readRDS(fig2_out2)
fig23 <- readRDS(fig2_out3)

cat("Available proteins NS1:\n")
print(names(fig21$truth))

cat("Available proteins TP7:\n")
print(names(fig22$truth))

cat("Available proteins TS5:\n")
print(names(fig23$truth))

# ----------------------------
# read Seurat objects
# ----------------------------
seu_prep1 <- readRDS(file.path(base1, "seu_final.rds"))
stopifnot("predicted.celltype.l1" %in% colnames(seu_prep1@meta.data))
celltype_vec1 <- setNames(
  as.character(seu_prep1@meta.data$predicted.celltype.l1),
  colnames(seu_prep1)
)

seu_prep2 <- readRDS(file.path(base2, "seu_final.rds"))
stopifnot("predicted.celltype.l1" %in% colnames(seu_prep2@meta.data))
celltype_vec2 <- setNames(
  as.character(seu_prep2@meta.data$predicted.celltype.l1),
  colnames(seu_prep2)
)

seu_prep3 <- readRDS(file.path(base3, "seu_final.rds"))
stopifnot("predicted.celltype.l1" %in% colnames(seu_prep3@meta.data))
celltype_vec3 <- setNames(
  as.character(seu_prep3@meta.data$predicted.celltype.l1),
  colnames(seu_prep3)
)

# ----------------------------
# read one donor
# ----------------------------
read_one_bench <- function(bench, fmle_ds) {
  
  base <- unname(bench_base[as.character(bench)])
  dsnm <- unname(bench_name[as.character(bench)])
  ctp  <- unname(ctp_base[as.character(bench)])
  
  fmle_path <- file.path(base, "FMLE", fmle_ds, "fmle_test_metrics.csv")
  scl_path  <- file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv")
  ctp_path  <- file.path(ctp,  "ctpnet_test_metrics_FMLEscale.csv")
  
  stopifnot(
    file.exists(fmle_path),
    file.exists(scl_path),
    file.exists(ctp_path)
  )
  
  bind_rows(
    read_csv(fmle_path, show_col_types = FALSE) %>%
      transmute(
        dataset = dsnm,
        method  = "FMLE",
        protein,
        Pearson
      ),
    
    read_csv(scl_path, show_col_types = FALSE) %>%
      transmute(
        dataset = dsnm,
        method  = "scLinear",
        protein,
        Pearson
      ),
    
    read_csv(ctp_path, show_col_types = FALSE) %>%
      transmute(
        dataset = dsnm,
        method  = "cTPnet",
        protein,
        Pearson
      )
  )
}

# ----------------------------
# build combined table
# ----------------------------
df_all <- bind_rows(
  read_one_bench(1, fmle_ds_by_bench["1"]),
  read_one_bench(2, fmle_ds_by_bench["2"]),
  read_one_bench(3, fmle_ds_by_bench["3"])
) %>%
  filter(is.finite(Pearson)) %>%
  mutate(
    dataset = factor(dataset, levels = c("NS1", "TP7", "TS5")),
    method  = factor(method, levels = c("scLinear", "cTPnet", "FMLE"))
  )

print(table(df_all$dataset, df_all$method))
head(df_all)

# -------- PANEL A PLOT --------
pA <- ggplot(df_all,
             aes(x=method, y=Pearson, color=method)) +
  
  geom_boxplot(outlier.shape=NA, width=0.6) +
  
  geom_jitter(width=0.15, alpha=0.7, size=2) +
  
  facet_wrap(~dataset, ncol=1) +
  
  theme_classic(base_size=12) +
  
  labs(
    x=NULL,
    y="Per-protein Pearson correlation"
  ) +
  
  theme(
    legend.position="none",
    strip.text=element_text(face="bold")
  )



# Panel B — FMLE vs scLinear per protien

scl_ct <- readr::read_csv(
  path.expand(scl3, "scLinear_within_celltype_pearson.csv"),
  show_col_types = FALSE
) %>%
  dplyr::rename(Pearson_scl = Pearson)

fmle_ct <- readr::read_csv(
  path.expand(fmle3, "fmle_within_celltype_pearson.csv"),
  show_col_types = FALSE
)

metric_col_fmle <- intersect(c("Pearson_fmle", "Pearson", "pearson", "cor"), names(fmle_ct))
stopifnot(length(metric_col_fmle) == 1)

fmle_ct <- fmle_ct %>%
  dplyr::rename(Pearson_fmle = dplyr::all_of(metric_col_fmle))

ctpnet_ct <- readr::read_csv(
  path.expand(ctp3, "ctpnet_within_celltype_pearson.csv"),
  show_col_types = FALSE
) %>%
  dplyr::rename(Pearson_ctp = Pearson)

cmp_fmle_scl <- dplyr::inner_join(
  fmle_ct,
  scl_ct,
  by = c("protein", "cell_type")
) %>%
  dplyr::filter(cell_type != "Unassigned")

cmp_fmle_ctpnet <- dplyr::inner_join(
  fmle_ct,
  ctpnet_ct,
  by = c("protein", "cell_type")
) %>%
  dplyr::filter(cell_type != "Unassigned")



proteins_show <- c(
  "anti-human-CD101-totalC",
  "anti-human-CD18-totalC",
  "anti-human-CD58-totalC",
  "anti-human-CD27-totalC",
  "anti-human-CD62P-totalC",
  "anti-human-CD7-totalC",
  "anti-human-CD45RA-totalC",
  "anti-human-CD71-totalC",
  "anti-human-CD95-totalC",
  "anti-human-HLA-DR-totalC"
)

protein_labels <- c(
  "anti-human-CD101-totalC"    = "CD101",
  "anti-human-CD18-totalC"     = "CD18",
  "anti-human-CD58-totalC"     = "CD58",
  "anti-human-CD27-totalC"     = "CD27",
  "anti-human-CD62P-totalC"    = "CD62P",
  "anti-human-CD7-totalC"      = "CD7",
  "anti-human-CD45RA-totalC"   = "CD45RA",
  "anti-human-CD71-totalC"     = "CD71",
  "anti-human-CD95-totalC"     = "CD95",
  "anti-human-HLA-DR-totalC"   = "HLA-DR"
)
celltypes_show <- c("B", "CD4 T", "CD8 T", "Mono", "NK")

celltype_shapes <- c(
  "B"     = 16,  # filled circle
  "CD4 T" = 17,  # filled triangle
  "CD8 T" = 15,  # filled square
  "Mono"  = 18,  # filled diamond
  "NK"    = 1    # open circle
)

method_cols <- c(
  "scLinear" = "#F8766D",
  "cTPnet"   = "#00BA38",
  "FMLE"     = "#619CFF"
)

sub_scl <- cmp_fmle_scl %>%
  dplyr::filter(protein %in% proteins_show, cell_type %in% celltypes_show) %>%
  dplyr::filter(is.finite(Pearson_scl), is.finite(Pearson_fmle)) %>%
  dplyr::mutate(
    protein   = factor(protein, levels = proteins_show, labels = protein_labels[proteins_show]),
    cell_type = factor(cell_type, levels = celltypes_show)
  )

sub_ctp <- cmp_fmle_ctpnet %>%
  dplyr::filter(protein %in% proteins_show, cell_type %in% celltypes_show) %>%
  dplyr::filter(is.finite(Pearson_ctp), is.finite(Pearson_fmle)) %>%
  dplyr::mutate(
    protein   = factor(protein, levels = proteins_show, labels = protein_labels[proteins_show]),
    cell_type = factor(cell_type, levels = celltypes_show)
  )

stopifnot(nrow(sub_scl) > 0, nrow(sub_ctp) > 0)

lim <- range(
  c(
    sub_scl$Pearson_scl,
    sub_scl$Pearson_fmle,
    sub_ctp$Pearson_ctp,
    sub_ctp$Pearson_fmle
  ),
  finite = TRUE
)

padx <- 0.06 * diff(lim)
pady <- 0.06 * diff(lim)

xlim2 <- lim + c(-padx, padx)
ylim2 <- lim + c(-pady, pady)


p_fmle_scl <- ggplot(
  sub_scl,
  aes(Pearson_scl, Pearson_fmle, color = protein, shape = cell_type)
) +
  geom_point(size = 2.1, alpha = 0.9, stroke = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.35) +
  coord_cartesian(xlim = xlim2, ylim = ylim2, expand = FALSE) +
  scale_shape_manual(values = celltype_shapes, drop = FALSE) +
  guides(
    shape = guide_legend(override.aes = list(size = 2.8, alpha = 1)),
    color = guide_legend(override.aes = list(size = 3.0, alpha = 1))
  ) +
  theme_classic(base_size = 12) +
  labs(x = "scLinear", y = "FMLE", color = "Protein", shape = "Cell type") +
  theme(
    axis.title = element_text(size = 8.2),
    axis.text  = element_text(size = 7.2),
    legend.title = element_text(size = 7.2),
    legend.text  = element_text(size = 6.6),
    plot.margin = margin(2, 2, 2, 2)
  )

p_fmle_ctp <- ggplot(
  sub_ctp,
  aes(Pearson_ctp, Pearson_fmle, color = protein, shape = cell_type)
) +
  geom_point(size = 2.1, alpha = 0.9, stroke = 0.25) +
  geom_abline(slope = 1, intercept = 0, color = "grey45", linewidth = 0.35) +
  coord_cartesian(xlim = xlim2, ylim = ylim2, expand = FALSE) +
  scale_shape_manual(values = celltype_shapes, drop = FALSE) +
  guides(
    shape = guide_legend(override.aes = list(size = 2.8, alpha = 1)),
    color = guide_legend(override.aes = list(size = 3.0, alpha = 1))
  ) +
  theme_classic(base_size = 12) +
  labs(x = "cTPnet", y = "FMLE", color = "Protein", shape = "Cell type") +
  theme(
    axis.title = element_text(size = 8.2),
    axis.text  = element_text(size = 7.2),
    legend.title = element_text(size = 7.2),
    legend.text  = element_text(size = 6.6),
    plot.margin = margin(2, 2, 2, 2)
  )

plot_panelD_method <- function(prot,
                               method = c("FMLE","scLinear","cTPnet"),
                               fig2,
                               celltype_vec,
                               n_points = 6000,
                               add_lm = FALSE) {
  
  method <- match.arg(method)
  
  pred_list <- switch(
    method,
    "FMLE"     = fig2$pred_fmle,
    "scLinear" = fig2$pred_sclinear,
    "cTPnet"   = fig2$pred_ctpnet
  )
  
  stopifnot(prot %in% names(fig2$truth),
            prot %in% names(pred_list),
            prot %in% names(fig2$cells))
  
  cells <- fig2$cells[[prot]]
  
  df <- data.frame(
    pred = as.numeric(pred_list[[prot]]),
    real = as.numeric(fig2$truth[[prot]]),
    cell = cells,
    stringsAsFactors = FALSE
  )
  
  df <- df[is.finite(df$pred) & is.finite(df$real), , drop = FALSE]
  
  df$celltype <- unname(celltype_vec[df$cell])
  df$celltype[is.na(df$celltype)] <- "Unassigned"
  df <- df[!df$celltype %in% c("Unassigned", "other", "other T"), , drop = FALSE]
  df$celltype <- droplevels(factor(df$celltype))
  
  if (nrow(df) > n_points) {
    set.seed(1)
    df <- df[sample.int(nrow(df), n_points), , drop = FALSE]
  }
  
  r   <- suppressWarnings(cor(df$pred, df$real, method = "pearson"))
  rho <- suppressWarnings(cor(df$pred, df$real, method = "spearman"))
  
  lim <- range(c(df$pred, df$real), na.rm = TRUE)
  
  prot_m <- prot
  prot_m <- gsub("^anti-human-", "", prot_m)
  prot_m <- gsub("^anti-mouse-human-", "", prot_m)
  prot_m <- gsub("-totalC$", "", prot_m)
  
  prot_m <- dplyr::recode(
    prot_m,
    "CD127-IL7Ra" = "CD127",
    "CD197-CCR7"  = "CCR7"
  )
  
  p <- ggplot(df, aes(x = pred, y = real, color = celltype)) +
    geom_point(alpha = 0.75, size = 1.2) +
    geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed",
      linewidth = 0.7,
      color = "grey25"
    ) +
    coord_equal(xlim = lim, ylim = lim) +
    theme_classic(base_size = 12) +
    labs(
      title = sprintf("%s  (r=%.2f, ρ=%.2f)", prot_m, r, rho),
      x = sprintf("%s predicted", method),
      y = "True protein",
      color = NULL
    ) +
    theme(
      plot.title = element_text(
        size = 10, hjust = 0, face = "bold",
        margin = margin(t = 12, b = -1)
      ),
      legend.position = "right"
    )
  
  if (add_lm) {
    p <- p + geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, color = "black")
  }
  
  p
}


make_donor_header <- function(tag, donor) {
  wrap_elements(full = grobTree(
    textGrob(tag,
             x = 0.02, y = 0.5,
             just = c("left", "center"),
             gp = gpar(fontface = "bold", fontsize = 12)),
    textGrob(donor,
             x = 0.50, y = 0.5,
             just = c("center", "center"),
             gp = gpar(fontface = "bold", fontsize = 12))
  ))
}

p1 <- plot_panelD_method("anti-human-CD58-totalC", "FMLE", fig21, celltype_vec1)
p2 <- plot_panelD_method("anti-human-CD58-totalC", "scLinear", fig21, celltype_vec1)
p3 <- plot_panelD_method("anti-human-CD58-totalC", "cTPnet", fig21, celltype_vec1)
p_CD58NS1 <- make_donor_header("(d)", "NS1") / p1 / p2 / p3 +
  plot_layout(heights = c(-0.10, 1, 1, 1))

p1 <- plot_panelD_method("anti-human-CD58-totalC", "FMLE", fig22, celltype_vec2)
p2 <- plot_panelD_method("anti-human-CD58-totalC", "scLinear", fig22, celltype_vec2)
p3 <- plot_panelD_method("anti-human-CD58-totalC", "cTPnet", fig22, celltype_vec2)
p_CD58TP7 <- make_donor_header("(e)", "TP7") / p1 / p2 / p3 +
  plot_layout(heights = c(-0.10, 1, 1, 1))

p1 <- plot_panelD_method("anti-human-CD58-totalC", "FMLE", fig23, celltype_vec3)
p2 <- plot_panelD_method("anti-human-CD58-totalC", "scLinear", fig23, celltype_vec3)
p3 <- plot_panelD_method("anti-human-CD58-totalC", "cTPnet", fig23, celltype_vec3)
p_CD58TS5 <- make_donor_header("(f)", "TS5") / p1 / p2 / p3 +
  plot_layout(heights = c(-0.10, 1, 1, 1))

rm(p1, p2, p3)


# Win percentage across 3 datasets

wins <- df_all %>%
  filter(is.finite(Pearson)) %>%
  group_by(dataset, protein) %>%
  mutate(best = max(Pearson, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(is_win = (Pearson == best)) %>%
  group_by(dataset, method) %>%
  summarise(
    n_prot = n_distinct(protein),
    win_pct = mean(is_win) * 100,
    .groups = "drop"
  )

p_win <- ggplot(wins, aes(x = method, y = win_pct, fill = method)) +
  geom_col(width = 0.75) +
  facet_wrap(~ dataset, nrow = 1) +
  labs(x = NULL, y = "Win % (best per protein)") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


p_win0 <- p_win + theme(legend.position = "none") + guides(fill = "none", color = "none")
p_fmle_scl <- p_fmle_scl +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))
p_fmle_ctp <- p_fmle_ctp +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))


pA <- pA +
  labs(tag="(a)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position = c(0.02,1.0))

pB <- (p_fmle_scl / p_fmle_ctp) +
  labs(tag="(b)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position = c(0.02, 2.080))

p_win0 <- p_win0 +
  labs(tag="(c)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position = c(0.02, 1.01)) +theme(
          strip.text = element_text(face = "bold")
        )

final <- ( pA
           | (pB / p_win0)
           | (p_CD58NS1 | p_CD58TP7 | p_CD58TS5)
) +
  plot_layout(
    guides = "collect",
    widths = c(0.7, 0.7, 1.70)
  ) &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 12)
  )


ggsave(file.path(base1, "Figure_5.pdf"),
       plot = final,
       device = cairo_pdf,
       width = 18, height = 10, units = "in",
       dpi = 300,
       limitsize = FALSE)


