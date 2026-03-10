suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(grid)
})

bench <- 1  # set 1 or 2 or 3
ds      <- "citeseq_v1"
base <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_%d", bench))
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ctp    <- file.path(base, "ctp") 
scl    <- file.path(base, "sclinear") 
fig_dir <- file.path(base, "paper_figures", ds)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
fig2_out <- file.path(fig_dir, "fig3_examples_same.rds")
fig2 <- readRDS(fig2_out)
available <- names(fig2$truth)
print(available)
seu_prep <- readRDS(file.path(base, ds, "seu_final.rds"))
stopifnot("cell_type" %in% colnames(seu_prep@meta.data))
celltype_vec <- setNames(as.character(seu_prep@meta.data$cell_type), colnames(seu_prep))

fmle_ds_by_bench <- c(
  "citeseq_v1_final",  # bench 1 FMLE folder
  "citeseq_v1_final",  # bench 2 FMLE folder (edit if different)
  "citeseq_v1_final"   # bench 3 FMLE folder (edit if different)
)


bench_name <- c(
  `1` = "Kaggle PBMC",
  `2` = "10x PBMC 10k",
  `3` = "TEA-seq PBMC"
)
read_one_bench <- function(bench, fmle_ds){
  
  base <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_%d", bench))
  
  fmle_path <- file.path(base, "FMLE", fmle_ds, "fmle_test_metrics.csv")
  scl_path  <- file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv")
  ctp_path  <- file.path(base, "ctp", "ctpnet_test_metrics_FMLEscale.csv")
  
  stopifnot(file.exists(fmle_path),
            file.exists(scl_path),
            file.exists(ctp_path))
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

# -------- BUILD DATA --------
df_all <- bind_rows(
  read_one_bench(1, fmle_ds_by_bench[1]),
  read_one_bench(2, fmle_ds_by_bench[2]),
  read_one_bench(3, fmle_ds_by_bench[3])
) %>% 
  filter(is.finite(Pearson))%>%
  mutate(
    dataset = factor(dataset, levels = c("Kaggle PBMC", "10x PBMC 10k", "TEA-seq PBMC")),
    method  = factor(method,  levels = c("scLinear","cTPnet","FMLE"))
  )

# -------- PANEL A PLOT --------
pA <- ggplot(df_all,
             aes(x=method, y=Pearson, color=method)) +
  
  geom_boxplot(outlier.shape=NA, width=0.6) +
  
  geom_jitter(width=0.15, alpha=0.7, size=2) +
  
  facet_wrap(~dataset, ncol=1) +
  
  theme_classic(base_size=15) +
  
  labs(
    x=NULL,
    y="Per-protein Pearson correlation"
  ) +
  
  theme(
    legend.position="none",
    strip.text=element_text(face="bold")
  )



# Panel B — FMLE vs scLinear per protien

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


# panel C D and E

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
  
  df <- df[is.finite(df$pred) & is.finite(df$real), , drop=FALSE]
  
  # attach cell types
  df$celltype <- unname(celltype_vec[df$cell])
  df$celltype[is.na(df$celltype)] <- "Unassigned"
  df <- df[df$celltype != "Unassigned", , drop=FALSE]
  df$celltype <- droplevels(factor(df$celltype))
  
  # subsample for visual clarity
  if (nrow(df) > n_points) {
    set.seed(1)
    df <- df[sample.int(nrow(df), n_points), , drop=FALSE]
  }
  
  r   <- suppressWarnings(cor(df$pred, df$real, method="pearson"))
  rho <- suppressWarnings(cor(df$pred, df$real, method="spearman"))
  
  lim <- range(c(df$pred, df$real), na.rm=TRUE)
  # lim <- range(fig2$truth[[prot]], na.rm=TRUE)
  prot_m <- dplyr::recode(prot,
                          "CD127-IL7Ra" = "CD127",
                          "CD197-CCR7"  = "CCR7",
                          .default      = prot
  )
  p <- ggplot(df, aes(x = pred, y = real, color = celltype)) +
    geom_point(alpha = 0.75, size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = 0.7, color = "grey25") +
    coord_equal(xlim = lim, ylim = lim) +
    theme_classic(base_size = 15) +
    labs(
      title = sprintf("%s  (r=%.2f, ρ=%.2f)", prot_m,  r, rho),
      x = sprintf("%s predicted", method),
      y = "True protein",
      color = NULL
    ) +
    theme(plot.title = element_text(size = 10, hjust = 0, face = "bold", margin = margin(t = 12, b = -1)),
          legend.position = "right")
  
  if (add_lm) {
    p <- p + geom_smooth(method="lm", se=FALSE, linewidth=0.6, color="black")
  }
  
  p
}


p1 <- plot_panelD_method("CD197-CCR7", "FMLE", fig2, celltype_vec)
p2 <- plot_panelD_method("CD197-CCR7", "scLinear", fig2, celltype_vec)
p3 <- plot_panelD_method("CD197-CCR7", "cTPnet", fig2, celltype_vec)
p_CD197 <- (p1 / p2 / p3)

p1 <- plot_panelD_method("CD14", "FMLE", fig2, celltype_vec)
p2 <- plot_panelD_method("CD14", "scLinear", fig2, celltype_vec)
p3 <- plot_panelD_method("CD14", "cTPnet", fig2, celltype_vec)
p_CD14 <- (p1 / p2 / p3)

p1 <- plot_panelD_method("CD127-IL7Ra", "FMLE", fig2, celltype_vec)
p2 <- plot_panelD_method("CD127-IL7Ra", "scLinear", fig2, celltype_vec)
p3 <- plot_panelD_method("CD127-IL7Ra", "cTPnet", fig2, celltype_vec)
p_CD127 <- (p1 / p2 / p3)

rm(p1, p2,p3)

# final <- (pA | (p_fmle_scl / p_fmle_ctp) | p_CD14 | p_CD127 | p_CD197) +
#   plot_layout(guides = "collect") &
#   theme(
#     legend.position = "right",
#     legend.box = "vertical",
#     legend.key.size = unit(0.5, "lines"),
#     legend.text = element_text(size = 8)
#   )
# 
# final

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
  theme_classic(base_size = 14) +
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
        plot.tag.position=c(0.5,1.0))

pB <- (p_fmle_scl / p_fmle_ctp) +
  labs(tag="(b)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position=c(0.5,2.20))

p_win0 <- p_win0 +
  labs(tag="(c)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position=c(0.2,1.05))

p_CD14 <- p_CD14 +
  labs(tag="(d)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position=c(0.5,3.010))

p_CD127 <- p_CD127 +
  labs(tag="(e)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position=c(0.5,3.010))

p_CD197 <- p_CD197 +
  labs(tag="(f)") +
  theme(plot.tag=element_text(face="bold", hjust=0.5),
        plot.tag.position=c(0.5,3.010))



final <- ( pA
           | (pB/ p_win0)
           | (p_CD14 | p_CD127 | p_CD197)
) +
  plot_layout(
    guides = "collect",
    widths = c(0.7, 0.7, 1.70)   # tweak if needed
  ) &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.key.size = unit(0.5, "lines"),
    legend.text = element_text(size = 12)
  )

final

ggsave(file.path(out_dir, "Fig2_FMLE_main.pdf"),
       plot = final,
       device = cairo_pdf,
       width = 21, height = 10, units = "in",
       dpi = 300,
       limitsize = FALSE)















