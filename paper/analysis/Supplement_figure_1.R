
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(grid)
  library(hexbin)
  library(forcats)
  library(scales)
})

#--------------------------------------------------------------------#
# Supplement figures
#--------------------------------------------------------------------#

safe_r2 <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]; yhat <- yhat[ok]
  if (length(y) < 3) return(NA_real_)
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}


# ----------------------------
#  Supplement
# 1) CD14 (myeloid/APC) 2) CD3 (pan-T) 3) CD127-IL7Ra (regulatory/noisy)
# 4) CD45RA or CD45RO (state biology) 5) CD56 (NK lineage) 6) CD19 (B cells) (optional but strong)

# prot1 <- "CD45RA"
# prot2 <- "CD45RO"
# rna_gene <- "PTPRC"
# 
# adt_raw_ra <- as.numeric(adt_mat[prot1, cells]); names(adt_raw_ra) <- cells
# adt_raw_ro <- as.numeric(adt_mat[prot2, cells]); names(adt_raw_ro) <- cells
# rna_raw    <- as.numeric(rna_mat[rna_gene, cells]); names(rna_raw) <- cells
# 
# c(
#   cor_rna_RA = cor(rna_raw, adt_raw_ra, method="spearman"),
#   cor_rna_RO = cor(rna_raw, adt_raw_ro, method="spearman"),
#   cor_RA_RO  = cor(adt_raw_ra, adt_raw_ro, method="spearman")
# )

# get_pair <- function(prot, rna_gene){
#   adt <- as.numeric(adt_mat[prot, cells]); names(adt) <- cells
#   rna <- as.numeric(rna_mat[rna_gene, cells]); names(rna) <- cells
#   c(spearman = cor(rna, adt, method="spearman"),
#     pearson  = cor(rna, adt, method="pearson"))
# }
# 
# get_pair("CD19", "CD19")
# get_pair("CD56", "NCAM1")
# ----------------------------
# prot <- "CD3"
# rna_gene <- "CD3E"
# prot <- "CD14"
# rna_gene <- "CD14"

prot     <- "CD127-IL7Ra"
rna_gene <- "IL7R"

# prot <- "CD56"  
# rna_gene <- "NCAM1"

# prot <- "CD19"  
# rna_gene <- "CD19"

bench <- 1  # set 1, 2, or 3
ds <- "citeseq_v1"
base <- file.path(cfg$out_root, sprintf("benchmarks_%d", bench))
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seu_prep <- readRDS(file.path(base, ds, "seu_final.rds")) 
stopifnot(file.exists(file.path(base, ds, "seu_final.rds")))
stopifnot(prot %in% rownames(seu_prep[["ADT"]]))
stopifnot(rna_gene %in% rownames(seu_prep[["RNA"]]))

obj   <- readRDS(file.path(out_dir, paste0("final_", prot, ".rds")))
alpha <- obj$alpha_test
stopifnot(!is.null(alpha), !is.null(rownames(alpha)))

cells <- rownames(alpha)          # definitive cell IDs
hard  <- max.col(alpha)
hard  <- factor(hard)

# regime counts for legend
tab_h <- table(hard)
print(tab_h)
# ------------------------------------------------------------------------------------
#  In case of CD3 
adt_mat <- GetAssayData(seu_prep, assay = "ADT", layer = "data")
stopifnot(prot %in% rownames(adt_mat))

adt_raw <- as.numeric(adt_mat[prot, cells])
names(adt_raw) <- cells
rna_mat <- GetAssayData(seu_prep, assay = "RNA", layer = "data")
stopifnot(rna_gene %in% rownames(rna_mat))

rna_raw <- as.numeric(rna_mat[rna_gene, cells])
names(rna_raw) <- cells
cells_r6 <- cells[hard == 6]

summary(rna_raw[cells_r6])
summary(adt_raw[cells_r6])
qc <- FetchData(seu_prep, vars=c("nCount_RNA","nFeature_RNA","percent.mt"), cells=cells)
summary(qc[cells_r6, ])

# compare to others
cells_not6 <- setdiff(cells, cells_r6)
summary(qc[cells_not6, ])

cd3d <- rna_mat["CD3D", cells_r6]
cd3g <- rna_mat["CD3G", cells_r6]

summary(cd3d)
summary(cd3g)

cd19 <- rna_mat["CD19", cells_r6]
cd19 <- rna_mat["CD19", cells_r6]
summary(cd19)
summary(cd19)

IL7R <- rna_mat["IL7R", cells_r6]
IL7R <- rna_mat["IL7R", cells_r6]

summary(IL7R)
summary(IL7R)

# ------------------------------------------------------------------------------------


hard_lab <- factor(hard, levels=levels(hard),
                   labels=paste0(levels(hard), " (n=", as.integer(tab_h[levels(hard)]), ")"))

# Align y/yhat to alpha rows (THIS is the key)
stopifnot(length(obj$y_test) == length(cells),
          length(obj$yhat_test) == length(cells))

rna_raw <- FetchData(seu_prep, vars=rna_gene, cells=cells)[,1]
rna_log <- log1p(as.numeric(rna_raw))

df <- data.frame(
  cell = cells,
  rna  = rna_log,
  y    = as.numeric(obj$y_test),
  yhat = as.numeric(obj$yhat_test),
  hard = hard_lab
)

df <- df[is.finite(df$rna) & is.finite(df$y) & is.finite(df$yhat), ]

df$rna_raw <- as.numeric(rna_raw)
df$hard <- factor(df$hard)

fit_int <- lm(y ~ rna * hard, data = df)
p_int <- anova(fit_int)[["Pr(>F)"]][match("rna:hard", rownames(anova(fit_int)))]
p_txt <- paste0("Interaction p=", format.pval(p_int, digits = 2))

pA <- ggplot(df, aes(rna, y, color = hard)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  theme_classic(base_size = 12) +
  labs(x = paste0(rna_gene, " RNA (log1p)"),
       y = paste0(prot, " protein"),
       color = "regime",
       title = "Regime-dependent RNA–protein coupling") +
  annotate("text", x = Inf, y = Inf, label = p_txt,
           hjust = 1.9, vjust = 1.2, size = 3.6) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# -----------------------------
# PANEL C — global failure (diagnostic, same cells)
# -----------------------------
m_global <- lm(y ~ rna, data=df)
df$resid_global <- df$y - predict(m_global, newdata=df)
# df_plot <- df[df$rna_raw > 0, ] 

pB <- ggplot(df, aes(rna, resid_global, color=hard)) +
  geom_point(alpha=0.45, size=0.7) +
  geom_hline(yintercept=0, linetype=2, linewidth=0.7) +
  theme_classic(base_size=12) +
  labs(x=paste0(rna_gene, " RNA (log1p)"),
       y="Global model residual",
       title="Residuals from a global mapping") +
  guides(color="none") +
  theme(plot.title = element_text(hjust = 0.9, face="bold"))
# labs(x=paste0(rna_gene, " RNA (log1p)"),
#      y="Global model residual",
#      color="FMLE regime",
#      title="Residuals under a global mapping")

# -----------------------------
# PANEL D — FMLE vs truth + (optional) global overlay
# -----------------------------
# global prediction from same model (optional overlay)
df$yhat_global <- predict(m_global, newdata=df)

df_plotD <- bind_rows(
  df %>% transmute(y, pred = yhat_global, model = "Global"),
  df %>% transmute(y, pred = yhat,        model = "FMLE")
) %>% filter(is.finite(y), is.finite(pred))

lims <- range(c(df_plotD$y, df_plotD$pred), finite = TRUE)
pad  <- 0.02 * diff(lims)
lims <- lims + c(-pad, pad)

metricsD <- df_plotD %>%
  group_by(model) %>%
  summarise(
    r2 = safe_r2(y, pred),
    r  = suppressWarnings(cor(y, pred, method="pearson", use="complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(lbl = sprintf("R²=%.3f; r=%.3f", r2, r),
         x = lims[2], y = lims[1])

cols <- c(Global = "#0072B2", FMLE = "#D55E00")

base_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none")

make_panelD <- function(m){
  ggplot(df_plotD %>% filter(model == m), aes(x = y, y = pred)) +
    geom_hex(bins = 55) +
    scale_fill_gradient(low = "grey80", high = cols[[m]], trans = "sqrt", guide = "none") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8, color="grey35") +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE)+
    labs(title = m, x = "True protein (cells with α)", y = "Prediction") +
    geom_text(
      data = metricsD %>% filter(model == m),
      aes(x = x, y = y, label = lbl),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = -3.0, size = 3.4
    ) +
    base_theme
}

# pD <- make_panelD("Global") + make_panelD("FMLE") +
#   plot_annotation(title = "Prediction vs truth (same cells)") &
#   theme(plot.title = element_text(hjust = 0.5))

pC <- make_panelD("Global")+
  plot_annotation(title = "Prediction vs truth (same cells)") &
  theme(plot.title = element_text(hjust = 0.5))

pD   <- make_panelD("FMLE")+
  plot_annotation(title = "Prediction vs truth (same cells)") &
  theme(plot.title = element_text(hjust = 0.5))

# Keep only cells that exist in Seurat AND in UMAP embeddings
cells <- rownames(alpha)
cells <- intersect(cells, colnames(seu_prep))

umap_all <- Embeddings(seu_prep, "umap")
cells <- intersect(cells, rownames(umap_all))

if (length(cells) < 200) stop("Too few overlapping cells between alpha_test and Seurat UMAP.")

# Reindex alpha + umap to SAME order
alpha <- alpha[cells, , drop = FALSE]
umap  <- umap_all[cells, , drop = FALSE]

# Hard regime + entropy
hard <- factor(max.col(alpha), levels = sort(unique(max.col(alpha))))
entropy <- -rowSums(alpha * log(alpha + 1e-12))
entropy <- entropy / log(ncol(alpha))
# ----------------------------
# PANEL A — UMAP colored by hard regime
# ----------------------------
dfA <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2], regime = hard)

pE <- ggplot(dfA, aes(UMAP1, UMAP2, color = regime)) +
  geom_point(size = 0.35, alpha = 0.85) +
  theme_classic(base_size = 12) +
  labs(title = "FMLE regimes (argmax α)", color = "regime") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# # ----------------------------
# # PANEL B — UMAP colored by entropy
# # (NO viridis dependency; keep ggplot default scale unless you insist)
# # ----------------------------
# dfB <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2], entropy = entropy)
# 
# pG <- ggplot(dfB, aes(UMAP1, UMAP2, color = entropy)) +
#   geom_point(size = 0.35, alpha = 0.85) +
#   theme_classic(base_size = 12) +
#   labs(title = "Regime uncertainty (α-entropy)", color = "Entropy") +
#   theme(plot.title = element_text(hjust = 0.5, face="bold"))

# cell composition


celltype_col <- "cell_type" 
entropy_vec <- function(p) {
  p <- p / sum(p)
  p <- p[p > 0]
  -sum(p * log(p))
}

# Normalize entropy to [0,1] by dividing by log(R)
entropy_norm <- function(p) {
  R <- length(p)
  if (R <= 1) return(0)
  entropy_vec(p) / log(R)
}


seu <- seu_prep
alpha <- as.matrix(alpha)
if (is.null(rownames(alpha))) {
  stop("alpha must have rownames = cell IDs.")
}
if (!all(rownames(alpha) %in% colnames(seu))) {
  stop("Some alpha cells not in Seurat object.")
}

# Subset and order to match alpha rows
seu2 <- seu[, rownames(alpha)]

if (!identical(colnames(seu2), rownames(alpha))) {
  stop("Cell ordering mismatch after subsetting.")
}



rs <- rowSums(alpha)
if (any(!is.finite(rs)) || any(rs <= 0)) stop("alpha has non-finite or non-positive row sums.")
alpha <- alpha / rs

R <- ncol(alpha)
regime_names <- colnames(alpha)
if (is.null(regime_names)) {
  regime_names <- paste0("Regime_", seq_len(R))
  colnames(alpha) <- regime_names
}

meta <- seu2@meta.data
if (!(celltype_col %in% colnames(meta))) {
  stop(paste0("Missing metadata column: ", celltype_col))
}

cell_ids <- colnames(seu2)
celltype <- meta[cell_ids, celltype_col, drop=TRUE] |> as.character()
celltype[is.na(celltype) | celltype == ""] <- "Unknown"

regime_id <- max.col(alpha, ties.method = "first")
regime <- factor(colnames(alpha)[regime_id], levels = colnames(alpha))

ent <- apply(alpha, 1, entropy_norm)

df_cell <- tibble(
  cell = cell_ids,
  celltype = factor(celltype),
  regime = regime,
  entropy = as.numeric(ent)
)

min_n <- 100
ct_counts <- table(df_cell$celltype)
keep_ct <- names(ct_counts)[ct_counts >= min_n]

df_cell <- df_cell %>%
  filter(as.character(celltype) %in% keep_ct) %>%
  mutate(celltype = droplevels(factor(celltype)))

# Order cell types by abundance (clean)
df_cell <- df_cell %>%
  mutate(celltype = fct_infreq(celltype))

df_comp <- df_cell %>%
  count(celltype, regime, name = "n") %>%
  group_by(celltype) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup()

tab <- df_comp %>%
  select(celltype, regime, n) %>%
  tidyr::pivot_wider(names_from=regime, values_from=n, values_fill=0)

chisq.test(as.matrix(tab[,-1]))
library(rcompanion)
cramerV(as.matrix(tab[,-1]))

df_comp$celltype <- factor(df_comp$celltype,
                           levels = c("Monocyte", "B", "NK cells", "T"))

pF <- ggplot(df_comp, aes(x = celltype, y = frac, fill = regime)) +
  geom_col(width = 0.85, color = NA) +
  scale_x_discrete(labels = c("Monocyte"="Monocyte", "B"="B", "NK cells"="NK", "T"="T"))  +
  labs(x = NULL, y = "Fraction of cells",
       title = "Regime within cell types") +
  guides(color="none") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold"),
    legend.position = "none" 
  )

# ----------------------------
# Assemble + save
# ----------------------------

base_theme <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.margin = margin(14, 6, 6, 6)
  )
tag_theme <- theme(plot.tag = element_text(face = "bold"),
                   plot.tag.position = c(0.5, 1.0))  

pA_std <- pA  + labs(tag="(a)") + tag_theme
pB_std <- pB + base_theme + labs(tag="(b)") + tag_theme
pC_std <- pC + base_theme + labs(tag="(c)") + tag_theme

pD_std  <- pD + base_theme + labs(tag="(d)") + tag_theme
pE_std  <- pE + base_theme + labs(tag="(e)") + tag_theme
pF_std  <- pF + base_theme + labs(tag="(f)") + tag_theme
# pF1_std <- p_b1 + base_theme + labs(tag="(g)") + tag_theme + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# pG_std  <- pG + base_theme + labs(tag="(h)") + tag_theme
# pH_std  <- pH + base_theme + labs(tag="(i)") + tag_theme
# pI_std  <- pI + base_theme + labs(tag="(j)") + tag_theme




row1 <- (pA_std | pB_std | pC_std) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

row2 <- row2 <- (pD_std | pE_std | pF_std )+
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


tag_title_fix <- theme(
  plot.title.position = "plot",
  plot.title = element_text(margin = margin(t = 6, b = 8)),  # adds space below title
  plot.tag = element_text(margin = margin(b = 10, r = 6))   # moves tag away from title                              # top-left
)

figure1_1 <- (row1 / row2) & tag_title_fix
figure1_1


ggsave(file.path(out_dir, "Fig1_FMLE_supp_CD19.pdf"),
       plot = figure1_1,
       device = cairo_pdf,
       width = 14,
       height = 12,
       units = "in")

