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
# Figure 1 Main
#--------------------------------------------------------------------#


safe_r2 <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]; yhat <- yhat[ok]
  if (length(y) < 3) return(NA_real_)
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

# -----------------------------
# PANEL A — schematic
# -----------------------------
pA <- ggplot() +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
  
  # Global row
  annotate("label", x=0.10, y=0.82, label="RNA", label.size=0.3, size=4) +
  annotate("label", x=0.88, y=0.82, label="Protein", label.size=0.3, size=4) +
  annotate("text",  x=0.50, y=0.92, label="Global: single mapping\nP = f(X)", size=4,hjust = 0.5) +
  geom_segment(aes(x=0.20, y=0.82, xend=0.75, yend=0.82),
               arrow=arrow(length=unit(0.02,"npc")), linewidth=0.8) +
  
  # FMLE row
  annotate("text", x=0.50, y=0.70, label="FMLE: mixture of coupling regimes", size=4) +
  annotate("label", x=0.10, y=0.48, label="RNA", label.size=0.3, size=4) +
  annotate("label", x=0.88, y=0.48, label="Protein", label.size=0.3, size=4) +
  
  geom_segment(aes(x=0.18, y=0.48, xend=0.33, yend=0.48),
               arrow=arrow(length=unit(0.02,"npc")), linewidth=0.8) +
  annotate("label", x=0.40, y=0.48, label=expression(alpha(X)),
           label.size=0.3, size=3.6) +
  
  geom_segment(aes(x=0.47, y=0.48, xend=0.60, yend=0.56),
               arrow=arrow(length=unit(0.02,"npc")), linewidth=0.8) +
  geom_segment(aes(x=0.47, y=0.48, xend=0.60, yend=0.40),
               arrow=arrow(length=unit(0.02,"npc")), linewidth=0.8) +
  annotate("text", x=0.50, y=0.60, label="regime 1:  f1(X)", hjust=0, size=3.6) +
  annotate("text", x=0.50, y=0.36, label="regime 2:  f2(X)", hjust=0, size=3.6) +
  
  geom_segment(aes(x=0.62, y=0.56, xend=0.77, yend=0.48),
               arrow=arrow(length=unit(0.02,"npc")), linewidth=0.8) +
  geom_segment(aes(x=0.62, y=0.40, xend=0.77, yend=0.48),
               arrow=arrow(length=unit(0.02,"npc")), linewidth=0.8) +
  
  # Equation row (small, correct operator)
  annotate("text", x=0.50, y=0.20,
           label=expression(P == alpha[1](X)*f[1](X) + alpha[2](X)*f[2](X)),
           size=3.4) 

pA <- pA +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    line       = element_blank()
  )

# ----------------------------
# Load FMLE object (TEST-only alpha)
# ----------------------------
prot <- "HLA.DR"
rna_gene <- "HLA-DRA"
bench <- 1  # set 1 or 2 or 3
base  <- path.expand(sprintf("~/Desktop/FMLE/benchmarks_%d", bench))

ds    <- "citeseq_v1"
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
seu_prep <- readRDS(file.path(base, ds, "seu_final.rds")) 

obj   <- readRDS(file.path(out_dir, paste0("final_", prot, ".rds")))
alpha <- obj$alpha_test
stopifnot(!is.null(alpha), !is.null(rownames(alpha)))

cells <- rownames(alpha)          # definitive cell IDs
hard  <- max.col(alpha)
hard  <- factor(hard)

# regime counts for legend
tab_h <- table(hard)
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

# -----------------------------
# PANEL B — slopes + interaction p
# -----------------------------
# df$rna_raw <- as.numeric(rna_raw)
# df_plot <- df[df$rna_raw > 0, ] 
# fit_int <- lm(y ~ rna * hard, data=df)
# p_int <- anova(fit_int)[["Pr(>F)"]][3]
# p_txt <- paste0("Interaction p=", format.pval(p_int, digits=2))
# pB <- ggplot(df_plot, aes(rna, y, color=hard)) +
#   geom_point(alpha=0.45, size=0.7) +
#   geom_smooth(method="lm", se=TRUE, linewidth=0.9) +
#   theme_classic(base_size=12) +
#   labs(x=paste0(rna_gene, " RNA (log1p)"),
#        y=paste0(prot, " protein"),
#        color="FMLE regime",
#        title="Regime-dependent RNA–protein coupling") +
#   annotate("text", x=Inf, y=Inf, label=p_txt, hjust=1.7, vjust=1.0, size=3.6) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face="bold")
#   )
df$rna_raw <- as.numeric(rna_raw)
df$hard <- factor(df$hard)

fit_int <- lm(y ~ rna * hard, data = df)
p_int <- anova(fit_int)[["Pr(>F)"]][match("rna:hard", rownames(anova(fit_int)))]
p_txt <- paste0("Interaction p=", format.pval(p_int, digits = 2))

pB <- ggplot(df, aes(rna, y, color = hard)) +
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

pC <- ggplot(df, aes(rna, resid_global, color=hard)) +
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

pD <- make_panelD("Global")+
  plot_annotation(title = "Prediction vs truth (same cells)") &
  theme(plot.title = element_text(hjust = 0.5))

pE   <- make_panelD("FMLE")+
  plot_annotation(title = "Prediction vs truth (same cells)") &
  theme(plot.title = element_text(hjust = 0.5))

# ----------------------------
prot <- "HLA.DR" # HLA.DR or CD14
obj <- readRDS(file.path(out_dir, paste0("final_", prot, ".rds")))
alpha <- obj$alpha_test
stopifnot(!is.null(alpha), !is.null(rownames(alpha)))

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

pF <- ggplot(dfA, aes(UMAP1, UMAP2, color = regime)) +
  geom_point(size = 0.35, alpha = 0.85) +
  theme_classic(base_size = 12) +
  labs(title = "FMLE regimes (argmax α)", color = "regime") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# ----------------------------
# PANEL B — UMAP colored by entropy
# (NO viridis dependency; keep ggplot default scale unless you insist)
# ----------------------------
dfB <- data.frame(UMAP1 = umap[,1], UMAP2 = umap[,2], entropy = entropy)

pG <- ggplot(dfB, aes(UMAP1, UMAP2, color = entropy)) +
  geom_point(size = 0.35, alpha = 0.85) +
  theme_classic(base_size = 12) +
  labs(title = "Regime uncertainty (α-entropy)", color = "Entropy") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# If you really want viridis, do:
# + scale_color_viridis_c()

# ----------------------------
# PANEL C — Signature scores by regime (robust gene filtering)
# ----------------------------
sig_apc <- c("HLA-DRA","HLA-DRB1","CD74","CIITA")
sig_inf <- c("S100A8","S100A9","LYZ","FCN1","LGALS3")

seu_sub <- subset(seu_prep, cells = cells)

# Filter to genes present; avoid AddModuleScore crashing / nonsense
sig_apc2 <- intersect(sig_apc, rownames(seu_sub))
sig_inf2 <- intersect(sig_inf, rownames(seu_sub))

if (length(sig_apc2) < 3) warning("APC signature has <3 genes present; interpretation weak.")
if (length(sig_inf2) < 3) warning("Inflam signature has <3 genes present; interpretation weak.")

if (length(sig_apc2) >= 3) seu_sub <- AddModuleScore(seu_sub, features = list(sig_apc2), name="APC", search=FALSE)
if (length(sig_inf2) >= 3) seu_sub <- AddModuleScore(seu_sub, features = list(sig_inf2), name="Inflam", search=FALSE)

md <- seu_sub@meta.data
md$cell <- rownames(md)
md$regime <- hard[match(md$cell, cells)]

score_cols <- intersect(c("APC1","Inflam1"), colnames(md))
if (length(score_cols) == 0) stop("No module scores computed (gene lists missing).")

dfC <- md %>%
  dplyr::select(regime, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols), names_to = "signature", values_to = "score")

pH <- ggplot(dfC, aes(regime, score, fill = regime)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.4) +
  facet_wrap(~ signature, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(title = "Biological signatures by regime",
       x = "FMLE regime", y = "Module score") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face="bold"))

# ----------------------------
# PANEL D — Regime gene programs (β weights) — FIXED ordering per facet
# ----------------------------
beta <- obj$beta
stopifnot(!is.null(beta), !is.null(rownames(beta)))
stopifnot(ncol(beta) >= 2)

top_beta <- function(b, k = 10){
  b <- b[is.finite(b)]
  pos <- head(sort(b, decreasing = TRUE), k)
  neg <- head(sort(b, decreasing = FALSE), k)
  c(pos, neg)
}

dfD <- bind_rows(lapply(seq_len(ncol(beta)), function(r){
  b <- beta[, r]; names(b) <- rownames(beta)
  tb <- top_beta(b, k = 5)
  tibble(
    gene = names(tb),
    beta = as.numeric(tb),
    regime = paste0("Regime ", r),
    sign = ifelse(tb >= 0, "pos", "neg")
  )
}))

# Critical: reorder genes WITHIN each (regime, sign) facet
dfD <- dfD %>%
  group_by(regime, sign) %>%
  mutate(gene = fct_reorder(gene, beta)) %>%
  ungroup()

pI <- ggplot(dfD, aes(gene, beta)) +
  geom_col() +
  coord_flip() +
  facet_grid(sign ~ regime, scales = "free_y", space = "free_y") +
  theme_classic(base_size = 7) +
  labs(title = "Regime gene programs (β weights)", x = NULL, y = "β") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))


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

df_comp$celltype <- factor(df_comp$celltype,
                           levels = c("Monocyte", "B", "NK cells", "T"))

p_b1 <- ggplot(df_comp, aes(x = celltype, y = frac, fill = regime)) +
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
pF1_std <- p_b1 + base_theme + labs(tag="(g)") + tag_theme + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
pG_std  <- pG + base_theme + labs(tag="(h)") + tag_theme
pH_std  <- pH + base_theme + labs(tag="(i)") + tag_theme
pI_std  <- pI + base_theme + labs(tag="(j)") + tag_theme




row1 <- (pA_std | pB_std | pC_std) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

row2 <- row2 <- (pD_std | pE_std | (pF_std | pF1_std))+
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

row3 <- (pG_std | pH_std | pI_std) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


tag_title_fix <- theme(
  plot.title.position = "plot",
  plot.title = element_text(margin = margin(t = 6, b = 8)),  # adds space below title
  plot.tag = element_text(margin = margin(b = 10, r = 6))   # moves tag away from title                              # top-left
)

figure1_1 <- (row1 / row2 / row3) & tag_title_fix
figure1_1


ggsave(file.path(out_dir, "Fig1_FMLE_main.pdf"),
       plot = figure1_1,
       device = cairo_pdf,
       width = 12,
       height = 12,
       units = "in")

