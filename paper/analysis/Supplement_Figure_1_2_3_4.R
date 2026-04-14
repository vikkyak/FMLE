suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(Seurat)
  library(hexbin)
  library(forcats)
  library(scales)
  library(rcompanion) 
})

#========================================================
# Helpers 
#========================================================

files <- c("entropy_error_stats.R", "delta_stats_one.R", "entropy_shift_stats.R", "lineage_pathway.R")
paths <- file.path(here::here(), "paper", "scripts", files)
stopifnot(all(file.exists(paths)))
invisible(lapply(paths, source))

safe_r2 <- function(y, yhat) {
  ok <- is.finite(y) & is.finite(yhat)
  y <- y[ok]; yhat <- yhat[ok]
  if (length(y) < 3) return(NA_real_)
  1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}

entropy_norm_rows <- function(alpha, eps=1e-12) {
  alpha <- as.matrix(alpha)
  rs <- rowSums(alpha)
  stopifnot(all(is.finite(rs)), all(rs > 0))
  alpha <- alpha / rs
  H <- -rowSums(alpha * log(pmax(alpha, eps)))
  H / log(ncol(alpha))
}

jsd <- function(p, q, base = 2, eps = 1e-12) {
  p <- p / sum(p); q <- q / sum(q)
  p <- pmax(p, eps); q <- pmax(q, eps)
  p <- p / sum(p); q <- q / sum(q)
  m <- 0.5 * (p + q)
  logb <- function(x) log(x) / log(base)
  kl <- function(a, b) sum(a * (logb(a) - logb(b)))
  0.5 * kl(p, m) + 0.5 * kl(q, m)
}

celltype_regime_comp <- function(alpha, seu, cell_ids, celltype_col="cell_type", min_n=100, split_label=NULL) {
  alpha <- as.matrix(alpha)
  stopifnot(!is.null(rownames(alpha)))
  stopifnot(all(cell_ids %in% rownames(alpha)))
  stopifnot(all(cell_ids %in% colnames(seu)))
  
  alpha <- alpha[cell_ids, , drop=FALSE]
  seu2  <- seu[, cell_ids]
  
  rs <- rowSums(alpha)
  stopifnot(all(is.finite(rs)), all(rs > 0))
  alpha <- alpha / rs
  
  reg_id <- max.col(alpha, ties.method="first")
  regime <- factor(colnames(alpha)[reg_id], levels = colnames(alpha))
  
  meta <- seu2@meta.data
  stopifnot(celltype_col %in% colnames(meta))
  ct <- as.character(meta[cell_ids, celltype_col])
  ct[is.na(ct) | ct == ""] <- "Unknown"
  ct <- factor(ct)
  
  df <- tibble(cell = cell_ids, celltype = ct, regime = regime)
  
  keep_ct <- names(which(table(df$celltype) >= min_n))
  df <- df %>% filter(as.character(celltype) %in% keep_ct) %>%
    mutate(celltype = droplevels(celltype))
  
  out <- df %>%
    count(celltype, regime, name="n") %>%
    group_by(celltype) %>%
    mutate(frac = n/sum(n)) %>%
    ungroup()
  
  if (!is.null(split_label)) out$split <- split_label
  out
}


#========================================================
# Load FMLE object (TEST-only alpha)
#========================================================

# ===============================================================
# 1) PBMC
# -------------------
# prot <- "CD14"
# rna_gene <- "CD14"
# -------------------
# -------------------
# prot <- "CD27"
# rna_gene <- "CD27"
# -------------------
bench <- 1 # set 1 or 2 or 3 here bench 1 is PBMC Kaggle
ds <- "citeseq_v1"
base <- file.path(cfg$out_root, sprintf("benchmarks_%d", bench))
# ===============================================================



#========================================================
# 2) BMMC
# -------------------
# prot <- "CD14"
# rna_gene <- "CD14"
# -------------------
# -------------------
# prot <- "CD27"
# rna_gene <- "CD27"
# -------------------
bench <- 1  # set 1,  or 2 here bench 1 is BMMC CITE-seq 
ds <- "citeseq_v1"
base <- file.path(cfg$out_root, sprintf("benchmarks_%d", bench))
# ========================================================
  
  
ds_dir <- file.path(base, ds)
out_dir <- file.path(base, "FMLE", paste0(ds, "_final"))
dir.create(ds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seu_prep <- readRDS(file.path(ds_dir, "seu_final.rds")) 

stopifnot(file.exists(file.path(ds_dir, "seu_final.rds")))
stopifnot(file.exists(file.path(out_dir, paste0("final_", prot, ".rds"))))

obj   <- readRDS(file.path(out_dir, paste0("final_", prot, ".rds")))
alpha <- obj$alpha_test
stopifnot(!is.null(alpha), !is.null(rownames(alpha)))

cells <- rownames(alpha)       
hard  <- max.col(alpha)
hard  <- factor(hard)

# regime counts for legend
tab_h <- table(hard)
hard_lab <- factor(hard, levels=levels(hard),
                   labels=paste0(levels(hard), " (n=", as.integer(tab_h[levels(hard)]), ")"))

stopifnot(length(obj$y_test) == length(cells),
          length(obj$yhat_test) == length(cells))

rna_raw <- FetchData(seu_prep, vars=rna_gene, cells=cells)[,1]
ct <- seu_prep$cell_type
names(ct) <- colnames(seu_prep)
celltype_vec <- ct[cells]

df <- data.frame(
  cell = cells,
  rna  = rna_raw,
  y    = as.numeric(obj$y_test),
  yhat = as.numeric(obj$yhat_test),
  hard = hard_lab,
  celltype = as.factor(celltype_vec),     # if cell annotation is available
  stringsAsFactors = FALSE
)

df <- df[is.finite(df$rna) & is.finite(df$y) & is.finite(df$yhat), ]

df <- df[is.finite(df$rna) & is.finite(df$y) & is.finite(df$yhat), ]
stopifnot(!anyNA(df$celltype))

H <- -rowSums(alpha_test * log(pmax(alpha_test, 1e-12))) / log(ncol(alpha_test))

entropy_error_stats(H, df$y, df$yha)

fit_ct <- lm(y ~ rna * celltype, data = df)
summary(fit_ct)
fit_ct_h <- lm(y ~ rna * celltype + hard + rna:hard, data=df)
anova(fit_ct, fit_ct_h)

anova(lm(y ~ hard, data = df))
anova(lm(y ~ rna * hard, data = df))

fit_int <- lm(y ~ rna * hard, data = df)
fit_h_c <- lm(y ~ rna * hard + celltype + rna:celltype, data=df)
anova(fit_int, fit_h_c)

delta_stats_one(df)


p_int <- anova(fit_int)[["Pr(>F)"]][match("rna:hard", rownames(anova(fit_int)))]
p_txt <- paste0("Interaction p=", format.pval(p_int, digits = 2))

p_txt <- if (p_int < 2.2e-16) {
  "Interaction p < 2.2 × 10^-16"
} else {
  paste0("Interaction p = ", format(signif(p_int, 3), scientific = TRUE))
}


#========================================================
# PANEL A — Regime-dependent RNA–protein coupling
#========================================================

pA <- ggplot(df, aes(rna, y, color = hard)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  theme_classic(base_size = 8) +
  labs(
    x = paste0(rna_gene, " RNA (log1p)"),
    y = paste0(prot, " protein"),
    color = "regime",
    title = "Regime-dependent RNA–protein coupling"
  ) +
  annotate("text", x = Inf, y = Inf, label = p_txt,
           hjust = 1.1, vjust = 1.2, size = 3.6) +
  coord_cartesian(xlim = c(0, 3.8), expand = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#========================================================
# PANEL B — global failure 
#========================================================
m_global <- lm(y ~ rna, data=df)
df$resid_global <- df$y - predict(m_global, newdata=df)

pB <- ggplot(df, aes(rna, resid_global, color = hard)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.7) +
  theme_classic(base_size = 8) +
  labs(
    x = paste0(rna_gene, " RNA (log1p)"),
    y = "Global model residual",
    title = "Residuals from a global mapping"
  ) +
  guides(color = "none") +
  coord_cartesian(xlim = c(0, 3.8), expand = FALSE) +  
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

#===============================================================
# PANEL C D E and F — Global, FMLE, cTPnet and scLiear  vs truth 
#===============================================================
# global prediction from same model (optional overlay)
df$yhat_global <- predict(m_global, newdata=df)

res_ctpnet <- readr::read_csv(file.path(base, "ctp", "ctpnet_test_metrics_FMLEscale.csv"))
pred_ctpnet <- readr::read_csv(file.path(base, "ctp", "ctpnet_test_predictions_FMLEscale.csv"))
res_sclinear  <- readr::read_csv(file.path(base, "sclinear", "scLinear_test_metrics_fmle_fair.csv"))
pred_sclinear <- readr::read_csv(file.path(base, "sclinear", "scLinear_test_predictions_fmle_fair.csv"))

cell_ids <- if ("cell" %in% colnames(df)) df$cell else rownames(df)
tmp_ctp <- pred_ctpnet %>%
  filter(protein == prot)

ctp_vec <- tmp_ctp$yhat_test_cal
names(ctp_vec) <- tmp_ctp$cell
df$ctp_yhat <- ctp_vec[cell_ids]

tmp_scl <- pred_sclinear %>%
  filter(protein == prot)
scl_vec <- tmp_scl$yhat_test_cal
names(scl_vec) <- tmp_scl$cell

df$sclinear_yhat <- scl_vec[cell_ids]

summary(df$ctp_yhat)
summary(df$sclinear_yhat)
sum(is.na(df$ctp_yhat))
sum(is.na(df$sclinear_yhat))
df$yhat_global <- predict(m_global, newdata=df)

df_plotD <- bind_rows(
  df %>% transmute(y, pred = yhat_global,   model = "Global"),
  df %>% transmute(y, pred = yhat,          model = "FMLE"),
  df %>% transmute(y, pred = ctp_yhat,      model = "cTPnet"),
  df %>% transmute(y, pred = sclinear_yhat, model = "scLinear")
) %>%
  filter(is.finite(y), is.finite(pred))

lims <- range(c(df_plotD$y, df_plotD$pred), finite = TRUE)
pad  <- 0.02 * diff(lims)
lims <- lims + c(-pad, pad)

metricsD <- df_plotD %>%
  group_by(model) %>%
  summarise(
    r   = suppressWarnings(cor(y, pred, method = "pearson",  use = "complete.obs")),
    rho = suppressWarnings(cor(y, pred, method = "spearman", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(
    lbl = sprintf("r=%.3f; \u03C1=%.3f", r, rho),
    x = lims[2], y = lims[1]
  )

base_theme <- theme_classic(base_size = 8) +
  theme(plot.title = element_text(hjust=0.5, face="bold"),
        legend.position = "none")
cols <- c(
  Global   = "#0072B2",
  FMLE     = "#D55E00",
  cTPnet   = "#009E73",
  scLinear = "#CC79A7"
)
make_panel_pred <- function(m){
  ggplot(df_plotD %>% filter(model == m), aes(x = y, y = pred)) +
    geom_hex(bins = 55) +
    scale_fill_gradient(low = "grey80", high = cols[[m]], trans = "sqrt", guide = "none") +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.8, color="grey35") +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE) +
    labs(title = m, x = "True protein (test data)", y = "Prediction") +
    geom_text(
      data = metricsD %>% filter(model == m),
      aes(x = x, y = y, label = lbl),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = -3.0, size = 3.4,
      fontface = "bold"
    ) +
    base_theme
}

pC <- make_panel_pred("Global")
pD <- make_panel_pred("FMLE")
p_ctpnet   <- make_panel_pred("cTPnet")
p_sclinear <- make_panel_pred("scLinear")


#========================================================
# UMAP df for train/test only (aligned set)
#========================================================

umap_all <- Embeddings(seu_prep, "umap")
stopifnot(!is.null(rownames(umap_all)))

cells_keep <- intersect(rownames(umap_all), c(train_cells, test_cells))
stopifnot(length(cells_keep) > 0)

umap_use <- umap_all[cells_keep, , drop=FALSE]

data_split <- rep(NA_character_, length(cells_keep))
names(data_split) <- cells_keep
data_split[names(data_split) %in% train_cells] <- "train_data"
data_split[names(data_split) %in% test_cells]  <- "test_data"
stopifnot(!anyNA(data_split))

umap_df <- data.frame(
  cell     = cells_keep,
  UMAP_1   = umap_use[,1],
  UMAP_2   = umap_use[,2],
  donor    = data_split,
  celltype = seu_prep$cell_type[cells_keep],  
  stringsAsFactors = FALSE
)

fit <- obj$fit
pred_tr <- FMLE::fmle_predict(fit, X_new = X_train, Z_new = Z_train, return_se = TRUE)
pred_te <- FMLE::fmle_predict(fit, X_new = X_test,  Z_new = Z_test,  return_se = TRUE)


alpha_te <- pred_te$alpha
alpha_tr <- pred_tr$alpha
if (is.null(rownames(alpha_te))) rownames(alpha_te) <- rownames(X_test)
if (is.null(rownames(alpha_tr))) rownames(alpha_tr) <- rownames(X_train)

stopifnot(setequal(rownames(alpha_te), test_cells))
stopifnot(setequal(rownames(alpha_tr), train_cells))

reg_te <- max.col(alpha_te, ties.method="first")
Hn_te  <- entropy_norm_rows(alpha_te)

tapply(Hn_te, hard, summary)
tapply(Hn_te, seu_prep$cell_type[cells], summary)

Hn_tr <- entropy_norm_rows(alpha_tr)
names(Hn_tr) <- rownames(alpha_tr)

hard_tr <- max.col(alpha_tr, ties.method="first")
names(hard_tr) <- rownames(alpha_tr)

ct_tr <- seu_prep@meta.data[rownames(alpha_tr), "cell_type", drop=TRUE]

tapply(Hn_tr[names(hard_tr)], factor(hard_tr), summary)
tapply(Hn_tr[rownames(alpha_tr)], factor(ct_tr), summary)

summary(Hn_tr)
summary(Hn_te)

entropy_shift_stats(Hn_tr, Hn_te)

pred_df_te <- data.frame(
  cell    = rownames(alpha_te),
  regime  = factor(reg_te),
  entropy = as.numeric(Hn_te),
  stringsAsFactors = FALSE
)

plot_df_te <- umap_df %>%
  filter(donor == "test_data") %>%
  left_join(pred_df_te, by="cell")

stopifnot(nrow(plot_df_te) == length(test_cells))
stopifnot(!anyNA(plot_df_te$regime), !anyNA(plot_df_te$entropy))

#========================================================
# PANEL G
#========================================================
p1 <- ggplot(plot_df_te, aes(UMAP_1, UMAP_2, color = regime)) +
  geom_point(size = 0.25, alpha = 0.7) +
  theme_classic(base_size = 8) +
  labs(title = "Test data UMAP: regime (argmax α)", color = "Regime")

#========================================================
# PANEL H
#========================================================

p2 <- ggplot(plot_df_te, aes(UMAP_1, UMAP_2, color = entropy)) +
  geom_point(size = 0.25, alpha = 0.7) +
  theme_classic(base_size = 8) +
  labs(title = "Test data UMAP: α-entropy", color = "Entropy")

# p3: regime composition train vs test + JSD
reg_tr <- max.col(alpha_tr, ties.method="first")

pred_df_tr <- data.frame(
  cell   = rownames(alpha_tr),
  regime = factor(reg_tr),
  stringsAsFactors = FALSE
)

comp_cells <- umap_df %>%
  filter(donor %in% c("train_data","test_data")) %>%
  select(cell, donor) %>%
  left_join(bind_rows(pred_df_tr, pred_df_te %>% select(cell, regime)), by="cell")

stopifnot(!anyNA(comp_cells$regime))

comp <- comp_cells %>%
  count(donor, regime, name="n") %>%
  group_by(donor) %>%
  mutate(frac = n/sum(n)) %>%
  ungroup()

# JSD on matched regime support
regs <- sort(unique(comp$regime))
pvec <- comp %>% filter(donor=="train_data") %>%
  right_join(data.frame(regime=regs), by="regime") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac)) %>% arrange(regime) %>% pull(frac)
qvec <- comp %>% filter(donor=="test_data") %>%
  right_join(data.frame(regime=regs), by="regime") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac)) %>% arrange(regime) %>% pull(frac)

JSD <- jsd(pvec, qvec, base=2)
cat(sprintf("JSD(train,test)=%.6g\n", JSD))


min_n <- 100
ct_col <- "cell_type" # for kaggle datasets

comp_tr_ct <- celltype_regime_comp(alpha_tr, seu_prep, cell_ids=train_cells,
                                   celltype_col=ct_col, min_n=min_n, split_label="train")
comp_te_ct <- celltype_regime_comp(alpha_te, seu_prep, cell_ids=test_cells,
                                   celltype_col=ct_col, min_n=min_n, split_label="test")
comp_all_ct <- bind_rows(comp_tr_ct, comp_te_ct)

p4_df <- comp_all_ct %>%
  dplyr::filter(frac > 0) %>%
  droplevels()

#========================================================
# PANEL M
#========================================================

p4 <- ggplot(p4_df, aes(x = celltype, y = frac, fill = regime)) +
  geom_col(width = 0.85) +
  facet_wrap(~split, nrow = 1, scales = "free_x") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey90", colour = "black")
  ) +
  labs(title = "Regime composition within \ncell types", x = NULL, y = "Fraction")

make_ctab <- function(alpha, seu, cell_ids, celltype_col="cell_type", min_n=100) {
  alpha <- as.matrix(alpha)[cell_ids, , drop=FALSE]
  reg   <- colnames(alpha)[max.col(alpha, ties.method="first")]
  ct    <- as.character(seu@meta.data[cell_ids, celltype_col])
  ct[is.na(ct) | ct==""] <- "Unknown"
  df <- tibble(celltype=ct, regime=reg)
  keep_ct <- names(which(table(df$celltype) >= min_n))
  df <- df %>% filter(celltype %in% keep_ct)
  table(df$celltype, df$regime)
}

xt_tr <- make_ctab(alpha_tr, seu_prep, train_cells, celltype_col=ct_col, min_n=min_n)
xt_te <- make_ctab(alpha_te, seu_prep, test_cells,  celltype_col=ct_col, min_n=min_n)

cat("\nTRAIN: chisq + Cramer's V\n")
print(chisq.test(xt_tr, simulate.p.value = TRUE, B = 20000))
print(cramerV(xt_tr))

cat("\nTEST: chisq + Cramer's V\n")
print(chisq.test(xt_te, simulate.p.value = TRUE, B = 20000))
print(cramerV(xt_te))


# =============================
# PBMC CD14
# =============================
lineage_use <- "Monocyte"
meta <- seu_prep@meta.data
mono_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_mono <- as.matrix(alpha_te[mono_cells, , drop = FALSE])

hard_mono <- colnames(A_mono)[max.col(A_mono, ties.method = "first")]
tab_mono  <- sort(table(hard_mono), decreasing = TRUE)
score_regimes <- names(tab_mono)[1:2]

A2 <- A_mono[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)

umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]

mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)

lims <- quantile(mono_df$score, c(0.02, 0.98), na.rm = TRUE)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)

xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)

p3_left <- ggplot() +
  geom_point(
    data  = mono_df,
    aes(UMAP_1, UMAP_2),
    color = "grey82",
    size  = 0.6,
    alpha = 0.7
  ) +
  geom_point(
    data  = mono_df,
    aes(UMAP_1, UMAP_2, color = score),
    size  = 0.5,
    alpha = 0.9
  ) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1,
    vjust       = 1.5,
    size        = 3.5,
    color       = "black",
    fontface    = "bold"
  ) +
  scale_color_gradientn(
    colors = c("#FEE0D2", "#FC9272", "#DE2D26", "#99000D"),
    limits = lims, oob = scales::squish,
    name   = "FMLE\nscore",
    guide  = guide_colorbar(barwidth = 0.5, barheight = 3)
  ) +
  coord_cartesian(
    xlim = xr + c(-0.1, 0.1),
    ylim = yr + c(-4, 4)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "Monocytes:\nclassical-to-FCGR3A-associated axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 8),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 8),
    plot.margin   = margin(3, 3, 3, 3)
  )
p3_left

# =============================

tag_theme <- theme(plot.tag = element_text(face="bold"))


pA <- pA + labs(tag="(a)") + tag_theme
pB <- pB + labs(tag="(b)") + tag_theme
pC <- pC + labs(tag="(c)") + tag_theme
p_ctpnet  <- p_ctpnet + labs(tag="(d)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(e)") + tag_theme
pD <- pD + labs(tag="(f)") + tag_theme
p1 <- p1 + labs(tag="(g)") + tag_theme
p2 <- p2 + labs(tag="(h)") + tag_theme
p3_left <- p3_left + labs(tag="(i)") + tag_theme
p4 <- p4 + labs(tag="(j)") + tag_theme

fig <- 
  (pA | pB | pC | p_ctpnet |p_sclinear) /
  (pD | p1 | p2 | p3_left | p4) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    legend.position = "right"
  ) 

fig

ggsave(file.path(out_dir, "Supplementary_Figure_1.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 13, height = 7, units = "in")

# ==============================================================================
# CD27 PBMC
# ==============================================================================
lineage_use <- "B" # or NK in CD27 PBMC
ct_col <- "cell_type"

lineages_to_test <- sort(unique(seu_prep@meta.data[[ct_col]]))
lineage_tab <- bind_rows(lapply(lineages_to_test, screen_one))
print(lineage_tab, n= Inf)

res_nk <- lineage_stats(
  seu_prep = seu_prep,
  alpha_tr = alpha_tr,
  alpha_te = alpha_te,
  train_cells = train_cells,
  test_cells = test_cells,
  ct_col = ct_col,
  lineage = lineage_use,
  score_regimes = NULL,  
  assay = "RNA", 
  layer = "data",
  min_cells = 200
)

res_nk$de_tr$summary
res_nk$de_te$summary
res_nk$hard_repro

res_nk$assoc_tr$summary
res_nk$assoc_te$summary
res_nk$cont_repro

head(res_nk$de_te$res, 20)
head(res_nk$assoc_te$assoc, 20)

res_mono <- res_nk   # because your current object holds Monocyte results

rank_te <- make_ranked_list(res_mono$assoc_te$assoc)
rank_tr <- make_ranked_list(res_mono$assoc_tr$assoc)

gsea_h_te  <- run_gsea_msig(rank_te, category = "H")
gsea_h_tr  <- run_gsea_msig(rank_tr, category = "H")

gsea_re_te <- run_gsea_reactome(rank_te)
gsea_re_tr <- run_gsea_reactome(rank_tr)

# top pathways
top_h_te  <- top_gsea_terms(gsea_h_te, 15)

top_h_tr  <- top_gsea_terms(gsea_h_tr, 15)

top_re_te <- top_gsea_terms(gsea_re_te, 15)
top_re_tr <- top_gsea_terms(gsea_re_tr, 15)

# concordant train-test pathways
hallmark_concordant <- concordant_gsea(gsea_h_tr, gsea_h_te, padj_cut = 0.05)
reactome_concordant <- concordant_gsea(gsea_re_tr, gsea_re_te, padj_cut = 0.05)

# inspect
top_h_te
top_re_te
hallmark_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)
reactome_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)
meta <- seu_prep@meta.data
mono_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_mono <- as.matrix(alpha_te[mono_cells, , drop = FALSE])

hard_mono <- colnames(A_mono)[max.col(A_mono, ties.method = "first")]
tab_mono  <- sort(table(hard_mono), decreasing = TRUE)
score_regimes <- names(tab_mono)[1:2]

A2 <- A_mono[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)

umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]

mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)
n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)

q <- quantile(mono_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)
xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)
p3_left <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 1, alpha = 1) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.5, vjust = 1.2,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#EAEAEA",
    high = "#B2182B",
    midpoint = 0,
    limits = lims,
    oob = scales::squish,
    name = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  coord_cartesian(
    xlim = xr + c(3, 7),
    ylim = yr + c(-0.5, 1)
  ) +
  labs(
    title = "B cells: naive-to-antigen-presentation-\nassociated axis",   # fix title
    subtitle = "PBMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right",
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 8),
    plot.margin   = margin(3, 3, 3, 3)
  )

p3_left
p3_left2 <- p3_left

lineage_use <- "NK cells"

mono_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_mono <- as.matrix(alpha_te[mono_cells, , drop = FALSE])

hard_mono <- colnames(A_mono)[max.col(A_mono, ties.method = "first")]
tab_mono  <- sort(table(hard_mono), decreasing = TRUE)
score_regimes <- names(tab_mono)[1:2]

A2 <- A_mono[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)

umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]

mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)
n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)

q <- quantile(mono_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)
xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)
p3_nk <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 1.0,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#EAEAEA",
    high = "#B2182B",
    midpoint = 0,
    limits = lims,
    oob = scales::squish,
    name = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  coord_cartesian(
    xlim = xr + c(-0.5, 0.5),
    ylim = yr + c(-0.5, 1)
  ) +
  labs(
    title = "NK cells:\nIL7R/KLRB1-to-cytotoxic-associated axis",   # fix title
    subtitle = "PBMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right",
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 8),
    plot.margin   = margin(3, 3, 3, 3)
  )

p3_nk

pathway_df_nk <- reactome_concordant %>%
  filter(ID %in% c(
    "REACTOME_TRANSLATION",
    "REACTOME_METABOLISM_OF_RNA",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_VIRAL_INFECTION_PATHWAYS"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "REACTOME_TRANSLATION" ~ "Translation",
      ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
      ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
      ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
      ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD",
      ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
      ID == "REACTOME_VIRAL_INFECTION_PATHWAYS" ~ "Viral infection pathways"
    ),
    NES     = NES_te,
    padj    = padj_te,
    fdr_lab = sprintf("FDR %.1e", padj_te)
  ) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))


xmin <- min(pathway_df_nk$NES, na.rm = TRUE)
xmax <- max(pathway_df_nk$NES, na.rm = TRUE)

p3_right_nk <- ggplot(pathway_df_nk, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x = NES / 2
    ),
    hjust = 0.5,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(-1.60, 1.740),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title    = "Reproducible NK-cell biosynthetic\npathway shifts",
    subtitle = "Concordant in train and test; bars show test NES",
    x        = "NES",
    y        = NULL
  ) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_nk



tag_theme <- theme(plot.tag = element_text(face="bold"))


pA <- pA + labs(tag="(a)") + tag_theme
pB <- pB + labs(tag="(b)") + tag_theme
pC <- pC + labs(tag="(c)") + tag_theme
p_ctpnet  <- p_ctpnet + labs(tag="(d)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(e)") + tag_theme
pD <- pD + labs(tag="(f)") + tag_theme
p1 <- p1 + labs(tag="(g)") + tag_theme
p2 <- p2 + labs(tag="(h)") + tag_theme
p4 <- p4 + labs(tag="(i)") + tag_theme

p3_left2 <- p3_left2 + labs(tag = "(j)") + tag_theme
p3_nk  <- p3_nk  + labs(tag = "k") + tag_theme  
p3_right_nk <- p3_right_nk + labs(tag = "l") + tag_theme  



fig <- 
  (pA | pB | pC | p_ctpnet) /
  (p_sclinear| pD | p1 | p2) /
  ( p4| p3_left2 |p3_nk|p3_right_nk )+ 
  plot_layout(heights = c(1, 1, 1)) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    legend.position = "right")

fig


ggsave(file.path(out_dir, "Supplementary_Figure_3.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 12.6, height = 8, units = "in")






# ==============================================================================
#  CD 14 BMMC
# ==============================================================================
meta <- seu_prep@meta.data
lineage_use  <- "CD14+ Mono" 

ct_col <- "cell_type"

lineages_to_test <- sort(unique(seu_prep@meta.data[[ct_col]]))
lineage_tab <- bind_rows(lapply(lineages_to_test, screen_one))
print(lineage_tab, n= Inf)

res_nk <- lineage_stats(
  seu_prep = seu_prep,
  alpha_tr = alpha_tr,
  alpha_te = alpha_te,
  train_cells = train_cells,
  test_cells = test_cells,
  ct_col = ct_col,
  lineage = lineage_use,
  score_regimes = NULL,  
  assay = "RNA", 
  layer = "data",
  min_cells = 200
)

res_nk$de_tr$summary
res_nk$de_te$summary
res_nk$hard_repro

res_nk$assoc_tr$summary
res_nk$assoc_te$summary
res_nk$cont_repro

head(res_nk$de_te$res, 20)
head(res_nk$assoc_te$assoc, 20)

res_mono <- res_nk   # because your current object holds Monocyte results

rank_te <- make_ranked_list(res_mono$assoc_te$assoc)
rank_tr <- make_ranked_list(res_mono$assoc_tr$assoc)

gsea_h_te  <- run_gsea_msig(rank_te, category = "H")
gsea_h_tr  <- run_gsea_msig(rank_tr, category = "H")

gsea_re_te <- run_gsea_reactome(rank_te)
gsea_re_tr <- run_gsea_reactome(rank_tr)

# top pathways
top_h_te  <- top_gsea_terms(gsea_h_te, 15)

top_h_tr  <- top_gsea_terms(gsea_h_tr, 15)

top_re_te <- top_gsea_terms(gsea_re_te, 15)
top_re_tr <- top_gsea_terms(gsea_re_tr, 15)

# concordant train-test pathways
hallmark_concordant <- concordant_gsea(gsea_h_tr, gsea_h_te, padj_cut = 0.05)
reactome_concordant <- concordant_gsea(gsea_re_tr, gsea_re_te, padj_cut = 0.05)

# inspect
top_h_te
top_re_te
hallmark_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)
reactome_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)



CD14_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)
A_cd14 <- as.matrix(alpha_te[CD14_cells, , drop = FALSE])
hard_cd14 <- colnames(A_cd14)[max.col(A_cd14, ties.method = "first")]
tab_cd14  <- sort(table(hard_cd14), decreasing = TRUE)
score_regimes_cd14 <- names(tab_cd14)[1:2]
A2 <- A_cd14[, score_regimes_cd14, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)
mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)
umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]
mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)
m <- max(abs(quantile(mono_df$score, c(0.02, 0.98), na.rm = TRUE)))
lims <- c(-m, m)
n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)
xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)
p3_left <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 1.5,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradientn(
    colors = c("#FEE0D2", "#FC9272", "#DE2D26", "#99000D"),
    limits = lims, oob = scales::squish,
    name   = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.4, 0.4),
    ylim = yr + c(-0.4, 0.4)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD14\u207a monocytes:\nTNF/NF\u03baB-to-IFN\u03b1-associated axis",
    subtitle = "BMMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    # plot.subtitle = element_text(size = 7, color = "grey40", face = "bold"),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 8),
    axis.line     = element_line(linewidth = 0.3),
    plot.margin   = margin(3, 3, 3, 3)
  )
p3_left_14 <- p3_left
p3_left_14

pathway_df_mono <- hallmark_concordant %>%
  filter(ID %in% c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"   ~ "Hallmark TNF/NF\u03baB",
      ID == "HALLMARK_INTERFERON_ALPHA_RESPONSE" ~ "Hallmark IFN-\u03b1"
    ),
    NES     = NES_te,
    padj    = padj_te,
    fdr_lab = sprintf("FDR %.1e", padj_te)
  ) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_mono$NES, na.rm = TRUE)
xmax <- max(pathway_df_mono$NES, na.rm = TRUE)

p3_right_mono <- ggplot(pathway_df_mono, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x = NES / 2
    ),
    hjust = 0.5,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(xmin + 0.08, xmax + 0.08),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title    = "Reproducible CD14\u207a monocyte\nTNF/NF\u03baB and IFN-\u03b1 pathways",
    subtitle = "Concordant in train and test; bars show test NES",
    x        = "NES",
    y        = NULL
  ) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_mono

# ======================
# CD16+ Mono
# ======================
lineage_use <- "CD16+ Mono"

ct_col <- "cell_type"

res_nk <- lineage_stats(
  seu_prep = seu_prep,
  alpha_tr = alpha_tr,
  alpha_te = alpha_te,
  train_cells = train_cells,
  test_cells = test_cells,
  ct_col = ct_col,
  lineage = lineage_use,
  score_regimes = NULL,  
  assay = "RNA", 
  layer = "data",
  min_cells = 200
)

res_nk$de_tr$summary
res_nk$de_te$summary
res_nk$hard_repro

res_nk$assoc_tr$summary
res_nk$assoc_te$summary
res_nk$cont_repro

head(res_nk$de_te$res, 20)
head(res_nk$assoc_te$assoc, 20)

res_mono <- res_nk   # because your current object holds Monocyte results

rank_te <- make_ranked_list(res_mono$assoc_te$assoc)
rank_tr <- make_ranked_list(res_mono$assoc_tr$assoc)

gsea_h_te  <- run_gsea_msig(rank_te, category = "H")
gsea_h_tr  <- run_gsea_msig(rank_tr, category = "H")

gsea_re_te <- run_gsea_reactome(rank_te)
gsea_re_tr <- run_gsea_reactome(rank_tr)

# top pathways
top_h_te  <- top_gsea_terms(gsea_h_te, 15)

top_h_tr  <- top_gsea_terms(gsea_h_tr, 15)

top_re_te <- top_gsea_terms(gsea_re_te, 15)
top_re_tr <- top_gsea_terms(gsea_re_tr, 15)

# concordant train-test pathways
hallmark_concordant <- concordant_gsea(gsea_h_tr, gsea_h_te, padj_cut = 0.05)
reactome_concordant <- concordant_gsea(gsea_re_tr, gsea_re_te, padj_cut = 0.05)

# inspect
top_h_te
top_re_te
hallmark_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)
reactome_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)


CD16_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)
A_cd16 <- as.matrix(alpha_te[CD16_cells, , drop = FALSE])
hard_cd16 <- colnames(A_cd16)[max.col(A_cd16, ties.method = "first")]
tab_cd16  <- sort(table(hard_cd16), decreasing = TRUE)
score_regimes_cd16 <- c("1", "3")
A2 <- A_cd16[, score_regimes_cd16, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)
mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)

# ── MISSING — rebuild mono_df with new score ──
umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]
mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)
# ── MISSING — rebuild lims and n_label from new mono_df ──
m <- quantile(abs(mono_df$score), 0.98, na.rm = TRUE)
lims <- c(-m, m)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)

xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)

mono_df2 <- mono_df[order(abs(mono_df$score)), ]

p3_left <- ggplot(mono_df2, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.9, alpha = 0.9) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 1.5,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "#F2F2F2",
    high = "#B2182B",
    midpoint = 0,
    limits = lims,
    oob = scales::squish,
    breaks = pretty(lims, n = 5),
    name = "FMLE\nscore"
  )  +
  coord_cartesian(
    xlim = xr + c(-0.4, 0.4),
    ylim = yr + c(-0.5, 0.5)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD16\u207a monocytes:\nTLR/NF\u03baB-associated axis",
    subtitle = "BMMC, test cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 8),
    plot.margin  = margin(3, 3, 3, 3)
  )

p3_left_16 <- p3_left
p3_left_16

pathway_df_mono <- bind_rows(
  hallmark_concordant %>%
    filter(ID %in% c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"   ~ "Hallmark TNF/NF\u03baB",
        ID == "HALLMARK_INTERFERON_ALPHA_RESPONSE" ~ "Hallmark IFN-\u03b1"
      ),
      NES     = NES_te,
      padj    = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    ),
  
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES",
      "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
      "REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE",
      "REACTOME_TOLL_LIKE_RECEPTOR_9_TLR9_CASCADE"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES"              ~ "TLR cascades",
        ID == "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE"    ~ "TLR1/TLR2 cascade",
        ID == "REACTOME_MYD88_INDEPENDENT_TLR4_CASCADE"          ~ "TLR4 cascade",
        ID == "REACTOME_TOLL_LIKE_RECEPTOR_9_TLR9_CASCADE"       ~ "TLR9 cascade"
      ),
      NES     = NES_te,
      padj    = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    )
) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_mono$NES, na.rm = TRUE)
xmax <- max(pathway_df_mono$NES, na.rm = TRUE)

p3_right_mono_16 <- ggplot(pathway_df_mono, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x = NES / 2
    ),
    hjust = 0.5,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(xmin + 0.15, xmax + 0.08),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title    = "Reproducible CD16\u207a monocyte TLR/NF\u03baB\nand interferon-associated pathways",
    subtitle = "Concordant in train and test; bars show test NES",
    x        = "NES",
    y        = NULL
  ) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_mono_16



tag_theme <- theme(plot.tag = element_text(face="bold"))


pA <- pA + labs(tag="(a)") + tag_theme
pB <- pB + labs(tag="(b)") + tag_theme
pC <- pC + labs(tag="(c)") + tag_theme
p_ctpnet  <- p_ctpnet + labs(tag="(d)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(e)") + tag_theme
pD <- pD + labs(tag="(f)") + tag_theme
p1 <- p1 + labs(tag="(g)") + tag_theme
p2 <- p2 + labs(tag="(h)") + tag_theme
p3_left_14 <- p3_left_14 + labs(tag="(i)") + tag_theme
p3_right <- p3_right_mono + labs(tag="(j)") + tag_theme
p3_left_16 <- p3_left_16 + labs(tag="(k)") + tag_theme
p3_right_mono_16 <- p3_right_mono_16 + labs(tag="(l)") + tag_theme
p4 <- p4 + labs(tag="(m)") + tag_theme
p4 <- p4 +
  facet_wrap(~split, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  )

fig <- 
  
  (pA | pB | pC | p_ctpnet)/
  (p_sclinear|pD |p1 | p2) /
  (p3_left_14 |p3_right | p3_left_16|p3_right_mono_16) /
  p4 + 
  plot_layout(heights = c(1.2, 1.2, 1.3, 1.2)) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    legend.position = "right"
  ) 

fig

ggsave(file.path(out_dir, "Supplementary_Figure_2.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 13.2, height = 11, units = "in")





# ==============================================================================
# CD27 BMMC
# ==============================================================================
meta <- seu_prep@meta.data
lineage_use <- "Naive CD20+ B IGKC+"
ct_col <- "cell_type"

lineages_to_test <- sort(unique(seu_prep@meta.data[[ct_col]]))
lineage_tab <- bind_rows(lapply(lineages_to_test, screen_one))
print(lineage_tab, n= Inf)

res_nk <- lineage_stats(
  seu_prep = seu_prep,
  alpha_tr = alpha_tr,
  alpha_te = alpha_te,
  train_cells = train_cells,
  test_cells = test_cells,
  ct_col = ct_col,
  lineage = lineage_use,
  score_regimes = NULL,  
  assay = "RNA", 
  layer = "data",
  min_cells = 200
)

res_nk$de_tr$summary
res_nk$de_te$summary
res_nk$hard_repro

res_nk$assoc_tr$summary
res_nk$assoc_te$summary
res_nk$cont_repro

head(res_nk$de_te$res, 20)
head(res_nk$assoc_te$assoc, 20)

res_mono <- res_nk   # because your current object holds Monocyte results

rank_te <- make_ranked_list(res_mono$assoc_te$assoc)
rank_tr <- make_ranked_list(res_mono$assoc_tr$assoc)

gsea_h_te  <- run_gsea_msig(rank_te, category = "H")
gsea_h_tr  <- run_gsea_msig(rank_tr, category = "H")

gsea_re_te <- run_gsea_reactome(rank_te)
gsea_re_tr <- run_gsea_reactome(rank_tr)

# top pathways
top_h_te  <- top_gsea_terms(gsea_h_te, 15)

top_h_tr  <- top_gsea_terms(gsea_h_tr, 15)

top_re_te <- top_gsea_terms(gsea_re_te, 15)
top_re_tr <- top_gsea_terms(gsea_re_tr, 15)

# concordant train-test pathways
hallmark_concordant <- concordant_gsea(gsea_h_tr, gsea_h_te, padj_cut = 0.05)
reactome_concordant <- concordant_gsea(gsea_re_tr, gsea_re_te, padj_cut = 0.05)

# inspect
top_h_te
top_re_te
hallmark_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)
reactome_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)

CD14_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)
A_cd14 <- as.matrix(alpha_te[CD14_cells, , drop = FALSE])
hard_cd14 <- colnames(A_cd14)[max.col(A_cd14, ties.method = "first")]
tab_cd14  <- sort(table(hard_cd14), decreasing = TRUE)
score_regimes_cd14 <- names(tab_cd14)[1:2]
A2 <- A_cd14[, score_regimes_cd14, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)
mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)
umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]
mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)
m <- max(abs(quantile(mono_df$score, c(0.02, 0.98), na.rm = TRUE)))
lims <- c(-m, m)
n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)
xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)
p3_left <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.2, vjust = 1.5,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradientn(
    colors = c("#FEE0D2", "#FC9272", "#DE2D26", "#99000D"),
    limits = lims, oob = scales::squish,
    name   = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.7, 0.4),
    ylim = yr + c(-0.4, 0.4)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "Naive CD20\u207a B cells:\nstress/inflammatory-to-biosynthetic axis",
    subtitle = "BMMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    # plot.subtitle = element_text(size = 7, color = "grey40", face = "bold"),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 8),
    axis.line     = element_line(linewidth = 0.3),
    plot.margin   = margin(3, 3, 3, 3)
  )
p3_left_14 <- p3_left
p3_left_14


pathway_df_B <- bind_rows(
  hallmark_concordant %>%
    filter(ID %in% c(
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_HYPOXIA",
      "HALLMARK_APOPTOSIS",
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_INFLAMMATORY_RESPONSE"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "HALLMARK_TNFA_SIGNALING_VIA_NFKB" ~ "TNF/NF\u03baB signaling",
        ID == "HALLMARK_HYPOXIA" ~ "Hypoxia",
        ID == "HALLMARK_APOPTOSIS" ~ "Apoptosis",
        ID == "HALLMARK_INTERFERON_GAMMA_RESPONSE" ~ "IFN-\u03b3 response",
        ID == "HALLMARK_INFLAMMATORY_RESPONSE" ~ "Inflammatory response"
      ),
      NES     = NES_te,
      padj    = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    ),
  
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response"
      ),
      NES     = NES_te,
      padj    = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    )
) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_B$NES, na.rm = TRUE)
xmax <- max(pathway_df_B$NES, na.rm = TRUE)

p3_right_B <- ggplot(pathway_df_B, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(label = fdr_lab, x = NES / 2),
    hjust = 0.5,
    size = 2.5,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(xmin + 0.15, xmax + 0.12),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "Reproducible naive CD20\u207a B-cell stress/\ninflammatory and biosynthetic pathways",
    subtitle = "Concordant in train and test; bars show test NES",
    x        = "NES",
    y        = NULL
  ) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_B




lineage_use <- "CD4+ T activated"
ct_col <- "cell_type"

res_nk <- lineage_stats(
  seu_prep = seu_prep,
  alpha_tr = alpha_tr,
  alpha_te = alpha_te,
  train_cells = train_cells,
  test_cells = test_cells,
  ct_col = ct_col,
  lineage = lineage_use,
  score_regimes = NULL,  
  assay = "RNA", 
  layer = "data",
  min_cells = 200
)

res_nk$de_tr$summary
res_nk$de_te$summary
res_nk$hard_repro

res_nk$assoc_tr$summary
res_nk$assoc_te$summary
res_nk$cont_repro

head(res_nk$de_te$res, 20)
head(res_nk$assoc_te$assoc, 20)

res_mono <- res_nk   # because your current object holds Monocyte results

rank_te <- make_ranked_list(res_mono$assoc_te$assoc)
rank_tr <- make_ranked_list(res_mono$assoc_tr$assoc)

gsea_h_te  <- run_gsea_msig(rank_te, category = "H")
gsea_h_tr  <- run_gsea_msig(rank_tr, category = "H")

gsea_re_te <- run_gsea_reactome(rank_te)
gsea_re_tr <- run_gsea_reactome(rank_tr)

# top pathways
top_h_te  <- top_gsea_terms(gsea_h_te, 15)

top_h_tr  <- top_gsea_terms(gsea_h_tr, 15)

top_re_te <- top_gsea_terms(gsea_re_te, 15)
top_re_tr <- top_gsea_terms(gsea_re_tr, 15)

# concordant train-test pathways
hallmark_concordant <- concordant_gsea(gsea_h_tr, gsea_h_te, padj_cut = 0.05)
reactome_concordant <- concordant_gsea(gsea_re_tr, gsea_re_te, padj_cut = 0.05)

# inspect
top_h_te
top_re_te
hallmark_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)
reactome_concordant %>% select(Description, NES_tr, padj_tr, NES_te, padj_te) %>% head(20)

CD14_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)
A_cd14 <- as.matrix(alpha_te[CD14_cells, , drop = FALSE])
hard_cd14 <- colnames(A_cd14)[max.col(A_cd14, ties.method = "first")]
tab_cd14  <- sort(table(hard_cd14), decreasing = TRUE)
score_regimes_cd14 <- names(tab_cd14)[1:2]
A2 <- A_cd14[, score_regimes_cd14, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)
mono_score <- A2[, 1] - A2[, 2]
names(mono_score) <- rownames(A2)
umap_mono <- Embeddings(seu_prep, "umap")[names(mono_score), , drop = FALSE]
mono_df <- data.frame(
  cell   = names(mono_score),
  UMAP_1 = umap_mono[, 1],
  UMAP_2 = umap_mono[, 2],
  score  = mono_score,
  stringsAsFactors = FALSE
)
m <- max(abs(quantile(mono_df$score, c(0.02, 0.98), na.rm = TRUE)))
lims <- c(-m, m)
n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(mono_df), big.mark = ","))
)
xr <- quantile(mono_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(mono_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)
p3_left <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.2, vjust = 1.5,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradientn(
    colors = c("#FEE0D2", "#FC9272", "#DE2D26", "#99000D"),
    limits = lims, oob = scales::squish,
    name   = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.2, 0.4),
    ylim = yr + c(-0.1, 1)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD4\u207a T activated:\neffector-like\u2013biosynthetic axis",
    subtitle = "BMMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    # plot.subtitle = element_text(size = 7, color = "grey40", face = "bold"),
    legend.title  = element_text(size = 8),
    legend.text   = element_text(size = 8),
    axis.line     = element_line(linewidth = 0.3),
    plot.margin   = margin(3, 3, 3, 3)
  )
p3_left_t <- p3_left
p3_left_t

pathway_df_4T <- bind_rows(
  hallmark_concordant %>%
    filter(ID %in% c(
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_COMPLEMENT",
      "HALLMARK_MYC_TARGETS_V1"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "HALLMARK_INTERFERON_GAMMA_RESPONSE" ~ "IFN-\u03b3 response",
        ID == "HALLMARK_COMPLEMENT" ~ "Complement",
        ID == "HALLMARK_MYC_TARGETS_V1" ~ "MYC targets"
      ),
      NES     = NES_te,
      padj    = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    ),
  
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response"
      ),
      NES     = NES_te,
      padj    = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    )
) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_4T$NES, na.rm = TRUE)
xmax <- max(pathway_df_4T$NES, na.rm = TRUE)

p3_right_4T <- ggplot(pathway_df_4T, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(label = fdr_lab, x = NES / 2),
    hjust = 0.5,
    size = 2.5,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(xmin + 0.12, xmax + 0.01),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title    = "Reproducible CD4\u207a T effector-like\nand biosynthetic pathways",
    subtitle = "Concordant in train and test; bars show test NES",
    x        = "NES",
    y        = NULL
  ) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_4T


tag_theme <- theme(plot.tag = element_text(face="bold"))


pA <- pA + labs(tag="(a)") + tag_theme
pB <- pB + labs(tag="(b)") + tag_theme
pC <- pC + labs(tag="(c)") + tag_theme
p_ctpnet  <- p_ctpnet + labs(tag="(d)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(e)") + tag_theme
pD <- pD + labs(tag="(f)") + tag_theme
p1 <- p1 + labs(tag="(g)") + tag_theme
p2 <- p2 + labs(tag="(h)") + tag_theme
p3_left_14 <- p3_left_14 + labs(tag="(i)") + tag_theme
p3_right_B <- p3_right_B + labs(tag="(j)") + tag_theme
p3_left_t <- p3_left_t + labs(tag="(k)") + tag_theme
p3_right_4T <- p3_right_4T + labs(tag="(l)") + tag_theme
p4 <- p4 + labs(tag="(m)") + tag_theme
p4 <- p4 +
  facet_wrap(~split, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    legend.position = "bottom"
  )

fig <- 
  
  (pA | pB | pC | p_ctpnet)/
  (p_sclinear|pD |p1 | p2) /
  (p3_left_14 |p3_right_B | p3_left_t |p3_right_4T) /
  p4 + 
  plot_layout(heights = c(1.5, 1.5, 1.5, 1.5)) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    legend.position = "right"
  ) 

fig

ggsave(file.path(out_dir, "Supplementary_Figure_4.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 15.2, height = 11, units = "in")





