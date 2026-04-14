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

# ----------------------------
# Helpers
# ----------------------------
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

# ----------------------------
# Inputs (edit here only)
# ----------------------------

prot <- "anti-human-CD18-totalC"
rna_gene <- "ITGB2"

# prot <- "anti-human-CD7-totalC"
# rna_gene <- "CD7"

base <- file.path(cfg$out_root, "COVID_PBMC_CITE_seq", "LODO_1_test_TP7")
# base <- file.path(cfg$out_root, "COVID_PBMC_CITE_seq", "LODO_1_test_TS5")
dir.create(base, recursive=TRUE, showWarnings=FALSE)
ds    <- "citeseq_v1"
out_dir <- file.path(base, paste0(ds, "_final"))

# ----------------------------
# Load data
# ----------------------------
seu_prep    <- readRDS(file.path(base,  "seu_final.rds"))
train_cells <- readRDS(file.path(base,  "train_cells.rds"))
test_cells  <- readRDS(file.path(base,  "test_cells.rds"))
# groups      <- readRDS(file.path(base, ds, "groups.rds"))
X           <- readRDS(file.path(base,  "X.rds"))
Z           <- readRDS(file.path(base,  "Z.rds"))

X_train <- X[train_cells, , drop=FALSE]
Z_train <- Z[train_cells, , drop=FALSE]
X_test  <- X[test_cells,  , drop=FALSE]
Z_test  <- Z[test_cells,  , drop=FALSE]


stopifnot(prot %in% rownames(seu_prep[["ADT"]]))
stopifnot(rna_gene %in% rownames(seu_prep[["RNA"]]))

# ----------------------------
# Load model object
# ----------------------------
obj   <- readRDS(file.path(base, paste0("final_", prot, ".rds")))
stopifnot(!is.null(obj$fit))
stopifnot(!is.null(obj$alpha_test), !is.null(obj$y_test), !is.null(obj$yhat_test))

alpha_test <- obj$alpha_test
cells <- rownames(alpha_test)

# Diagnostics 
stopifnot(nrow(alpha_test) == length(obj$y_test),
          nrow(alpha_test) == length(obj$yhat_test))
stopifnot(setequal(cells, test_cells))

y    <- as.numeric(obj$y_test);    names(y)    <- cells
yhat <- as.numeric(obj$yhat_test); names(yhat) <- cells

# ----------------------------
# (A–D) Test donor: coupling + diagnostics + prediction
# ----------------------------
hard <- max.col(alpha_test, ties.method="first")
hard <- factor(hard)

tab_h <- table(hard)
print(tab_h)
# ------------------------------------------------------------------------------------

adt_mat <- GetAssayData(seu_prep, assay = "ADT", layer = "data")
stopifnot(prot %in% rownames(adt_mat))

adt_raw <- as.numeric(adt_mat[prot, cells])
names(adt_raw) <- cells
rna_mat <- GetAssayData(seu_prep, assay = "RNA", layer = "data")
stopifnot(rna_gene %in% rownames(rna_mat))
rna_raw <- as.numeric(rna_mat[rna_gene, cells])
names(rna_raw) <- cells

# ------------------------------------------------------------------------------------

hard_lab <- factor(hard, levels=levels(hard),
                   labels=paste0(levels(hard), " (n=", as.integer(tab_h[levels(hard)]), ")"))

rna_raw <- FetchData(seu_prep, vars=rna_gene, cells=cells)[,1]
df <- data.frame(
  cell = cells,
  rna  = rna_raw,
  y    = y[cells],
  yhat = yhat[cells],
  hard = hard_lab,
  stringsAsFactors = FALSE
)
df <- df[is.finite(df$rna) & is.finite(df$y) & is.finite(df$yhat), ]
df$hard <- factor(df$hard)
ct <- seu_prep$predicted.celltype.l1
names(ct) <- colnames(seu_prep)
df$celltype <- ct[df$cell]
stopifnot(!anyNA(df$celltype))
df$celltype <- factor(df$celltype)
fit_ct <- lm(y ~ rna * celltype, data = df)
summary(fit_ct)
# ============================================
# # Add sclienar and ctpnet test predict
# ============================================

res_ctpnet <- readr::read_csv(file.path(base, "ctpnet_test_metrics_FMLEscale.csv"))
pred_ctpnet <- readr::read_csv(file.path(base, "ctpnet_test_predictions_FMLEscale.csv"))
res_sclinear  <- readr::read_csv(file.path(base, "scLinear_test_metrics_fmle_fair.csv"))
pred_sclinear <- readr::read_csv(file.path(base, "scLinear_test_predictions_fmle_fair.csv"))

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

# ============================================
# pA: regime-dependent coupling (test donor)
anova(lm(y ~ hard, data = df))
anova(lm(y ~ rna * hard, data = df))
fit_int <- lm(y ~ rna * hard, data = df)
p_int <- anova(fit_int)[["Pr(>F)"]][match("rna:hard", rownames(anova(fit_int)))]
p_txt <- paste0("Interaction p=", format.pval(p_int, digits = 2))

p_txt <- if (p_int < 2.2e-16) {
  "Interaction p < 2.2 × 10^-16"
} else {
  paste0("Interaction p = ", format(signif(p_int, 3), scientific = TRUE))
}

pA <- ggplot(df, aes(rna, y, color = hard)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  # theme_classic(base_size = 12) +
  theme_classic(base_size = 8) +
  labs(
    x = paste0(rna_gene, " RNA (log1p)"),
    # y = paste0(prot, " protein"),
    y = paste0("CD7", " protein"),
    color = "regime",
    title = "Test donor TS5: \nregime-dependent RNA–protein coupling"
  ) +
  annotate("text", x = Inf, y = Inf, label = p_txt,
           hjust = 1.2, vjust = 1.2, size = 3.6) +
  coord_cartesian(xlim = c(0, 3.8), expand = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# pB: residuals under global mapping (test donor)
m_global <- lm(y ~ rna, data=df)
df$resid_global <- df$y - predict(m_global, newdata=df)

pB <- ggplot(df, aes(rna, resid_global, color=hard)) +
  geom_point(alpha=0.45, size=0.7) +
  geom_hline(yintercept=0, linetype=2, linewidth=0.7) +
  theme_classic(base_size=12) +
  labs(
    x = paste0(rna_gene, " RNA (log1p)"),
    y = "Global model residual",
    title = "Test donor NS1: \nresiduals from a global mapping"
  ) +
  guides(color="none") +
  coord_cartesian(xlim = c(0, 3.8), expand = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# pC/pD: prediction vs truth (test donor cells)
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
    labs(title = m, x = "True protein (test donor)", y = "Prediction") +
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
# ----------------------------
# Build UMAP df for train/test only (aligned set)
# ----------------------------
umap_all <- Embeddings(seu_prep, "ref.umap")   # in case of cross donor
stopifnot(!is.null(rownames(umap_all)))

cells_keep <- intersect(rownames(umap_all), c(train_cells, test_cells))
stopifnot(length(cells_keep) > 0)

umap_use <- umap_all[cells_keep, , drop=FALSE]

donor_split <- rep(NA_character_, length(cells_keep))
names(donor_split) <- cells_keep
donor_split[names(donor_split) %in% train_cells] <- "train_donors"
donor_split[names(donor_split) %in% test_cells]  <- "test_donors"
stopifnot(!anyNA(donor_split))

umap_df <- data.frame(
  cell     = cells_keep,
  UMAP_1   = umap_use[,1],
  UMAP_2   = umap_use[,2],
  donor    = donor_split,
  # celltype = seu_prep$cell_type[cells_keep],  # for kaggle
  celltype = seu_prep$predicted.celltype.l1[cells_keep],
  stringsAsFactors = FALSE
)

# ----------------------------
# (p1–p4) Cross-donor generalization panels (UMAP + entropy + compositions)
# ----------------------------
fit <- obj$fit
pred_tr <- FMLE::fmle_predict(fit, X_new = X_train, Z_new = Z_train, return_se = TRUE)
pred_te <- FMLE::fmle_predict(fit, X_new = X_test,  Z_new = Z_test,  return_se = TRUE)


alpha_te <- pred_te$alpha
alpha_tr <- pred_tr$alpha
if (is.null(rownames(alpha_te))) rownames(alpha_te) <- rownames(X_test)
if (is.null(rownames(alpha_tr))) rownames(alpha_tr) <- rownames(X_train)

stopifnot(setequal(rownames(alpha_te), test_cells))
stopifnot(setequal(rownames(alpha_tr), train_cells))

# p1/p2: test donor UMAP colored by regime + entropy

reg_te <- max.col(alpha_te, ties.method="first")
Hn_te  <- entropy_norm_rows(alpha_te)

tapply(Hn_te, hard, summary)
tapply(Hn_te, seu_prep$predicted.celltype.l1[cells], summary)

Hn_tr  <- entropy_norm_rows(alpha_tr)

names(Hn_tr) <- rownames(alpha_tr)

hard_tr <- max.col(alpha_tr, ties.method="first")
names(hard_tr) <- rownames(alpha_tr)

ct_tr <- seu_prep@meta.data[rownames(alpha_tr), "predicted.celltype.l1", drop=TRUE]

tapply(Hn_tr[names(hard_tr)], factor(hard_tr), summary)
tapply(Hn_tr[rownames(alpha_tr)], factor(ct_tr), summary)

summary(Hn_tr)
summary(Hn_te)

# =======================================
# Train vs test entropy

entropy_shift_stats(Hn_tr, Hn_te)

pred_df_te <- data.frame(
  cell    = rownames(alpha_te),
  regime  = factor(reg_te),
  entropy = as.numeric(Hn_te),
  stringsAsFactors = FALSE
)

plot_df_te <- umap_df %>%
  filter(donor == "test_donors") %>%
  left_join(pred_df_te, by="cell")

stopifnot(nrow(plot_df_te) == length(test_cells))
stopifnot(!anyNA(plot_df_te$regime), !anyNA(plot_df_te$entropy))

p1 <- ggplot(plot_df_te, aes(UMAP_1, UMAP_2, color = regime)) +
  geom_point(size = 0.25, alpha = 0.7) +
  theme_classic(base_size = 8) +
  labs(title = "Test donor TS5 UMAP: regime (argmax α)", color = "Regime") +        # NS1, TP7 and TS5
  theme(plot.title = element_text(hjust = 0.5,  face = "bold"))

p2 <- ggplot(plot_df_te, aes(UMAP_1, UMAP_2, color = entropy)) +
  geom_point(size = 0.25, alpha = 0.7) +
  theme_classic(base_size = 8) +
  labs(title = "Test donor TP7 UMAP: \nα-entropy", color = "Entropy")

# p3: regime composition train vs test + JSD
reg_tr <- max.col(alpha_tr, ties.method="first")

pred_df_tr <- data.frame(
  cell   = rownames(alpha_tr),
  regime = factor(reg_tr),
  stringsAsFactors = FALSE
)

comp_cells <- umap_df %>%
  filter(donor %in% c("train_donors","test_donors")) %>%
  select(cell, donor) %>%
  left_join(bind_rows(pred_df_tr, pred_df_te %>% select(cell, regime)), by="cell")

stopifnot(!anyNA(comp_cells$regime))

comp <- comp_cells %>%
  count(donor, regime, name="n") %>%
  group_by(donor) %>%
  mutate(frac = n/sum(n)) %>%
  ungroup()

p3 <- ggplot(comp, aes(donor, frac, fill = regime)) +
  geom_col(width=0.8) +
  theme_classic(base_size = 8) +
  labs(title="Regime composition: train donors vs test donors", x=NULL, y="Fraction")+
  coord_cartesian(xlim = c(0, 3.8), expand = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold"))

# JSD on matched regime support
regs <- sort(unique(comp$regime))
pvec <- comp %>% filter(donor=="train_donors") %>%
  right_join(data.frame(regime=regs), by="regime") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac)) %>% arrange(regime) %>% pull(frac)
qvec <- comp %>% filter(donor=="test_donors") %>%
  right_join(data.frame(regime=regs), by="regime") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac)) %>% arrange(regime) %>% pull(frac)

JSD <- jsd(pvec, qvec, base=2)
cat(sprintf("JSD(train,test)=%.6g\n", JSD))

# p4: regime within celltypes, faceted train vs test
min_n <- 100
# ct_col <- "cell_type" # for kaggle datasets
ct_col <- "predicted.celltype.l1"   # for real cross donor datasets
comp_tr_ct <- celltype_regime_comp(alpha_tr, seu_prep, cell_ids=train_cells,
                                   celltype_col=ct_col, min_n=min_n, split_label="train donors")
comp_te_ct <- celltype_regime_comp(alpha_te, seu_prep, cell_ids=test_cells,
                                   celltype_col=ct_col, min_n=min_n, split_label="test donors TS5")
comp_all_ct <- bind_rows(comp_tr_ct, comp_te_ct)

p4 <- ggplot(comp_all_ct, aes(x = celltype, y = frac, fill = regime)) +
  geom_col(width = 0.85) +
  facet_wrap(~split, nrow=1) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  labs(title="Regime composition within cell types", x=NULL, y="Fraction")+
  theme(plot.title = element_text(hjust = 0.5,  face = "bold"))

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


# ==============================================================================
#  CD18 donor TP7
# =============================================================================
lineage_use <- "Mono"
ct_col <- "predicted.celltype.l1" 

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
umap_mono <- Embeddings(seu_prep, "ref.umap")[names(mono_score), , drop = FALSE]
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
print(q)
print(m)
summary(mono_df$score)

library(FNN)

# mono_df already has UMAP_1, UMAP_2, score
xy <- as.matrix(mono_df[, c("UMAP_1", "UMAP_2")])

kn <- FNN::get.knn(xy, k = 25)$nn.index
score_sm <- sapply(seq_len(nrow(xy)), function(i) {
  median(mono_df$score[kn[i, ]], na.rm = TRUE)
})

mono_df$score_sm <- score_sm
q <- quantile(mono_df$score_sm, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

p3_left_mono <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score_sm)) +
  geom_point(size = 1.0, alpha = 0.95) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1,
    vjust       = 2.0,
    size        = 2.8,
    color       = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = lims, oob = scales::squish,
    name = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "Monocytes:IFN-activated versus \nresting/biosynthetic axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  )+
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_mono

pathway_df_mono <- bind_rows(
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" ~ "Cytokine signaling",
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
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

p3_right_mono <- ggplot(pathway_df_mono, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    # Monocytes: positive NES = biosynthetic/resting = BLUE
    #            negative NES = IFN-activated        = RED
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x = ifelse(NES > 0, NES / 2, NES / 2)
    ),
    hjust    = 0.5,
    size     = 2.0,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(min(0, xmin) + 0.14, xmax + 0.15),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    x        = "Test-donor TP7",
    y        = NULL,
    title = "Reproducible monocyte IFN-activated\nand biosynthetic pathways",
    subtitle = "Concordant in train and test donors;\nbars show test-donor TP7"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )
p3_right_mono




lineage_use <- "CD4 T"
ct_col <- "predicted.celltype.l1" 
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

cd4_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_cd4 <- as.matrix(alpha_te[cd4_cells, , drop = FALSE])
hard_cd4 <- colnames(A_cd4)[max.col(A_cd4, ties.method = "first")]
tab_cd4  <- sort(table(hard_cd4), decreasing = TRUE)

score_regimes <- names(tab_cd4)[1:2]

A2 <- A_cd4[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

cd4_score <- A2[, 1] - A2[, 2]
names(cd4_score) <- rownames(A2)

umap_cd4 <- Embeddings(seu_prep, "ref.umap")[names(cd4_score), , drop = FALSE]

cd4_df <- data.frame(
  cell   = names(cd4_score),
  UMAP_1 = umap_cd4[, 1],
  UMAP_2 = umap_cd4[, 2],
  score  = cd4_score,
  stringsAsFactors = FALSE
)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(cd4_df), big.mark = ","))
)

q <- quantile(cd4_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

p3_left_4T <- ggplot(cd4_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 5.0,
    size        = 2.8, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD4\u207a T cells:\nactivation-associated\u2013biosynthetic axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_4T

pathway_df_4T <- bind_rows(
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM",
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_RNA",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" ~ "Cytokine signaling",
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
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
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x     = NES * 0.55
    ),
    hjust    = 0.5,
    size     = 2.0,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(0, xmax + 0.15),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x        = "Test-donor TS5",
    y        = NULL,
    title = "Reproducible CD4\u207a T activated\nand biosynthetic pathways",
    subtitle = "Concordant in train donor and test donor; \nbars show test donor TS5"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_4T


lineage_use <- "CD8 T"

ct_col <- "predicted.celltype.l1" 
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

cd4_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_cd4 <- as.matrix(alpha_te[cd4_cells, , drop = FALSE])
hard_cd4 <- colnames(A_cd4)[max.col(A_cd4, ties.method = "first")]
tab_cd4  <- sort(table(hard_cd4), decreasing = TRUE)

score_regimes <- names(tab_cd4)[1:2]

A2 <- A_cd4[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

cd4_score <- A2[, 1] - A2[, 2]
names(cd4_score) <- rownames(A2)

umap_cd4 <- Embeddings(seu_prep, "ref.umap")[names(cd4_score), , drop = FALSE]

cd4_df <- data.frame(
  cell   = names(cd4_score),
  UMAP_1 = umap_cd4[, 1],
  UMAP_2 = umap_cd4[, 2],
  score  = cd4_score,
  stringsAsFactors = FALSE
)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(cd4_df), big.mark = ","))
)

q <- quantile(cd4_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

p3_left_8T <- ggplot(cd4_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 5.0,
    size        = 2.8, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD8\u207a T cells:\ncytotoxic\u2013biosynthetic axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_8T

pathway_df_8T <- bind_rows(
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TRANSLATION",
      "REACTOME_METABOLISM_OF_RNA",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
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

p3_right_8T <- ggplot(pathway_df_8T, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x     = NES * 0.55
    ),
    hjust    = 0.5,
    size     = 2.0,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(min(0, xmin) - 1.55, xmax + 0.15),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    x        = "Test-donor TS5",
    y        = NULL,
    title = "Reproducible CD8\u207a T biosynthetic pathways",
    subtitle = "Concordant in train donor and test donor; \nbars show test donor TS5"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_8T

tag_theme <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.98)
)

tag_theme1 <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.9)
)


pA         <- pA         + labs(tag="(a)") + tag_theme
pC         <- pC         + labs(tag="(b)") + tag_theme
p_ctpnet   <- p_ctpnet   + labs(tag="(c)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(d)") + tag_theme
pD         <- pD         + labs(tag="(e)") + tag_theme
p1         <- p1         + labs(tag="(f)") + tag_theme1
p4          <- p4          + labs(tag="(g)") + tag_theme1
p3_left_4T  <- p3_left_4T  + labs(tag="(h)") + tag_theme1
p3_right_4T <- p3_right_4T + labs(tag="(i)") + tag_theme1
p3_left_8T  <- p3_left_8T  + labs(tag="(j)") + tag_theme1
p3_right_8T <- p3_right_8T + labs(tag="(k)") + tag_theme1
p3_left_mono  <- p3_left_mono  + labs(tag="(l)") + tag_theme1
p3_right_mono <- p3_right_mono + labs(tag="(m)") + tag_theme1


fig <-
  (pA | pC | p_ctpnet | p_sclinear | pD ) /
  (p1|p4| p3_left_4T | p3_right_4T )/ (p3_left_8T | p3_right_8T| p3_left_mono | p3_right_mono) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    # legend.position = "right",
    # legend.box = "vertical"
  )

fig
dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Supplementary_Figure_6.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 13.6, height = 8, units = "in")
# ========================================================================
# CD18 donor  TS5
# ===================================================================
lineage_use <- "Mono"
ct_col <- "predicted.celltype.l1" 

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
umap_mono <- Embeddings(seu_prep, "ref.umap")[names(mono_score), , drop = FALSE]
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
print(q)
print(m)
summary(mono_df$score)

library(FNN)

# mono_df already has UMAP_1, UMAP_2, score
xy <- as.matrix(mono_df[, c("UMAP_1", "UMAP_2")])

kn <- FNN::get.knn(xy, k = 25)$nn.index
score_sm <- sapply(seq_len(nrow(xy)), function(i) {
  median(mono_df$score[kn[i, ]], na.rm = TRUE)
})

mono_df$score_sm <- score_sm
q <- quantile(mono_df$score_sm, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

p3_left_mono <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score_sm)) +
  geom_point(size = 1.0, alpha = 0.95) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1,
    vjust       = 2.0,
    size        = 2.8,
    color       = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = lims, oob = scales::squish,
    name = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "Monocytes:\nIFN-associated–biosynthetic axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1", y = "UMAP 2"
  )+
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_mono

pathway_df_mono <- bind_rows(
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_RNA",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
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

p3_right_mono <- ggplot(pathway_df_mono, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    # Monocytes: positive NES = biosynthetic/resting = BLUE
    #            negative NES = IFN-activated        = RED
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x     = NES + 0.1
    ),
    hjust    = 0.0,
    size     = 2.0,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(xmin + 0.05, 0),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    x        = "Test-donor TS5",
    y        = NULL,
    title    = "Reproducible monocyte \nbiosynthetic pathways",
    subtitle = "Concordant in train and test donors;\nbars show test-donor TS5"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )
p3_right_mono




lineage_use <- "CD4 T"

ct_col <- "predicted.celltype.l1" 

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

cd4_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_cd4 <- as.matrix(alpha_te[cd4_cells, , drop = FALSE])
hard_cd4 <- colnames(A_cd4)[max.col(A_cd4, ties.method = "first")]
tab_cd4  <- sort(table(hard_cd4), decreasing = TRUE)

score_regimes <- names(tab_cd4)[1:2]

A2 <- A_cd4[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

cd4_score <- A2[, 1] - A2[, 2]
names(cd4_score) <- rownames(A2)

umap_cd4 <- Embeddings(seu_prep, "ref.umap")[names(cd4_score), , drop = FALSE]

cd4_df <- data.frame(
  cell   = names(cd4_score),
  UMAP_1 = umap_cd4[, 1],
  UMAP_2 = umap_cd4[, 2],
  score  = cd4_score,
  stringsAsFactors = FALSE
)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(cd4_df), big.mark = ","))
)

q <- quantile(cd4_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

p3_left_4T <- ggplot(cd4_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 5.0,
    size        = 2.8, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD4\u207a T cells:\nactivation-associated\u2013biosynthetic axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_4T

pathway_df_4T <- bind_rows(
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_RNA",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
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
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x     = NES * 0.55
    ),
    hjust    = 0.5,
    size     = 2.0,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(0, xmax + 0.15),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  scale_x_continuous(expand = c(0, 0))+
  labs(
    x        = "Test-donor TS5",
    y        = NULL,
    title = "Reproducible CD4\u207a T biosynthetic \npathways",
    subtitle = "Concordant in train donor and test donor; \nbars show test donor TS5"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_4T


lineage_use <- "CD8 T"

ct_col <- "predicted.celltype.l1" 

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

cd4_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_cd4 <- as.matrix(alpha_te[cd4_cells, , drop = FALSE])
hard_cd4 <- colnames(A_cd4)[max.col(A_cd4, ties.method = "first")]
tab_cd4  <- sort(table(hard_cd4), decreasing = TRUE)

score_regimes <- names(tab_cd4)[1:2]

A2 <- A_cd4[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

cd4_score <- A2[, 1] - A2[, 2]
names(cd4_score) <- rownames(A2)

umap_cd4 <- Embeddings(seu_prep, "ref.umap")[names(cd4_score), , drop = FALSE]

cd4_df <- data.frame(
  cell   = names(cd4_score),
  UMAP_1 = umap_cd4[, 1],
  UMAP_2 = umap_cd4[, 2],
  score  = cd4_score,
  stringsAsFactors = FALSE
)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(cd4_df), big.mark = ","))
)

q <- quantile(cd4_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

p3_left_8T <- ggplot(cd4_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 5.0,
    size        = 2.8, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD8\u207a T cells:\ncytotoxic\u2013biosynthetic axis",
    subtitle = "PBMC, test cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_8T

pathway_df_8T <- bind_rows(
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_HEMOSTASIS",
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_RNA",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_HEMOSTASIS" ~ "Hemostasis",
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
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

p3_right_8T <- ggplot(pathway_df_8T, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x     = NES * 0.55
    ),
    hjust    = 0.5,
    size     = 2.0,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(min(0, xmin) - 1.55, xmax + 0.15),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    x        = "Test-donor TS5",
    y        = NULL,
    title = "Reproducible CD8\u207a T biosynthetic \npathways",
    subtitle = "Concordant in train donor and test donor; \nbars show test donor TS5"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_8T

tag_theme <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.98)
)

tag_theme1 <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.9)
)


pA         <- pA         + labs(tag="(a)") + tag_theme
pC         <- pC         + labs(tag="(b)") + tag_theme
p_ctpnet   <- p_ctpnet   + labs(tag="(c)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(d)") + tag_theme
pD         <- pD         + labs(tag="(e)") + tag_theme
p1         <- p1         + labs(tag="(f)") + tag_theme1
p4          <- p4          + labs(tag="(g)") + tag_theme1
p3_left_4T  <- p3_left_4T  + labs(tag="(h)") + tag_theme1
p3_right_4T <- p3_right_4T + labs(tag="(i)") + tag_theme1
p3_left_8T  <- p3_left_8T  + labs(tag="(j)") + tag_theme1
p3_right_8T <- p3_right_8T + labs(tag="(k)") + tag_theme1
p3_left_mono  <- p3_left_mono  + labs(tag="(l)") + tag_theme1
p3_right_mono <- p3_right_mono + labs(tag="(m)") + tag_theme1


fig <-
  (pA | pC | p_ctpnet | p_sclinear | pD ) /
  (p1|p4| p3_left_4T | p3_right_4T )/ (p3_left_8T | p3_right_8T| p3_left_mono | p3_right_mono) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    # legend.position = "right",
    # legend.box = "vertical"
  )

fig
dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Supplementary_Figure_8.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 13.6, height = 8, units = "in")



# ==============================================================================
# CD7 # Donor TP7
# ==============================================================================
lineage_use <- "CD4 T"
ct_col <- "predicted.celltype.l1" 

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

cd4_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_cd4 <- as.matrix(alpha_te[cd4_cells, , drop = FALSE])
hard_cd4 <- colnames(A_cd4)[max.col(A_cd4, ties.method = "first")]
tab_cd4  <- sort(table(hard_cd4), decreasing = TRUE)

score_regimes <- names(tab_cd4)[1:2]
A2 <- A_cd4[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

cd4_score <- A2[, 1] - A2[, 2]
names(cd4_score) <- rownames(A2)

umap_cd4 <- Embeddings(seu_prep, "ref.umap")[names(cd4_score), , drop = FALSE]

cd4_df <- data.frame(
  cell   = names(cd4_score),
  UMAP_1 = umap_cd4[, 1],
  UMAP_2 = umap_cd4[, 2],
  score  = cd4_score,
  stringsAsFactors = FALSE
)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(cd4_df), big.mark = ","))
)

q <- quantile(cd4_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

xr <- quantile(cd4_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(cd4_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)

p3_left_4T <- ggplot(cd4_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 2.0,
    vjust       = 2.0,
    size        = 3.5,
    color       = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.5, 1),
    ylim = yr + c(-0.5, 3)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD4\u207a T cells:\ncytokine-associated\u2013biosynthetic axis",
    subtitle = "PBMC, test donor cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_4T


# ----------------------------
pathway_df_4T <- reactome_concordant %>%
  filter(ID %in% c(
    "REACTOME_TRANSLATION",
    "REACTOME_METABOLISM_OF_RNA",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
    "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "REACTOME_TRANSLATION" ~ "Translation",
      ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
      ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
      ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
      ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
      ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD",
      ID == "REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM" ~ "Cytokine signaling"
    ),
    NES     = NES_te,
    padj    = padj_te,
    fdr_lab = sprintf("FDR %.1e", padj_te)
  ) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_4T$NES, na.rm = TRUE)
xmax <- max(pathway_df_4T$NES, na.rm = TRUE)

p3_right_4T <- ggplot(pathway_df_4T, aes(x = NES, y = pathway)) +
  geom_col(
    aes(fill = NES > 0),
    width       = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      label = fdr_lab,
      x     = NES * 0.55
    ),
    hjust    = 0.5,
    size     = 2.8,
    color    = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(
    xlim = c(xmin + 0.1, xmax + 0.15),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    x        = "Test-donor TP7",
    y        = NULL,
    title = "Reproducible CD4\u207a T cytokine-\nsignaling and biosynthetic pathways",
    subtitle = "Concordant in train and test donors;\nbars show test-donor TP7"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_4T


lineage_use <- "CD8 T"

ct_col <- "predicted.celltype.l1" 

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
umap_mono <- Embeddings(seu_prep, "ref.umap")[names(mono_score), , drop = FALSE]
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
p3_left_8T <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 1.1, vjust = 2.0,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.5, 0.5),
    ylim = yr + c(-0.5, 0.5)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD8\u207a T cells:\neffector-associated\u2013biosynthetic axis",
    subtitle = "PBMC, test donor cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )
p3_left_8T

# ----------------------------
pathway_df_8T <- reactome_concordant %>%
  filter(ID %in% c(
    "REACTOME_TRANSLATION",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_METABOLISM_OF_RNA",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "REACTOME_TRANSLATION" ~ "Translation",
      ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
      ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
      ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
      ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
      ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
    ),
    NES     = NES_te,
    padj    = padj_te,
    fdr_lab = sprintf("FDR %.1e", padj_te)
  ) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_8T$NES, na.rm = TRUE)
xmax <- max(pathway_df_8T$NES, na.rm = TRUE)
left_vis <- xmin - 0.23
pathway_df_8T <- pathway_df_8T %>%
  mutate(
    fdr_x = left_vis + 0.5 * (NES - left_vis)
  )

p3_right_8T <- ggplot(pathway_df_8T) +
  geom_rect(
    aes(
      xmin = left_vis,
      xmax = NES,
      ymin = as.numeric(pathway) - 0.31,
      ymax = as.numeric(pathway) + 0.31,
      fill = NES > 0
    ),
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      x = fdr_x,
      y = pathway,
      label = fdr_lab
    ),
    hjust = 0.5,
    size = 2.6,
    color = "white",
    fontface = "bold"
  ) +
  coord_cartesian(
    xlim = c(-2.540, -1.9),
    clip = "off"
  ) +
  theme_classic(base_size = 8) +
  labs(
    x = "Test-donor TP7",
    y = NULL,
    title = "Reproducible CD8\u207a T biosynthetic \npathways",
    subtitle = "Concordant in train and test donors; \nbars show test-donor TP7"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_8T


library(patchwork)
library(ggplot2)

tag_theme <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.98)
)

pA         <- pA         + labs(tag="(a)") + tag_theme
pC         <- pC         + labs(tag="(b)") + tag_theme
p_ctpnet   <- p_ctpnet   + labs(tag="(c)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(d)") + tag_theme
pD         <- pD         + labs(tag="(e)") + tag_theme
p1         <- p1         + labs(tag="(f)") + tag_theme
p4          <- p4          + labs(tag="(g)") + tag_theme
p3_left_4T  <- p3_left_4T  + labs(tag="(h)") + tag_theme
p3_right_4T <- p3_right_4T + labs(tag="(i)") + tag_theme
p3_left_8T  <- p3_left_8T  + labs(tag="(j)") + tag_theme
p3_right_8T <- p3_right_8T + labs(tag="(k)") + tag_theme


fig <-
  (pA | pC | p_ctpnet | p_sclinear) /
  (pD | p1| p4)/ (p3_left_4T | p3_right_4T|p3_left_8T | p3_right_8T)  + 
  plot_layout(heights = c(1, 1, 1)) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    # legend.position = "right",
    # legend.box = "vertical"
  )

fig
dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Supplementary_Figure_9.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 14, height = 8, units = "in")




# =================
# donor 3 TS5
# =================
lineage_use <- "CD4 T"
ct_col <- "predicted.celltype.l1" 

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

cd4_cells <- intersect(
  test_cells,
  colnames(seu_prep)[meta[[ct_col]] == lineage_use]
)

A_cd4 <- as.matrix(alpha_te[cd4_cells, , drop = FALSE])
hard_cd4 <- colnames(A_cd4)[max.col(A_cd4, ties.method = "first")]
tab_cd4  <- sort(table(hard_cd4), decreasing = TRUE)

score_regimes <- names(tab_cd4)[1:2]
A2 <- A_cd4[, score_regimes, drop = FALSE]
A2 <- A2[rowSums(A2) > 0, , drop = FALSE]
A2 <- A2 / rowSums(A2)

cd4_score <- A2[, 1] - A2[, 2]
names(cd4_score) <- rownames(A2)

umap_cd4 <- Embeddings(seu_prep, "ref.umap")[names(cd4_score), , drop = FALSE]

cd4_df <- data.frame(
  cell   = names(cd4_score),
  UMAP_1 = umap_cd4[, 1],
  UMAP_2 = umap_cd4[, 2],
  score  = cd4_score,
  stringsAsFactors = FALSE
)

n_label <- data.frame(
  x     = Inf,
  y     = Inf,
  label = paste0("n = ", format(nrow(cd4_df), big.mark = ","))
)

q <- quantile(cd4_df$score, c(0.02, 0.98), na.rm = TRUE)
m <- max(abs(q))
lims <- c(-m, m)

xr <- quantile(cd4_df$UMAP_1, c(0.01, 0.99), na.rm = TRUE)
yr <- quantile(cd4_df$UMAP_2, c(0.01, 0.99), na.rm = TRUE)

p3_left_4T <- ggplot(cd4_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 2.0,
    vjust       = 2.0,
    size        = 3.5,
    color       = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.5, 1),
    ylim = yr + c(-0.5, 5)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD4\u207a T cells:\nactivation-associated\u2013biosynthetic axis",
    subtitle = "PBMC, test donor cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_4T

pathway_df_4T <- reactome_concordant %>%
  filter(ID %in% c(
    "REACTOME_TRANSLATION",
    "REACTOME_METABOLISM_OF_RNA",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "REACTOME_TRANSLATION" ~ "Translation",
      ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
      ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
      ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
      ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
      ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
    ),
    NES     = NES_te,
    padj    = padj_te,
    fdr_lab = sprintf("FDR %.2e", padj_te)
  ) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin_4T <- min(pathway_df_4T$NES, na.rm = TRUE)
xmax_4T <- max(pathway_df_4T$NES, na.rm = TRUE)
left_vis_4T <- xmin_4T - 0.45

pathway_df_4T <- pathway_df_4T %>%
  mutate(
    NES_mid = (left_vis_4T + NES) / 2
  )

p3_right_4T <- ggplot(pathway_df_4T) +
  geom_rect(
    aes(
      xmax = NES,
      ymin = as.numeric(pathway) - 0.31,
      ymax = as.numeric(pathway) + 0.31
    ),
    xmin = left_vis_4T,
    fill = "#4393C3"
  ) +
  geom_text(
    aes(x = NES_mid, y = pathway, label = fdr_lab),
    hjust = 0.5, size = 2.6, color = "white", fontface = "bold"
  ) +
  coord_cartesian(xlim = c(left_vis_4T, xmax_4T + 0.2), clip = "off") +
  theme_classic(base_size = 8) +
  labs(
    x = "Test-donor TS5",
    y = NULL,
    title = "Reproducible CD4\u207a T biosynthetic\npathways",
    subtitle = "Concordant in train and test donors;\nbars show test-donor TS5"
  ) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 35)
  )

p3_right_4T


lineage_use <- "CD8 T"
ct_col <- "predicted.celltype.l1" 

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
umap_mono <- Embeddings(seu_prep, "ref.umap")[names(mono_score), , drop = FALSE]
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
p3_left_8T <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data        = n_label,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust       = 2, vjust = 2.0,
    size        = 3.5, color = "black",
    fontface    = "bold"
  ) +
  scale_color_gradient2(
    low      = "#2166AC",
    mid      = "white",
    high     = "#B2182B",
    midpoint = 0,
    limits   = lims,
    oob      = scales::squish,
    name     = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.5, 0.5),
    ylim = yr + c(-0.5, 2)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "CD8\u207a T cells:\neffector-associated\u2013biosynthetic axis",
    subtitle = "PBMC, test donor cells",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )
p3_left_8T

pathway_df_8T <- reactome_concordant %>%
  filter(ID %in% c(
    "REACTOME_TRANSLATION",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_METABOLISM_OF_RNA",
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "REACTOME_TRANSLATION" ~ "Translation",
      ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
      ID == "REACTOME_METABOLISM_OF_RNA" ~ "RNA metabolism",
      ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
      ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
      ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD"
    ),
    NES     = NES_te,
    padj    = padj_te,
    fdr_lab = sprintf("FDR %.1e", padj_te)
  ) %>%
  arrange(NES) %>%
  mutate(pathway = factor(pathway, levels = pathway))



xmin_8T <- min(pathway_df_8T$NES, na.rm = TRUE)
xmax_8T <- max(pathway_df_8T$NES, na.rm = TRUE)
left_vis_8T <- xmin_8T - 0.2

pathway_df_8T <- pathway_df_8T %>%
  mutate(
    fdr_x = left_vis_8T + 0.5 * (NES - left_vis_8T)
  )

p3_right_8T <- ggplot(pathway_df_8T) +
  geom_rect(
    aes(
      xmax = NES,
      ymin = as.numeric(pathway) - 0.31,
      ymax = as.numeric(pathway) + 0.31,
      fill = NES > 0
    ),
    xmin = left_vis_8T,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  ) +
  geom_text(
    aes(
      x = fdr_x,
      y = pathway,
      label = fdr_lab
    ),
    hjust = 0.5,
    size = 2.6,
    color = "white",
    fontface = "bold"
  ) +
  coord_cartesian(xlim = c(left_vis_8T, xmax_8T + 0.14), clip = "off") +
  theme_classic(base_size = 8) +
  labs(
    x = "Test-donor TS5",
    y = NULL,
    title = "Reproducible CD8\u207a T biosynthetic\npathways",
    subtitle = "Concordant in train and test donors;\nbars show test-donor TS5"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 8),
    plot.subtitle = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_8T


library(patchwork)
library(ggplot2)

tag_theme <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.98)
)

pA         <- pA         + labs(tag="(a)") + tag_theme
pC         <- pC         + labs(tag="(b)") + tag_theme
p_ctpnet   <- p_ctpnet   + labs(tag="(c)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(d)") + tag_theme
pD         <- pD         + labs(tag="(e)") + tag_theme
p1         <- p1         + labs(tag="(f)") + tag_theme
p4          <- p4          + labs(tag="(g)") + tag_theme
p3_left_4T  <- p3_left_4T  + labs(tag="(h)") + tag_theme
p3_right_4T <- p3_right_4T + labs(tag="(i)") + tag_theme
p3_left_8T  <- p3_left_8T  + labs(tag="(j)") + tag_theme
p3_right_8T <- p3_right_8T + labs(tag="(k)") + tag_theme


fig <-
  (pA | pC | p_ctpnet | p_sclinear) /
  (pD | p1| p4)/ (p3_left_4T | p3_right_4T|p3_left_8T | p3_right_8T)  + 
  plot_layout(heights = c(1, 1, 1)) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    # legend.position = "right",
    # legend.box = "vertical"
  )
fig
dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Supplementary_Figure_10.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 13.7, height = 8, units = "in")






