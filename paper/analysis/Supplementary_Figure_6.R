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
# Load data
# ----------------------------
#  BMMC CITE-seq → healthy young BM 

# for BMMC CITE-seq  train
prot <- "CD33"
rna_gene <- "CD33"

base1 <- file.path(cfg$out_root, "Trasnfer_BMMC_CITE_to_healthy_young_BM")
dir.create(base1, recursive = TRUE, showWarnings = FALSE)
ds    <- "citeseq_v1"
seu_prepC     <- readRDS(file.path(base1, ds, "seu_final.rds"))
train_cellsC  <- readRDS(file.path(base1, ds, "train_cells.rds"))
XC           <- readRDS(file.path(base1, ds, "X.rds"))
ZC           <- readRDS(file.path(base1, ds, "Z.rds"))
XC_train  <- XC[train_cellsC,  , drop=FALSE]
ZC_train  <- ZC[train_cellsC,  , drop=FALSE]

# Young healthy BM test
base2 <- file.path(cfg$out_root, "Trasnfer_healthy_young_BM_to_BMMC_CITE")
seu_prepA    <- readRDS(file.path(base2, ds, "seu_final.rds"))
test_cellsA  <- readRDS(file.path(base2, ds, "test_cells.rds"))
XA           <- readRDS(file.path(base2, ds, "X.rds"))
ZA           <- readRDS(file.path(base2, ds, "Z.rds"))
XA_test <- XA[test_cellsA, , drop=FALSE]
ZA_test <- ZA[test_cellsA, , drop=FALSE]
# ------------------------------------------------------------------------------------

stopifnot(prot %in% rownames(seu_prepC[["ADT"]]))
stopifnot(rna_gene %in% rownames(seu_prepC[["RNA"]]))

canon_prot_global <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  
  # remove antibody suffix in both raw styles first
  x <- sub("\\.AB$", "", x)
  x <- sub("-AB$", "", x)
  
  # normalize punctuation
  x <- gsub("\\.", "-", x)
  x <- gsub("--+", "-", x)
  x <- gsub("-+$", "", x)
  
  # remove suffix again after dot->dash conversion
  x <- sub("-AB$", "", x)
  
  alias_map <- c(
    "HLA-ABC"   = "HLA-A-B-C",
    "HLA-A-B-C" = "HLA-A-B-C",
    
    "TCRab" = "TCR",
    "TCR"   = "TCR",
    
    "TCRgd" = "TCRgd",
    
    "Tim3"  = "TIM-3",
    "TIM3"  = "TIM-3",
    "TIM-3" = "TIM-3",
    
    "B7-H4" = "B7-H4",
    "B7H4"  = "B7-H4",
    
    "IL21R"  = "IL-21R",
    "IL-21R" = "IL-21R"
  )
  
  idx <- match(x, names(alias_map))
  hit <- !is.na(idx)
  x[hit] <- unname(alias_map[idx[hit]])
  
  x
}


new <- canon_prot_global(rownames(seu_prepA[["ADT"]]))
stopifnot(!anyDuplicated(new))
m <- GetAssayData(seu_prepA, assay = "ADT", layer = "counts")
rownames(m) <- new
seu_prepA[["ADT"]] <- CreateAssayObject(counts = m)

stopifnot(prot %in% rownames(seu_prepA[["ADT"]]))
stopifnot(rna_gene %in% rownames(seu_prepA[["RNA"]]))

# ----------------------------
# Load model object train
# ----------------------------
# A = healthy young BM and B = BMMC CITE-seq
base1 <- file.path(cfg$out_root, "Trasnfer_BMMC_CITE_to_healthy_young_BM")
out_dir <- file.path(base1, "Transfer_BA")
obj   <- readRDS(file.path(out_dir, paste0("final_", prot, ".rds")))
stopifnot(!is.null(obj$fit))
stopifnot(!is.null(obj$alpha_test), !is.null(obj$y_test))

yhat_vec <- NULL
if (!is.null(obj$yhat_test_cal)) {
  yhat_vec <- obj$yhat_test_cal
} else if (!is.null(obj$yhat_test)) {
  yhat_vec <- obj$yhat_test
} else {
  stop("No yhat_test(_cal) found in obj.")
}
stopifnot(length(yhat_vec) == length(obj$y_test))



alpha_test <- obj$alpha_test
cells <- rownames(alpha_test)

# Diagnostics (exactly what you wanted)
stopifnot(setequal(cells, test_cellsA))
stopifnot(nrow(alpha_test) == length(obj$y_test),
          nrow(alpha_test) == length(yhat_vec))

y    <- as.numeric(obj$y_test);    names(y)    <- cells
yhat <- as.numeric(yhat_vec); names(yhat) <- cells

# ----------------------------
# Test donor: coupling + diagnostics + prediction
# ----------------------------
hard <- max.col(alpha_test, ties.method="first")
hard <- factor(hard)

tab_h <- table(hard)
print(tab_h)

# ------------------------------------------------------------------------------------
adt_matA <- GetAssayData(seu_prepA, assay="ADT", layer="data")
rna_matA <- GetAssayData(seu_prepA, assay="RNA", layer="data")

adt_raw <- as.numeric(adt_matA[prot, cells]); names(adt_raw) <- cells

stopifnot(rna_gene %in% rownames(rna_matA))

rna_raw <- as.numeric(rna_matA[rna_gene, cells]); names(rna_raw) <- cells

hard_lab <- factor(hard, levels=levels(hard),
                   labels=paste0(levels(hard), " (n=", as.integer(tab_h[levels(hard)]), ")"))

ct <- seu_prepA$predicted.celltype.l1

names(ct) <- colnames(seu_prepA)
celltype_vec <- ct[cells]
rna_log <- log1p(as.numeric(rna_raw))

df <- data.frame(
  cell = cells,
  rna  = rna_raw,
  y    = y[cells],
  yhat = yhat[cells],
  hard = hard_lab,
  celltype = as.factor(celltype_vec),  
  stringsAsFactors = FALSE
)
df <- df[is.finite(df$rna) & is.finite(df$y) & is.finite(df$yhat), ]
stopifnot(!anyNA(df$celltype))

H <- -rowSums(alpha_test * log(pmax(alpha_test, 1e-12))) / log(ncol(alpha_test))
entropy_error_stats(H, df$y, df$yhat)

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

# ======================================
# Panel A
# ======================================

pA <- ggplot(df, aes(rna, y, color = hard)) +
  geom_point(alpha = 0.45, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  theme_classic(base_size = 8) +                                 
  labs(
    x = paste0(rna_gene, " RNA (log1p)"),
    y = paste0(prot, " protein"),
    color = "regime",
    title="Test (Young healthy BM):\nregime-dependent RNA–protein coupling") +
  annotate("text", x = Inf, y = Inf, label = p_txt,
           hjust = 1.2, vjust = 1.2, size = 3.6) +
  coord_cartesian(xlim = c(0, 3.8), expand = FALSE) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

out_dir <- file.path(base1, "Transfer_BA")

res_sclinear  <- readr::read_csv(file.path(out_dir, "res_scl_BA.csv"))
pred_sclinear <- readr::read_csv(file.path(out_dir, "scLinear_test_predictions_FMLEscale_B_to_A.csv"))

res_ctpnet <- readr::read_csv(file.path(out_dir, "res_ctp_BA.csv"))
pred_ctpnet <- readr::read_csv(file.path(out_dir, "cTPnet_test_predictions_FMLEscale_B_to_A.csv"))

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


# pB: residuals under global mapping (test donor)
m_global <- lm(y ~ rna, data=df)
df$resid_global <- df$y - predict(m_global, newdata=df)

# ======================================
# Panel B
# ======================================

pB <- ggplot(df, aes(rna, resid_global, color=hard)) +
  geom_point(alpha=0.45, size=0.7) +
  geom_hline(yintercept=0, linetype=2, linewidth=0.7) +
  theme_classic(base_size=8) +
  labs(
    x = paste0(rna_gene, " RNA (log1p)"),
    y = "Global model residual",
    title="Test (Young healthy BM):\nresiduals under a global mapping") +
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
    labs(title = m, x = "True protein (Young healthy BM)", y = "Prediction") +    #test
    geom_text(
      data = metricsD %>% filter(model == m),
      aes(x = x, y = y, label = lbl),
      inherit.aes = FALSE,
      hjust = 1.05, vjust = -3.0, size = 3.4,
      fontface = "bold"
    ) +
    base_theme
}

# ======================================
# Panel C D E and F
# ======================================

pC <- make_panel_pred("Global")
pD <- make_panel_pred("FMLE")
p_ctpnet   <- make_panel_pred("cTPnet")
p_sclinear <- make_panel_pred("scLinear")
# ----------------------------
# Build UMAP df for train/test only (aligned set)
# ----------------------------

umapA <- Embeddings(seu_prepA, "ref.umap")
stopifnot(!is.null(rownames(umapA)))
stopifnot(all(test_cellsA %in% rownames(umapA)))

umap_dfA <- data.frame(
  cell     = test_cellsA,
  UMAP_1   = umapA[test_cellsA, 1],
  UMAP_2   = umapA[test_cellsA, 2],
  celltype = seu_prepA$predicted.celltype.l1[test_cellsA],
  stringsAsFactors = FALSE
)

# ----------------------------
# (p1–p4) Cross-data generalization panels (UMAP + entropy + compositions)
# ----------------------------
fit <- obj$fit

# Predict regimes on C-train and A-test using the SAME fit

stopifnot(!is.null(obj$X_cols))
rm(obj)
gc()
pred_tr <- FMLE::fmle_predict(fit, X_new=XC_train, Z_new=ZC_train, return_se=TRUE)
pred_te <- FMLE::fmle_predict(fit, X_new=XA_test,  Z_new=ZA_test,  return_se=TRUE)

alpha_te <- pred_te$alpha
alpha_tr <- pred_tr$alpha
if (is.null(rownames(alpha_tr))) rownames(alpha_tr) <- rownames(XC_train)
if (is.null(rownames(alpha_te))) rownames(alpha_te) <- rownames(XA_test)

stopifnot(setequal(rownames(alpha_tr), train_cellsC))
stopifnot(setequal(rownames(alpha_te), test_cellsA))

# p1/p2: test donor UMAP colored by regime + entropy
reg_te <- max.col(alpha_te, ties.method="first")
Hn_te  <- entropy_norm_rows(alpha_te)

tapply(Hn_te, hard, summary)
tapply(Hn_te, seu_prepA$predicted.celltype.l1[cells], summary)

Hn_tr  <- entropy_norm_rows(alpha_tr)

names(Hn_tr) <- rownames(alpha_tr)

hard_tr <- max.col(alpha_tr, ties.method="first")
names(hard_tr) <- rownames(alpha_tr)

ct_tr <- seu_prepC@meta.data[rownames(alpha_tr), "cell_type", drop=TRUE]

tapply(Hn_tr[names(hard_tr)], factor(hard_tr), summary)
tapply(Hn_tr[rownames(alpha_tr)], factor(ct_tr), summary)

summary(Hn_tr)
summary(Hn_te)

# =======================================
# Train vs test entropy
entropy_shift_stats(Hn_tr, Hn_te)
# =======================================

pred_df_te <- data.frame(
  cell    = rownames(alpha_te),
  regime  = factor(reg_te),
  entropy = as.numeric(Hn_te),
  stringsAsFactors = FALSE
)

plot_df_te <- umap_dfA %>%
  left_join(pred_df_te, by="cell")

stopifnot(nrow(plot_df_te) == length(test_cellsA))
stopifnot(!anyNA(plot_df_te$regime), !anyNA(plot_df_te$entropy))

# ======================================
# Panel G
# ======================================
p1 <- ggplot(plot_df_te, aes(UMAP_1, UMAP_2, color = regime)) +
  geom_point(size = 0.25, alpha = 0.7) +
  theme_classic(base_size = 8) +
  labs(title="Test on Young healthy BM UMAP:\ntransferred regime (argmax α)", color="Regime")+
  theme(plot.title = element_text(hjust = 0.5, vjust = -2, face = "bold"))

reg_tr <- max.col(alpha_tr, ties.method="first")

comp <- bind_rows(
  data.frame(split="BMMC CITE-seq", regime=factor(reg_tr), stringsAsFactors=FALSE),
  data.frame(split="Young healthy BM",  regime=factor(reg_te), stringsAsFactors=FALSE)
) %>%
  count(split, regime, name="n") %>%
  group_by(split) %>%
  mutate(frac = n/sum(n)) %>%
  ungroup()

comp$split <- factor(comp$split, levels = c("BMMC CITE-seq", "Young healthy BM"))

# JSD on matched regime support
regs <- sort(unique(comp$regime))
pvec <- comp %>% filter(split=="BMMC CITE-seq") %>%        #train
  right_join(data.frame(regime=regs), by="regime") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac)) %>% arrange(regime) %>% pull(frac)
qvec <- comp %>% filter(split=="Young healthy BM") %>%            #test 
  right_join(data.frame(regime=regs), by="regime") %>%
  mutate(frac = ifelse(is.na(frac), 0, frac)) %>% arrange(regime) %>% pull(frac)

JSD <- jsd(pvec, qvec, base=2)
cat(sprintf("JSD(train: BMMC CITE-seq,test: Young healthy BM)=%.6g\n", JSD))

# p4: regime within celltypes, faceted train vs test
min_n <- 100

map_A_l2_to_shared <- function(x) {
  dplyr::case_when(
    x %in% c("Memory B", "Naive B", "transitional B", "Plasma", "pre B", "pro B") ~ "B",
    
    x %in% c("CD4 Effector", "CD4 Memory", "CD4 Naive") ~ "CD4 T",
    
    x %in% c("CD8 Effector_1", "CD8 Effector_2", "CD8 Effector_3",
             "CD8 Memory", "CD8 Naive") ~ "CD8 T",
    
    x %in% c("CD14 Mono", "CD16 Mono") ~ "Mono",
    
    x %in% c("cDC2", "pDC", "pre-mDC", "pre-pDC", "ASDC") ~ "DC",
    
    x %in% c("HSC", "CLP", "LMPP", "GMP", "EMP", "Prog Mk") ~ "HSPC",
    
    x %in% c("Early Eryth", "Late Eryth") ~ "Erythroid",
    
    x %in% c("NK", "NK CD56+") ~ "NK",
    
    x %in% c("MAIT", "T Proliferating") ~ "other T",
    
    x %in% c("ILC") ~ "other innate",
    
    TRUE ~ "other"
  )
}


map_C_to_shared <- function(x) {
  dplyr::case_when(
    x %in% c(
      "B1 B IGKC-", "B1 B IGKC+",
      "Naive CD20+ B IGKC-", "Naive CD20+ B IGKC+",
      "Transitional B",
      "Plasma cell IGKC-", "Plasma cell IGKC+",
      "Plasmablast IGKC-", "Plasmablast IGKC+"
    ) ~ "B",
    
    x %in% c(
      "CD4+ T activated", "CD4+ T activated integrinB7+",
      "CD4+ T CD314+ CD45RA+", "CD4+ T naive", "T reg"
    ) ~ "CD4 T",
    
    x %in% c(
      "CD8+ T CD49f+", "CD8+ T CD57+ CD45RA+",
      "CD8+ T CD57+ CD45RO+", "CD8+ T CD69+ CD45RA+",
      "CD8+ T CD69+ CD45RO+", "CD8+ T naive",
      "CD8+ T naive CD127+ CD26- CD101-",
      "CD8+ T TIGIT+ CD45RA+", "CD8+ T TIGIT+ CD45RO+"
    ) ~ "CD8 T",
    
    x %in% c("CD14+ Mono", "CD16+ Mono") ~ "Mono",
    
    x %in% c("cDC1", "cDC2", "pDC") ~ "DC",
    
    x %in% c("HSC", "Lymph prog", "G/M prog", "MK/E prog", "T prog cycling") ~ "HSPC",
    
    x %in% c("Erythroblast", "Normoblast", "Proerythroblast", "Reticulocyte") ~ "Erythroid",
    
    x %in% c("NK", "NK CD158e1+") ~ "NK",
    
    x %in% c("dnT", "gdT CD158b+", "gdT TCRVD2+", "MAIT") ~ "other T",
    
    x %in% c("ILC", "ILC1") ~ "other innate",
    
    TRUE ~ "other"
  )
}

seu_prepA$celltype_shared <- map_A_l2_to_shared(seu_prepA$predicted.celltype.l2)
seu_prepC$celltype_shared <- map_C_to_shared(seu_prepC$cell_type)

ct_col <- "celltype_shared"

stopifnot(ct_col %in% colnames(seu_prepA@meta.data))
stopifnot(ct_col %in% colnames(seu_prepC@meta.data))

comp_tr_ct <- celltype_regime_comp(alpha_tr, seu_prepC, cell_ids=train_cellsC,
                                   celltype_col=ct_col, min_n=min_n, split_label="train: BMMC CITE-seq")
comp_te_ct <- celltype_regime_comp(alpha_te, seu_prepA, cell_ids=test_cellsA,
                                   celltype_col=ct_col, min_n=min_n, split_label="test: Young healthy BM")
comp_all_ct <- bind_rows(comp_tr_ct, comp_te_ct)

comp_all_ct2 <- comp_all_ct %>%
  dplyr::filter(!celltype %in% c("unknown", "Unknown", "other"), !is.na(celltype))



regime_colors <- c("1" = "#E8806A",   # salmon
                   "2" = "#6DAF4B",   # green  
                   "3" = "#3BBFBF",   # teal
                   "4" = "#B07FD4")   # purple

# ======================================
# Panel H
# ======================================


cd33_celltypes <- c("HSPC", "Mono", "DC", "NK")
comp_cd33 <- comp_all_ct2 %>%
  dplyr::filter(celltype %in% cd33_celltypes)


regime_colors <- c("1" = "#E8806A",   # salmon
                   "2" = "#6DAF4B",   # green  
                   "3" = "#3BBFBF",   # teal
                   "4" = "#B07FD4")   # purple

p_cd33 <- ggplot(comp_cd33, 
                 aes(x = celltype, y = frac, fill = regime)) +
  geom_col(width = 0.75) +
  facet_wrap(~ split, nrow = 1) +
  scale_fill_manual(values = regime_colors,
                    name = "Regime") +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  scale_x_discrete(limits = cd33_celltypes) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Cross-dataset CD33 regime composition\nin shared myeloid and progenitor cell types",
    x = NULL,
    y = "Fraction"
  )


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

xt_tr <- make_ctab(alpha_tr, seu_prepC, train_cellsC, celltype_col=ct_col, min_n=min_n)
xt_te <- make_ctab(alpha_te, seu_prepA, test_cellsA,  celltype_col=ct_col, min_n=min_n)

cat("\nTRAIN_BMMC CITE-seq: chisq + Cramer's V\n")
print(chisq.test(xt_tr, simulate.p.value=TRUE, B=20000))
print(cramerV(xt_tr))

cat("\nTEST_Young healthy BM: chisq + Cramer's V\n")
print(chisq.test(xt_te, simulate.p.value=TRUE, B=20000))
print(cramerV(xt_te))


# =======================================
# Panel I UMAP
# =======================================


lineage_use <- "HSPC" 

ct_col <- "celltype_shared"

lineages_to_test <- sort(unique(seu_prepA@meta.data[[ct_col]]))
lineage_tab <- bind_rows(lapply(lineages_to_test, screen_one_cross))
print(lineage_tab, n= Inf)

res_nk  <- lineage_stats_cross(
  seu_train   = seu_prepC,
  seu_test    = seu_prepA,
  alpha_tr    = alpha_tr,
  alpha_te    = alpha_te,
  train_cells = train_cellsC,
  test_cells  = test_cellsA,
  ct_col      = "celltype_shared",
  lineage     = lineage_use,
  assay       = "RNA",
  layer       = "data",
  min_cells   = 200
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

seu_prep <- seu_prepA
test_cells <- test_cellsA

meta <- seu_prepA@meta.data
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

p3_left_HSPC <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data = n_label,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.0, vjust = 2.0,
    size = 3.5, color = "black", fontface = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = lims, oob = scales::squish,
    name = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.4, 0.4),
    ylim = yr + c(-0.5, 0.5)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "HSPCs: CD33 progenitor–granulocytic \ncontinuum",
    subtitle = "Healthy young BM, test cells",
    x = "UMAP 1", y = "UMAP 2"
  )+ theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_HSPC

# =======================================
# Panel J Pathway
# =======================================

pathway_df_HSPC <- bind_rows(
  hallmark_concordant %>%
    filter(ID %in% c(
      "HALLMARK_MYC_TARGETS_V1",
      "HALLMARK_COMPLEMENT"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "HALLMARK_MYC_TARGETS_V1" ~ "Hallmark MYC targets",
        ID == "HALLMARK_COMPLEMENT"     ~ "Hallmark complement"
      ),
      NES_tr  = NES_tr,
      NES_te  = NES_te,
      padj_tr = padj_tr,
      padj_te = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    ),
  
  reactome_concordant %>%
    filter(ID %in% c(
      "REACTOME_TRANSLATION",
      "REACTOME_RRNA_PROCESSING",
      "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
      "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
      "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
      "REACTOME_NEUTROPHIL_DEGRANULATION"
    )) %>%
    transmute(
      pathway = case_when(
        ID == "REACTOME_TRANSLATION" ~ "Translation",
        ID == "REACTOME_RRNA_PROCESSING" ~ "rRNA processing",
        ID == "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES" ~ "Amino acid metabolism",
        ID == "REACTOME_CELLULAR_RESPONSE_TO_STARVATION" ~ "Starvation response",
        ID == "REACTOME_NONSENSE_MEDIATED_DECAY_NMD" ~ "NMD",
        ID == "REACTOME_NEUTROPHIL_DEGRANULATION" ~ "Neutrophil degranulation"
      ),
      NES_tr  = NES_tr,
      NES_te  = NES_te,
      padj_tr = padj_tr,
      padj_te = padj_te,
      fdr_lab = sprintf("FDR %.1e", padj_te)
    )
) %>%
  arrange(NES_te) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_HSPC$NES_te, na.rm = TRUE)
xmax <- max(pathway_df_HSPC$NES_te, na.rm = TRUE)

p3_right_HSPC <- ggplot(pathway_df_HSPC, aes(x = NES_te, y = pathway)) +
  geom_col(
    aes(fill = NES_te > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  )+
  geom_text(
    aes(
      label = fdr_lab,
      x = NES_te/2
    ),
    hjust = 0.5,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(xlim = c(min(0, xmin) + 0.14, xmax + 0.15), clip = "off") +
  theme_classic(base_size = 8) +
  labs(
    title = "HSPC biosynthetic and granulocytic\npathway shifts",
    subtitle = "Test dataset: healthy young BM",
    x        = "Test NES (healthy young BM)",
    y        = NULL
  )+
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_HSPC

# =======================================
# Panel K UMAP
# =======================================

lineage_use <- "Mono"
ct_col <- "celltype_shared"

lineages_to_test <- sort(unique(seu_prepA@meta.data[[ct_col]]))
lineage_tab <- bind_rows(lapply(lineages_to_test, screen_one))
print(lineage_tab, n= Inf)

res_nk  <- lineage_stats_cross(
  seu_train   = seu_prepC,
  seu_test    = seu_prepA,
  alpha_tr    = alpha_tr,
  alpha_te    = alpha_te,
  train_cells = train_cellsC,
  test_cells  = test_cellsA,
  ct_col      = "celltype_shared",
  lineage     = lineage_use,
  assay       = "RNA",
  layer       = "data",
  min_cells   = 200
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
seu_prep <- seu_prepA
test_cells <- test_cellsA

meta <- seu_prepA@meta.data
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

p3_left_mono <- ggplot(mono_df, aes(UMAP_1, UMAP_2, color = score)) +
  geom_point(size = 0.8, alpha = 0.85) +
  geom_text(
    data = n_label,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.40, vjust = 2.0,
    size = 3.5, color = "black", fontface = "bold"
  ) +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = lims, oob = scales::squish,
    name = "FMLE\nscore"
  ) +
  coord_cartesian(
    xlim = xr + c(-0.4, 0.4),
    ylim = yr + c(-0.5, 0.5)
  ) +
  theme_classic(base_size = 8) +
  labs(
    title = "Monocytes: inflammatory–antigen-\npresentation CD33 axis",
    subtitle = "Healthy young BM, test cells",
    x = "UMAP 1", y = "UMAP 2"
  )+ theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 8),
    plot.margin     = margin(3, 3, 3, 3)
  )

p3_left_mono

# =======================================
# Panel L Pathway
# =======================================

pathway_df_mono <- reactome_concordant %>%
  filter(ID %in% c(
    "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE",
    "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES"
  )) %>%
  transmute(
    pathway = case_when(
      ID == "REACTOME_TOLL_LIKE_RECEPTOR_TLR1_TLR2_CASCADE" ~ "TLR1/2 cascade",
      ID == "REACTOME_TOLL_LIKE_RECEPTOR_CASCADES" ~ "TLR cascades"
    ),
    NES_tr  = NES_tr,
    NES_te  = NES_te,
    padj_tr = padj_tr,
    padj_te = padj_te,
    fdr_lab = sprintf("FDR %.1e", padj_te)
  ) %>%
  arrange(desc(NES_te)) %>%
  mutate(pathway = factor(pathway, levels = pathway))

xmin <- min(pathway_df_mono$NES_te, na.rm = TRUE)
xmax <- max(pathway_df_mono$NES_te, na.rm = TRUE)

p3_right_mono <- ggplot(pathway_df_mono, aes(x = NES_te, y = pathway)) +
  geom_col(
    aes(fill = NES_te > 0),
    width = 0.65,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D")
  )+
  geom_text(
    aes(
      label = fdr_lab,
      x = NES_te/2
    ),
    hjust = 0.5,
    size = 2.8,
    color = "white",
    fontface = "bold"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  coord_cartesian(xlim = c(min(0, xmin) + 0.09, xmax + 0.15), clip = "off") +
  theme_classic(base_size = 8) +
  labs(
    title = "Monocyte innate inflammatory\nTLR signalling activation",
    subtitle = "Test dataset: healthy young BM",
    x        = "Test NES (healthy young BM)",
    y        = NULL
  )+
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

p3_right_mono



tag_theme <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.02, 0.98)
)
tag_theme1 <- theme(
  plot.tag = element_text(face = "bold"),
  plot.tag.position = c(0.01, 0.97)
)

pA         <- pA         + labs(tag="(a)") + tag_theme
pB         <- pB         + labs(tag="(b)") + tag_theme
pC         <- pC         + labs(tag="(c)") + tag_theme
p_ctpnet   <- p_ctpnet   + labs(tag="(d)") + tag_theme
p_sclinear <- p_sclinear + labs(tag="(e)") + tag_theme1
pD         <- pD         + labs(tag="(f)") + tag_theme1
p1         <- p1         + labs(tag="(g)") + tag_theme1
p_cd33          <- p_cd33          + labs(tag="(h)") + tag_theme1
p3_left_HSPC  <- p3_left_HSPC + labs(tag="(i)") + tag_theme1
p3_right_HSPC <- p3_right_HSPC + labs(tag="(j)") + tag_theme1
p3_left_mono  <- p3_left_mono + labs(tag="(k)") + tag_theme1
p3_right_mono <- p3_right_mono + labs(tag="(l)") + tag_theme1


fig <- 
  (pA | pB | pC | p_ctpnet) /
  (p_sclinear |pD | p1 | p_cd33) /
  (p3_left_HSPC | p3_right_HSPC| p3_left_mono | p3_right_mono) + 
  plot_layout(heights = c(1, 1, 1)) &
  theme(
    plot.title = element_text(face = "bold"),
    plot.tag = element_text(face = "bold"),
    legend.position = "right"
  ) 

fig

dir.create(out_dir, recursive = TRUE)
ggsave(file.path(out_dir, "Supplementary_Figure_6.pdf"),
       plot = fig,
       device = cairo_pdf,
       width = 13.7, height = 9, units = "in")








