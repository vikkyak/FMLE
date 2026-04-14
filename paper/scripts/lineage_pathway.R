suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(clusterProfiler)
  library(msigdbr)
})


#=======================================
# 1) Add hard regime labels to metadata
#=======================================
add_regime_meta <- function(seu, alpha, prefix = c("tr", "te")) {
  prefix <- match.arg(prefix)
  alpha <- as.matrix(alpha)
  
  cells <- intersect(colnames(seu), rownames(alpha))
  alpha <- alpha[cells, , drop = FALSE]
  
  hard_id  <- max.col(alpha, ties.method = "first")
  hard_lab <- colnames(alpha)[hard_id]
  
  seu[[paste0(prefix, "_regime")]]    <- setNames(hard_lab, cells)
  seu[[paste0(prefix, "_regime_id")]] <- setNames(hard_id,  cells)
  seu
}

# --------------------------------------------------
# 2) Hard-regime DE within one lineage
# --------------------------------------------------
run_lineage_de <- function(
    seu_prep, split_cells, lineage, regime_col, ct_col,
    assay = "RNA", layer = "data", min_cells = 200,
    min.pct = 0.05, logfc.threshold = 0
) {
  md <- seu_prep@meta.data
  
  cells_lin <- intersect(
    split_cells,
    rownames(md)[md[[ct_col]] == lineage]
  )
  
  tab <- table(md[cells_lin, regime_col, drop = TRUE])
  
  # auto-adjust cutoff
  min_cells_use <- min(min_cells, min(tab))
  min_cells_use <- max(10, min_cells_use)
  
  keep_reg <- names(tab)[tab >= min_cells_use]
  
  if (length(keep_reg) < 2) {
    message(
      "Skipping lineage='", lineage, "' for ", regime_col,
      " : need >=2 regimes with at least min_cells=", min_cells_use,
      ". Found: ", paste(paste(names(tab), as.integer(tab), sep=":"), collapse=", ")
    )
    return(NULL)
  }
  
  keep_cells <- cells_lin[md[cells_lin, regime_col, drop = TRUE] %in% keep_reg]
  
  obj <- subset(seu_prep, cells = keep_cells)
  Idents(obj) <- factor(obj@meta.data[[regime_col]], levels = keep_reg)
  
  regs <- levels(Idents(obj))
  out <- list()
  
  for (i in seq_along(regs)) {
    for (j in seq_along(regs)) {
      if (i >= j) next
      
      r1 <- regs[i]
      r2 <- regs[j]
      
      mk <- FindMarkers(
        object = obj,
        ident.1 = r1,
        ident.2 = r2,
        assay = assay,
        slot = layer,
        test.use = "wilcox",
        min.pct = min.pct,
        logfc.threshold = logfc.threshold,
        verbose = FALSE
      )
      
      if (nrow(mk) == 0) next
      
      mk <- mk %>%
        rownames_to_column("gene") %>%
        # filter(
        #   pmax(pct.1, pct.2) >= 0.10,
        #   abs(pct.1 - pct.2) >= 0.10,
        #   abs(avg_log2FC) >= 0.25
        # ) %>%
        filter(
          pct.1 >= 0.10,
          pct.2 >= 0.10,
          abs(pct.1 - pct.2) >= 0.10,
          abs(avg_log2FC) >= 0.25
        )%>%
        mutate(
          contrast = paste0(r1, "_vs_", r2),
          p_bh = p.adjust(p_val, method = "BH")
        )
      
      if (nrow(mk) == 0) next
      
      out[[paste0(r1, "_vs_", r2)]] <- mk
    }
  }
  
  res <- if (length(out) == 0) tibble() else bind_rows(out)
  
  list(
    counts = tab,
    min_cells_used = min_cells_use,
    kept_regimes = keep_reg,
    res = res,
    summary = if (nrow(res) == 0) {
      tibble(n = 0L, min_p = NA_real_, min_bh = NA_real_, n_fdr = 0L)
    } else {
      res %>% summarise(
        n = n(),
        min_p = min(p_val, na.rm = TRUE),
        min_bh = min(p_bh, na.rm = TRUE),
        n_fdr = sum(p_bh < 0.05, na.rm = TRUE)
      )
    }
  )
}


#====================================================
# 3) Continuous score association within one lineage
#====================================================
run_lineage_score_assoc <- function(
    seu_prep, alpha, split_cells, lineage, ct_col,
    score_regimes = NULL,
    assay = "RNA", layer = "data", min_pct = 0.05
) {
  md <- seu_prep@meta.data
  
  cells_lin <- intersect(
    split_cells,
    rownames(md)[md[[ct_col]] == lineage]
  )
  
  if (length(cells_lin) == 0) {
    message("Skipping lineage='", lineage, "' : no cells found")
    return(NULL)
  }
  
  A <- as.matrix(alpha[cells_lin, , drop = FALSE])
  
  if (nrow(A) == 0 || ncol(A) == 0) {
    message("Skipping lineage='", lineage, "' : no alpha values found")
    return(NULL)
  }
  
  if (is.null(score_regimes)) {
    hard <- colnames(A)[max.col(A, ties.method = "first")]
    tab  <- sort(table(hard), decreasing = TRUE)
    
    if (length(tab) < 2) {
      message("Skipping lineage='", lineage, "' : only 1 regime present")
      return(NULL)
    }
    
    score_regimes <- names(tab)[1:2]
  }
  
  if (!all(score_regimes %in% colnames(A))) {
    message("Skipping lineage='", lineage, "' : requested regimes not present")
    return(NULL)
  }
  
  A <- A[, score_regimes, drop = FALSE]
  A <- A[rowSums(A) > 0, , drop = FALSE]
  
  if (nrow(A) < 20) {
    message("Skipping lineage='", lineage, "' : too few cells after regime filtering")
    return(NULL)
  }
  
  A <- A / rowSums(A)
  
  score <- A[, 1] - A[, 2]
  names(score) <- rownames(A)
  
  X <- GetAssayData(seu_prep, assay = assay, layer = layer)[, names(score), drop = FALSE]
  keep_genes <- rownames(X)[Matrix::rowMeans(X > 0) >= min_pct]
  X <- X[keep_genes, , drop = FALSE]
  
  assoc <- t(apply(X, 1, function(x) {
    ok <- is.finite(x) & is.finite(score)
    x <- as.numeric(x[ok])
    y <- score[ok]
    
    if (length(x) < 20 || sd(x) == 0 || sd(y) == 0) {
      return(c(rho = NA_real_, p = NA_real_))
    }
    
    ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }))
  
  assoc <- as.data.frame(assoc) %>%
    rownames_to_column("gene") %>%
    as_tibble() %>%
    mutate(p_bh = p.adjust(p, method = "BH")) %>%
    arrange(p_bh, p)
  
  list(
    assoc = assoc,
    summary = assoc %>% summarise(
      n = n(),
      min_p = min(p, na.rm = TRUE),
      min_bh = min(p_bh, na.rm = TRUE),
      n_fdr = sum(p_bh < 0.05, na.rm = TRUE)
    )
  )
}



#===================================
# 4) Train-vs-test reproducibility
#===================================

compare_hard_de <- function(de_tr, de_te) {
  if (is.null(de_tr) || is.null(de_te) || nrow(de_tr) == 0 || nrow(de_te) == 0) {
    return(tibble(n = NA_integer_, cor_lfc = NA_real_, prop_same_sign = NA_real_))
  }
  
  cmp <- de_tr %>%
    select(gene, contrast, lfc_tr = avg_log2FC) %>%
    inner_join(
      de_te %>% select(gene, contrast, lfc_te = avg_log2FC),
      by = c("gene", "contrast")
    ) %>%
    mutate(same_sign = sign(lfc_tr) == sign(lfc_te))
  
  if (nrow(cmp) == 0) {
    return(tibble(n = 0L, cor_lfc = NA_real_, prop_same_sign = NA_real_))
  }
  
  tibble(
    n = nrow(cmp),
    cor_lfc = cor(cmp$lfc_tr, cmp$lfc_te, method = "spearman"),
    prop_same_sign = mean(cmp$same_sign)
  )
}

compare_score_assoc <- function(a_tr, a_te) {
  if (is.null(a_tr) || is.null(a_te) || nrow(a_tr) == 0 || nrow(a_te) == 0) {
    return(tibble(n = NA_integer_, cor_rho = NA_real_, prop_same_sign = NA_real_))
  }
  
  cmp <- a_tr %>%
    select(gene, rho_tr = rho) %>%
    inner_join(
      a_te %>% select(gene, rho_te = rho),
      by = "gene"
    ) %>%
    mutate(same_sign = sign(rho_tr) == sign(rho_te))
  
  if (nrow(cmp) == 0) {
    return(tibble(n = 0L, cor_rho = NA_real_, prop_same_sign = NA_real_))
  }
  
  tibble(
    n = nrow(cmp),
    cor_rho = cor(cmp$rho_tr, cmp$rho_te, method = "spearman"),
    prop_same_sign = mean(cmp$same_sign)
  )
}


lineage_stats <- function(
    seu_prep, alpha_tr, alpha_te, train_cells, test_cells,
    ct_col, lineage, score_regimes = NULL,
    assay = "RNA", layer = "data", min_cells = 200
) {
  seu_prep <- add_regime_meta(seu_prep, alpha_tr, "tr")
  seu_prep <- add_regime_meta(seu_prep, alpha_te, "te")
  
  # ----------------------------
  # hard DE
  # ----------------------------
  de_tr <- run_lineage_de(
    seu_prep, train_cells, lineage, "tr_regime", ct_col,
    assay = assay, layer = layer, min_cells = min_cells
  )
  
  de_te <- run_lineage_de(
    seu_prep, test_cells, lineage, "te_regime", ct_col,
    assay = assay, layer = layer, min_cells = min_cells
  )
  
  hard_repro <- if (!is.null(de_tr) && !is.null(de_te)) {
    compare_hard_de(de_tr$res, de_te$res)
  } else {
    tibble(n = NA_integer_, cor_lfc = NA_real_, prop_same_sign = NA_real_)
  }
  
  # ----------------------------
  # choose score_regimes ONCE from TRAIN if NULL
  # ----------------------------
  score_regimes_used <- score_regimes
  
  if (is.null(score_regimes_used)) {
    md <- seu_prep@meta.data
    
    cells_lin_tr <- intersect(
      train_cells,
      rownames(md)[md[[ct_col]] == lineage]
    )
    
    if (length(cells_lin_tr) > 0) {
      A_tr <- as.matrix(alpha_tr[cells_lin_tr, , drop = FALSE])
      A_tr <- A_tr[rowSums(A_tr) > 0, , drop = FALSE]
      
      if (nrow(A_tr) > 0 && ncol(A_tr) >= 2) {
        hard_tr <- colnames(A_tr)[max.col(A_tr, ties.method = "first")]
        tab_tr <- sort(table(hard_tr), decreasing = TRUE)
        
        if (length(tab_tr) >= 2) {
          score_regimes_used <- names(tab_tr)[1:2]
        }
      }
    }
  }
  
  # ----------------------------
  # continuous score association
  # ----------------------------
  assoc_tr <- NULL
  assoc_te <- NULL
  
  if (!is.null(score_regimes_used) && length(score_regimes_used) == 2) {
    assoc_tr <- run_lineage_score_assoc(
      seu_prep, alpha_tr, train_cells, lineage, ct_col,
      score_regimes = score_regimes_used, assay = assay, layer = layer
    )
    
    assoc_te <- run_lineage_score_assoc(
      seu_prep, alpha_te, test_cells, lineage, ct_col,
      score_regimes = score_regimes_used, assay = assay, layer = layer
    )
  }
  
  cont_repro <- if (!is.null(assoc_tr) && !is.null(assoc_te)) {
    compare_score_assoc(assoc_tr$assoc, assoc_te$assoc)
  } else {
    tibble(n = NA_integer_, cor_rho = NA_real_, prop_same_sign = NA_real_)
  }
  
  list(
    lineage = lineage,
    score_regimes = score_regimes_used,
    de_tr = de_tr,
    de_te = de_te,
    hard_repro = hard_repro,
    assoc_tr = assoc_tr,
    assoc_te = assoc_te,
    cont_repro = cont_repro
  )
}


lineage_stats_cross <- function(
    seu_train, seu_test,
    alpha_tr, alpha_te,
    train_cells, test_cells,
    ct_col, lineage, score_regimes = NULL,
    assay = "RNA", layer = "data", min_cells = 200
) {
  seu_train <- add_regime_meta(seu_train, alpha_tr, "tr")
  seu_test  <- add_regime_meta(seu_test,  alpha_te, "te")
  
  de_tr <- run_lineage_de(
    seu_train, train_cells, lineage, "tr_regime", ct_col,
    assay = assay, layer = layer, min_cells = min_cells
  )
  
  de_te <- run_lineage_de(
    seu_test, test_cells, lineage, "te_regime", ct_col,
    assay = assay, layer = layer, min_cells = min_cells
  )
  
  hard_repro <- if (!is.null(de_tr) && !is.null(de_te)) {
    compare_hard_de(de_tr$res, de_te$res)
  } else {
    tibble(n = NA_integer_, cor_lfc = NA_real_, prop_same_sign = NA_real_)
  }
  
  score_regimes_used <- score_regimes
  
  if (is.null(score_regimes_used)) {
    md <- seu_train@meta.data
    cells_lin_tr <- intersect(train_cells, rownames(md)[md[[ct_col]] == lineage])
    
    if (length(cells_lin_tr) > 0) {
      A_tr <- as.matrix(alpha_tr[cells_lin_tr, , drop = FALSE])
      A_tr <- A_tr[rowSums(A_tr) > 0, , drop = FALSE]
      
      if (nrow(A_tr) > 0 && ncol(A_tr) >= 2) {
        hard_tr <- colnames(A_tr)[max.col(A_tr, ties.method = "first")]
        tab_tr <- sort(table(hard_tr), decreasing = TRUE)
        if (length(tab_tr) >= 2) score_regimes_used <- names(tab_tr)[1:2]
      }
    }
  }
  
  assoc_tr <- NULL
  assoc_te <- NULL
  
  if (!is.null(score_regimes_used) && length(score_regimes_used) == 2) {
    assoc_tr <- run_lineage_score_assoc(
      seu_train, alpha_tr, train_cells, lineage, ct_col,
      score_regimes = score_regimes_used, assay = assay, layer = layer
    )
    
    assoc_te <- run_lineage_score_assoc(
      seu_test, alpha_te, test_cells, lineage, ct_col,
      score_regimes = score_regimes_used, assay = assay, layer = layer
    )
  }
  
  cont_repro <- if (!is.null(assoc_tr) && !is.null(assoc_te)) {
    compare_score_assoc(assoc_tr$assoc, assoc_te$assoc)
  } else {
    tibble(n = NA_integer_, cor_rho = NA_real_, prop_same_sign = NA_real_)
  }
  
  list(
    lineage = lineage,
    score_regimes = score_regimes_used,
    de_tr = de_tr,
    de_te = de_te,
    hard_repro = hard_repro,
    assoc_tr = assoc_tr,
    assoc_te = assoc_te,
    cont_repro = cont_repro
  )
}



# ==============================================================================

screen_one <- function(lineage_name) {
  z <- lineage_stats(
    seu_prep = seu_prep,
    alpha_tr = alpha_tr,
    alpha_te = alpha_te,
    train_cells = train_cells,
    test_cells = test_cells,
    ct_col = ct_col,
    lineage = lineage_name,
    score_regimes = NULL, 
    assay = "RNA",
    layer = "data",
    min_cells = 200
  )
  
  # =====================
  #  Cross dataset
  # z<- lineage_stats_cross(
  #   seu_train   = seu_prepC,
  #   seu_test    = seu_prepA,
  #   alpha_tr    = alpha_tr,
  #   alpha_te    = alpha_te,
  #   train_cells = train_cellsC,
  #   test_cells  = test_cellsA,
  #   ct_col      = "celltype_shared",
  #   lineage     = lineage_name,
  #   assay       = "RNA",
  #   layer       = "data",
  #   min_cells   = 200
  # )
  
  tibble(
    lineage = lineage_name,
    
    hard_ok_tr = !is.null(z$de_tr),
    hard_ok_te = !is.null(z$de_te),
    
    hard_min_bh_tr = if (!is.null(z$de_tr)) z$de_tr$summary$min_bh else NA_real_,
    hard_min_bh_te = if (!is.null(z$de_te)) z$de_te$summary$min_bh else NA_real_,
    hard_cor_lfc  = z$hard_repro$cor_lfc,
    
    cont_min_bh_tr = z$assoc_tr$summary$min_bh,
    cont_min_bh_te = z$assoc_te$summary$min_bh,
    cont_cor_rho   = z$cont_repro$cor_rho
  )
}

# =====================
#  Cross dataset

screen_one_cross <- function(lineage_name) {

  z<- lineage_stats_cross(
    seu_train   = seu_prepC,
    seu_test    = seu_prepA,
    alpha_tr    = alpha_tr,
    alpha_te    = alpha_te,
    train_cells = train_cellsC,
    test_cells  = test_cellsA,
    ct_col      = "celltype_shared",
    lineage     = lineage_name,
    assay       = "RNA",
    layer       = "data",
    min_cells   = 200
  )

  tibble(
    lineage = lineage_name,
    
    hard_ok_tr = !is.null(z$de_tr),
    hard_ok_te = !is.null(z$de_te),
    
    hard_min_bh_tr = if (!is.null(z$de_tr)) z$de_tr$summary$min_bh else NA_real_,
    hard_min_bh_te = if (!is.null(z$de_te)) z$de_te$summary$min_bh else NA_real_,
    hard_cor_lfc  = z$hard_repro$cor_lfc,
    
    cont_min_bh_tr = z$assoc_tr$summary$min_bh,
    cont_min_bh_te = z$assoc_te$summary$min_bh,
    cont_cor_rho   = z$cont_repro$cor_rho
  )
}




# Pathway

# --------------------------------------------------
# 1) Build ranked gene list from continuous-score association
#    Input: res_mono$assoc_te$assoc or res_mono$assoc_tr$assoc
# --------------------------------------------------
make_ranked_list <- function(assoc_tbl) {
  assoc_tbl %>%
    filter(is.finite(rho), is.finite(p), !is.na(gene), gene != "") %>%
    mutate(rank_score = rho * (-log10(p + 1e-300))) %>%
    group_by(gene) %>%
    slice_max(order_by = abs(rank_score), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(rank_score)) %>%
    { setNames(.$rank_score, .$gene) }
}

# --------------------------------------------------
# 2) Run GSEA from MSigDB gene sets
# --------------------------------------------------
run_gsea_msig <- function(geneList, category = "H") {
  msig <- msigdbr(species = "Homo sapiens", category = category)
  
  term2gene <- msig %>%
    select(gs_name, gene_symbol)
  
  GSEA(
    geneList = sort(geneList, decreasing = TRUE),
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
}

run_gsea_reactome <- function(geneList) {
  msig <- msigdbr(
    species = "Homo sapiens",
    category = "C2",
    subcategory = "CP:REACTOME"
  )
  
  term2gene <- msig %>%
    select(gs_name, gene_symbol)
  
  GSEA(
    geneList = sort(geneList, decreasing = TRUE),
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    verbose = FALSE
  )
}

# --------------------------------------------------
# 3) Extract top positive / negative pathways
# --------------------------------------------------
top_gsea_terms <- function(gsea_obj, n = 15) {
  as.data.frame(gsea_obj) %>%
    arrange(p.adjust, desc(abs(NES))) %>%
    mutate(direction = ifelse(NES > 0, "positive_score", "negative_score")) %>%
    group_by(direction) %>%
    slice_head(n = n) %>%
    ungroup()
}

# --------------------------------------------------
# 4) Concordant pathways across train and test
# --------------------------------------------------
concordant_gsea <- function(gsea_tr, gsea_te, padj_cut = 0.05) {
  tr <- as.data.frame(gsea_tr) %>%
    filter(p.adjust < padj_cut) %>%
    select(ID, Description, NES_tr = NES, padj_tr = p.adjust)
  
  te <- as.data.frame(gsea_te) %>%
    filter(p.adjust < padj_cut) %>%
    select(ID, Description, NES_te = NES, padj_te = p.adjust)
  
  inner_join(tr, te, by = c("ID", "Description")) %>%
    mutate(
      same_sign = sign(NES_tr) == sign(NES_te),
      mean_abs_NES = (abs(NES_tr) + abs(NES_te)) / 2
    ) %>%
    filter(same_sign) %>%
    arrange(desc(mean_abs_NES), padj_te, padj_tr)
}


