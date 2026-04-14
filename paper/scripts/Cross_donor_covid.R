suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(clusterProfiler)
  library(IRanges)
  library(purrr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(FMLE)
  library(Azimuth)
})

source(file.path(here::here(), "paper", "scripts", "_config.R"))

# -------------------------
# 0) paths
# -------------------------
base_dir <- file.path(cfg$data_root, "GSE155224_RAW")
out_base <- file.path(cfg$out_root, "citeseq_v1")
out_ctp  <- file.path(cfg$out_root, "ctp")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ctp, recursive = TRUE, showWarnings = FALSE)

# Use all six metadata files if present. The current validation below is strongest for CITE1-CITE4.
meta_files <- c(
  CITE1 = file.path(base_dir, "CITE1_cite_metadata.txt"),
  CITE2 = file.path(base_dir, "CITE2_cite_metadata.txt"),
  CITE3 = file.path(base_dir, "CITE3_cite_metadata.txt"),
  CITE4 = file.path(base_dir, "CITE4_cite_metadata.txt"),
  CITE5 = file.path(base_dir, "CITE5_cite_metadata.txt"),
  CITE6 = file.path(base_dir, "CITE6_cite_metadata.txt")
)
meta_files <- meta_files[file.exists(meta_files)]
stopifnot(length(meta_files) > 0)

pooled_paths <- setNames(
  file.path(base_dir, paste0(names(meta_files), "_cite/filtered_feature_bc_matrix")),
  names(meta_files)
)
stopifnot(all(dir.exists(unname(pooled_paths))))

# Full HHT set used in your later corrected code.
hht_files <- c(
  NS0A = file.path(base_dir, "GSM4697611_NS0A_HHT_cellranger.tar.gz"),
  NS0B = file.path(base_dir, "GSM4697612_NS0B_HHT_cellranger.tar.gz"),
  NS1A = file.path(base_dir, "GSM4697613_NS1A_HHT_cellranger.tar.gz"),
  NS1B = file.path(base_dir, "GSM4697614_NS1B_HHT_cellranger.tar.gz"),
  TP6A = file.path(base_dir, "GSM4697615_TP6A_HHT_cellranger.tar.gz"),
  TP6B = file.path(base_dir, "GSM4697616_TP6B_HHT_cellranger.tar.gz"),
  TP7A = file.path(base_dir, "GSM4697617_TP7A_HHT_cellranger.tar.gz"),
  TP7B = file.path(base_dir, "GSM4697618_TP7B_HHT_cellranger.tar.gz"),
  TP9B = file.path(base_dir, "GSM4697620_TP9B_HHT_cellranger.tar.gz"),
  TS2A = file.path(base_dir, "GSM4697621_TS2A_HHT_cellranger.tar.gz"),
  TS2B = file.path(base_dir, "GSM4697622_TS2B_HHT_cellranger.tar.gz"),
  TS3A = file.path(base_dir, "GSM4697623_TS3A_HHT_cellranger.tar.gz"),
  TS3B = file.path(base_dir, "GSM4697624_TS3B_HHT_cellranger.tar.gz"),
  TS4A = file.path(base_dir, "GSM4697625_TS4A_HHT_cellranger.tar.gz"),
  TS4B = file.path(base_dir, "GSM4697626_TS4B_HHT_cellranger.tar.gz"),
  TS5A = file.path(base_dir, "GSM4697627_TS5A_HHT_cellranger.tar.gz"),
  TS5B = file.path(base_dir, "GSM4697628_TS5B_HHT_cellranger.tar.gz")
)
hht_files <- hht_files[file.exists(hht_files)]
stopifnot(length(hht_files) > 0)

subject_from_sample <- function(x) sub("[AB]$", "", x)

# -------------------------
# 1) metadata-driven pool truth
# -------------------------
read_pool_samples <- function(meta_file) {
  md <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  nm <- names(md)
  sample_col <- nm[nm %in% c("Sample name", "Sample.name")]
  stopifnot(length(sample_col) == 1)
  x <- as.character(md[[sample_col]])
  x <- x[!is.na(x) & nzchar(x)]
  x
}

pool_samples <- lapply(meta_files, read_pool_samples)
pool_expected_n <- sapply(pool_samples, length)
pool_subject_slots <- lapply(pool_samples, function(samps) table(subject_from_sample(samps)))

cat("\n=== METADATA-DERIVED POOL STRUCTURE ===\n")
for (pool in names(pool_samples)) {
  cat("\n", pool, "\n", sep = "")
  print(pool_samples[[pool]])
  print(pool_subject_slots[[pool]])
}

# -------------------------
# 2) HHT references
# -------------------------
load_hht_rna <- function(path, sample_name) {
  exdir <- tempfile(pattern = "hht_")
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  untar(path, exdir = exdir)

  subdirs <- list.dirs(exdir, recursive = FALSE, full.names = TRUE)
  stopifnot(length(subdirs) == 1)

  m <- Read10X(file.path(subdirs[1], "filtered_feature_bc_matrix"))
  if (is.list(m)) {
    stopifnot("Gene Expression" %in% names(m))
    rna <- m[["Gene Expression"]]
  } else {
    rna <- m
  }

  seu <- CreateSeuratObject(counts = rna, project = sample_name)
  seu$sample_true  <- sample_name
  seu$subject_true <- subject_from_sample(sample_name)
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
  seu
}

ref_all <- lapply(names(hht_files), function(s) load_hht_rna(hht_files[[s]], s))
names(ref_all) <- names(hht_files)

# -------------------------
# 3) pool-specific donor references
# -------------------------
make_pool_subject_ref_avg <- function(ref_all, candidate_samples) {
  stopifnot(all(candidate_samples %in% names(ref_all)))
  ref_sub <- ref_all[candidate_samples]
  sample_names <- names(ref_sub)
  subject_names <- sort(unique(subject_from_sample(sample_names)))

  out <- list()
  for (subj in subject_names) {
    samps <- sample_names[subject_from_sample(sample_names) == subj]
    mats <- lapply(samps, function(s) GetAssayData(ref_sub[[s]], assay = "RNA", layer = "data"))
    common_genes <- Reduce(intersect, lapply(mats, rownames))
    stopifnot(length(common_genes) > 100)
    avg_vecs <- lapply(mats, function(m) Matrix::rowMeans(m[common_genes, , drop = FALSE]))
    out[[subj]] <- Reduce(`+`, avg_vecs) / length(avg_vecs)
  }
  out
}

# -------------------------
# 4) pooled CITE loader + HTO demux
# -------------------------
load_pooled_cite_auto <- function(path, project, expected_n,
                                  min_hto_feature_total = 1000,
                                  min_hto_total = 5,
                                  hto_kfunc = "clara") {
  m <- Read10X(path)
  stopifnot("Gene Expression" %in% names(m))
  stopifnot("Antibody Capture" %in% names(m))

  rna <- m[["Gene Expression"]]
  adt <- m[["Antibody Capture"]]

  hto_rows_all <- grep("Hashtag", rownames(adt), value = TRUE)
  stopifnot(length(hto_rows_all) > 0)

  hto_totals <- Matrix::rowSums(adt[hto_rows_all, , drop = FALSE])
  hto_ranked <- names(sort(hto_totals, decreasing = TRUE))
  hto_ranked <- hto_ranked[hto_totals[hto_ranked] > min_hto_feature_total]
  stopifnot(length(hto_ranked) >= expected_n)

  hto_keep <- head(hto_ranked, expected_n)

  adt_rows <- setdiff(rownames(adt), hto_rows_all)
  hto <- adt[hto_keep, , drop = FALSE]

  keep_cells <- Matrix::colSums(hto) >= min_hto_total

  seu <- CreateSeuratObject(counts = rna[, keep_cells, drop = FALSE], project = project)
  seu$pool_id <- project
  seu[["ADT"]] <- CreateAssayObject(counts = adt[adt_rows, keep_cells, drop = FALSE])
  seu[["HTO"]] <- CreateAssayObject(counts = hto[, keep_cells, drop = FALSE])

  seu <- NormalizeData(seu, assay = "HTO", normalization.method = "CLR", margin = 2, verbose = FALSE)
  seu <- HTODemux(seu, assay = "HTO", init = nrow(seu[["HTO"]]) + 1, kfunc = hto_kfunc, verbose = FALSE)

  attr(seu, "hto_keep") <- hto_keep
  attr(seu, "hto_totals") <- sort(hto_totals, decreasing = TRUE)
  seu
}

subset_candidate_singlets <- function(seu, n_slots, min_cells = 30) {
  tab <- sort(table(as.character(seu$hash.ID[seu$HTO_classification.global == "Singlet"])), decreasing = TRUE)
  tab <- tab[tab >= min_cells]
  stopifnot(length(tab) >= n_slots)

  keep_hto <- names(tab)[seq_len(n_slots)]
  keep_idx <- seu$HTO_classification.global == "Singlet" &
    as.character(seu$hash.ID) %in% keep_hto

  seu_sub <- subset(seu, cells = colnames(seu)[keep_idx])
  attr(seu_sub, "candidate_hto_counts") <- tab
  attr(seu_sub, "kept_hto") <- keep_hto
  seu_sub
}

# -------------------------
# 5) HTO-group vs donor correlation matrix
# -------------------------
compute_subject_cor <- function(seu_pool, ref_subject_avg) {
  seu_pool <- NormalizeData(seu_pool, assay = "RNA", verbose = FALSE)
  seu_pool <- FindVariableFeatures(seu_pool, assay = "RNA", nfeatures = 3000, verbose = FALSE)

  hto_ids <- sort(unique(as.character(seu_pool$hash.ID)))
  candidate_subjects <- names(ref_subject_avg)

  mat <- GetAssayData(seu_pool, assay = "RNA", layer = "data")
  pool_avg <- lapply(hto_ids, function(h) {
    cells_h <- colnames(seu_pool)[as.character(seu_pool$hash.ID) == h]
    Matrix::rowMeans(mat[, cells_h, drop = FALSE])
  })
  names(pool_avg) <- hto_ids

  common_genes <- Reduce(intersect, c(
    list(VariableFeatures(seu_pool)),
    lapply(pool_avg, names),
    lapply(ref_subject_avg, names)
  ))
  if (length(common_genes) < 200) {
    common_genes <- Reduce(intersect, c(lapply(pool_avg, names), lapply(ref_subject_avg, names)))
  }
  stopifnot(length(common_genes) >= 200)

  cor_mat <- matrix(NA_real_,
                    nrow = length(pool_avg),
                    ncol = length(candidate_subjects),
                    dimnames = list(names(pool_avg), candidate_subjects))

  for (h in names(pool_avg)) {
    for (subj in candidate_subjects) {
      cor_mat[h, subj] <- cor(pool_avg[[h]][common_genes],
                              ref_subject_avg[[subj]][common_genes],
                              method = "spearman",
                              use = "pairwise.complete.obs")
    }
  }
  cor_mat
}

# -------------------------
# 6) slot-constrained donor assignment, but keep QC flags
# -------------------------
assign_subject_slots_qc <- function(cor_mat, slot_counts, ambiguous_gap = 0.01) {
  htos <- rownames(cor_mat)
  donors <- colnames(cor_mat)
  stopifnot(all(names(slot_counts) %in% donors))

  slots <- unlist(mapply(rep, names(slot_counts), as.integer(slot_counts), SIMPLIFY = FALSE))
  stopifnot(length(slots) == length(htos))

  best_score <- -Inf
  best_assign <- NULL
  used <- rep(FALSE, length(slots))
  current_assign <- character(length(htos))

  recurse <- function(i, score) {
    if (i > length(htos)) {
      if (score > best_score) {
        best_score <<- score
        best_assign <<- current_assign
      }
      return(invisible(NULL))
    }
    h <- htos[i]
    for (j in seq_along(slots)) {
      if (!used[j]) {
        used[j] <<- TRUE
        current_assign[i] <<- slots[j]
        recurse(i + 1, score + cor_mat[h, slots[j]])
        used[j] <<- FALSE
      }
    }
  }

  recurse(1, 0)

  out <- data.frame(
    hto = htos,
    subject_true = best_assign,
    hto_ncells = NA_integer_,
    cor_assigned = mapply(function(h, s) cor_mat[h, s], htos, best_assign),
    cor_best = apply(cor_mat[htos, , drop = FALSE], 1, max),
    cor_second = apply(cor_mat[htos, , drop = FALSE], 1, function(z) sort(z, decreasing = TRUE)[2]),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      cor_gap = cor_best - cor_second,
      loss_vs_best = cor_best - cor_assigned,
      forced = loss_vs_best > 0,
      ambiguous = cor_gap < ambiguous_gap,
      accept = is.finite(cor_assigned)
    )

  rownames(out) <- NULL
  out
}

attach_group_qc <- function(seu_pool, assign_df) {
  key_subj <- setNames(assign_df$subject_true, assign_df$hto)
  key_force <- setNames(assign_df$forced, assign_df$hto)
  key_amb <- setNames(assign_df$ambiguous, assign_df$hto)
  key_loss <- setNames(assign_df$loss_vs_best, assign_df$hto)
  key_gap <- setNames(assign_df$cor_gap, assign_df$hto)
  key_assigned <- setNames(assign_df$cor_assigned, assign_df$hto)

  h <- as.character(seu_pool$hash.ID)
  seu_pool$subject_true <- unname(key_subj[h])
  seu_pool$subject_forced <- unname(key_force[h])
  seu_pool$subject_ambiguous <- unname(key_amb[h])
  seu_pool$subject_loss_vs_best <- unname(key_loss[h])
  seu_pool$subject_cor_gap <- unname(key_gap[h])
  seu_pool$subject_cor_assigned <- unname(key_assigned[h])
  seu_pool
}

# -------------------------
# 7) diagnostics
# -------------------------
diagnose_pool <- function(pool, seu_pool, seu_sub, cor_mat, assign_df, slot_counts) {
  cat("\n====================\n")
  cat("POOL:", pool, "\n")
  cat("====================\n")

  cat("\nMetadata samples:\n")
  print(pool_samples[[pool]])

  cat("\nMetadata slot counts:\n")
  print(slot_counts)

  cat("\nHTO channels kept:\n")
  print(attr(seu_pool, "hto_keep"))

  cat("\nHTO classification:\n")
  print(table(seu_pool$HTO_classification.global))

  cat("\nCandidate singlet HTO counts:\n")
  print(attr(seu_sub, "candidate_hto_counts"))

  cat("\nKept HTO groups:\n")
  print(attr(seu_sub, "kept_hto"))

  cat("\nCorrelation matrix:\n")
  print(round(cor_mat, 3))

  cat("\nSlot-constrained assignment with QC:\n")
  print(assign_df)

  seu_lab <- attach_group_qc(seu_sub, assign_df)
  cat("\nRecovered subject counts:\n")
  print(table(seu_lab$subject_true, useNA = "ifany"))
}

# -------------------------
# 8) run pools
# -------------------------
run_one_pool <- function(pool_name,
                         min_hto_feature_total = 1000,
                         min_hto_total = 5,
                         min_cells = 30,
                         ambiguous_gap = 0.01) {
  candidate_samples <- pool_samples[[pool_name]]
  slot_counts <- pool_subject_slots[[pool_name]]
  n_slots <- sum(slot_counts)

  missing_refs <- setdiff(candidate_samples, names(ref_all))
  if (length(missing_refs) > 0) {
    stop(sprintf("Missing HHT refs for pool %s: %s", pool_name, paste(missing_refs, collapse = ", ")))
  }

  seu_pool <- load_pooled_cite_auto(
    path = pooled_paths[[pool_name]],
    project = pool_name,
    expected_n = n_slots,
    min_hto_feature_total = min_hto_feature_total,
    min_hto_total = min_hto_total
  )

  seu_sub <- subset_candidate_singlets(seu_pool, n_slots = n_slots, min_cells = min_cells)

  ref_subject_avg <- make_pool_subject_ref_avg(ref_all = ref_all, candidate_samples = candidate_samples)
  cor_mat <- compute_subject_cor(seu_pool = seu_sub, ref_subject_avg = ref_subject_avg)
  assign_df <- assign_subject_slots_qc(cor_mat = cor_mat, slot_counts = slot_counts, ambiguous_gap = ambiguous_gap)

  # attach HTO group sizes
  hto_tab <- table(as.character(seu_sub$hash.ID))
  assign_df$hto_ncells <- as.integer(hto_tab[assign_df$hto])

  diagnose_pool(pool_name, seu_pool, seu_sub, cor_mat, assign_df, slot_counts)
  seu_lab <- attach_group_qc(seu_sub, assign_df)

  list(
    seu_pool = seu_pool,
    seu_sub = seu_sub,
    cor_mat = cor_mat,
    assign_df = assign_df,
    seu_lab = seu_lab,
    slot_counts = slot_counts,
    candidate_samples = candidate_samples
  )
}

# -------------------------
# CITE1 + CITE2 only
# -------------------------
pool_names_12 <- names(pooled_paths)[1:2]
pool_results_12 <- setNames(lapply(pool_names_12, run_one_pool), pool_names_12)

collect_labeled_cells <- function(pool_results,
                                  max_loss_strict = 0.02,
                                  drop_ambiguous_strict = FALSE) {
  seu_list_all <- lapply(pool_results, function(x) {
    keep <- !is.na(x$seu_lab$subject_true)
    subset(x$seu_lab, cells = colnames(x$seu_lab)[keep])
  })
  
  seu_list_strict <- lapply(pool_results, function(x) {
    keep <- !is.na(x$seu_lab$subject_true) &
      !x$seu_lab$subject_forced &
      is.finite(x$seu_lab$subject_loss_vs_best) &
      x$seu_lab$subject_loss_vs_best <= max_loss_strict
    
    if (drop_ambiguous_strict) {
      keep <- keep & !x$seu_lab$subject_ambiguous
    }
    
    subset(x$seu_lab, cells = colnames(x$seu_lab)[keep])
  })
  
  out <- list()
  
  if (length(seu_list_all) == 1) {
    out$all <- seu_list_all[[1]]
  } else if (length(seu_list_all) > 1) {
    out$all <- merge(seu_list_all[[1]], y = seu_list_all[-1], add.cell.ids = names(seu_list_all))
  }
  
  if (length(seu_list_strict) == 1) {
    out$strict <- seu_list_strict[[1]]
  } else if (length(seu_list_strict) > 1) {
    out$strict <- merge(seu_list_strict[[1]], y = seu_list_strict[-1], add.cell.ids = names(seu_list_strict))
  }
  
  out
}

merged_12 <- collect_labeled_cells(
  pool_results_12,
  max_loss_strict = 0.02,
  drop_ambiguous_strict = FALSE
)

seu_12_all    <- merged_12$all
seu_12_strict <- merged_12$strict

cat("=== CITE1 + CITE2 ONLY ===\n")
print(table(seu_12_all$pool_id, seu_12_all$subject_true))
print(table(seu_12_all$subject_true))
print(table(seu_12_all$subject_forced, useNA = "ifany"))
print(table(seu_12_all$subject_ambiguous, useNA = "ifany"))

seu_all <- seu_12_all 

seu_all <- JoinLayers(seu_all, assay = "RNA")
Layers(seu_all[["RNA"]])

seu_all <- JoinLayers(seu_all, assay = "RNA")

seu_all[["nFeature_RNA"]] <- colSums(GetAssayData(seu_all, assay="RNA", layer="counts") > 0)
seu_all[["nCount_RNA"]]   <- colSums(GetAssayData(seu_all, assay="RNA", layer="counts"))
seu_all[["percent.mt"]]   <- PercentageFeatureSet(seu_all, pattern="^MT-")

print(summary(seu_all$nFeature_RNA))
print(summary(seu_all$nCount_RNA))
print(summary(seu_all$percent.mt))



f_lo <- 200      # minimum genes
f_hi <- 3619     # your 99.5th percentile
c_hi <- 25809    # your 99.5th percentile
mt_hi <- 10      # standard cutoff for PBMCs


seu_all <- subset(
  seu_all,
  nFeature_RNA >= f_lo &
    nFeature_RNA <= f_hi &
    nCount_RNA   <= c_hi &
    percent.mt   <= mt_hi 
)


adt_counts0 <- GetAssayData(seu_all, assay="ADT", layer="counts")
cs <- Matrix::colSums(adt_counts0)
cat("ADT total summary:\n")
cat("BEFORE filtering\n")
print(summary(cs))
cat("n cells:", length(cs), "\n")

thr <- as.numeric(quantile(cs, 0.995))

seu_all <- subset(seu_all, cells = names(cs)[cs <= thr])
adt_counts1 <- GetAssayData(seu_all, assay="ADT", layer="counts")
cs1 <- Matrix::colSums(adt_counts1)
cat("ADT total summary:\n")
cat("After filtering\n")
print(summary(cs1))
cat("n cells:", length(cs1), "\n")

stopifnot(identical(colnames(seu_all[["RNA"]]), colnames(seu_all[["ADT"]])))

set.seed(42)
seu_covid <- preprocess(
  object = seu_all,
  remove_doublets = TRUE,
  low_qc_cell_removal = TRUE,
  remove_empty_droplets = FALSE,
  resolution = 0.8,
  do_clustering = FALSE,
  seed = 42,
  return_plots = FALSE,
  print_plots = FALSE,
  species = "Hs"
)

cat("After doublet removal:", ncol(seu_covid), "\n")
print(table(seu_covid$subject_true))


stopifnot("ADT" %in% Assays(seu_covid))
stopifnot(identical(colnames(seu_covid[["RNA"]]), colnames(seu_covid[["ADT"]])))

seu_covid <- RunAzimuth(seu_covid, reference="pbmcref")

print(summary(seu_covid$predicted.celltype.l1.score))
print(table(seu_covid$predicted.celltype.l1))
hist(seu_covid$predicted.celltype.l1.score, breaks=50)

table(seu_covid$predicted.celltype.l1.score >= 0.5)
table(seu_covid$predicted.celltype.l1.score >= 0.5, 
      seu_covid$predicted.celltype.l1)

keep <- seu_covid$predicted.celltype.l1.score >= 0.5
seu_covid <- subset(seu_covid, cells=colnames(seu_covid)[keep])

cat("L1 table:\n")
print(table(seu_covid$predicted.celltype.l1))

seu_covid <- NormalizeData(seu_covid, assay="RNA", normalization.method="LogNormalize", scale.factor=1e4)
seu_covid <- FindVariableFeatures(seu_covid, assay="RNA", selection.method="vst", nfeatures=2000, verbose=FALSE)
hvg <- VariableFeatures(seu_covid)
seu_covid <- ScaleData(seu_covid, assay="RNA", features=hvg, verbose=FALSE)

seu_covid <- RunPCA(seu_covid, assay="RNA", features=hvg, npcs=20, verbose=FALSE)
ElbowPlot(seu_covid, ndims = 30)
Zfull <- Embeddings(seu_covid, "pca")
Z <- Zfull[, seq_len(min(20, ncol(Zfull))), drop=FALSE]
X <- t(as.matrix(GetAssayData(seu_covid, assay="RNA", layer="data")[hvg, , drop=FALSE]))
rownames(X) <- colnames(seu_covid)

Y_counts <- GetAssayData(seu_covid, assay="ADT", layer="counts")  # proteins x cells

# Align again (paranoia)
common <- Reduce(intersect, list(rownames(X), rownames(Z), colnames(Y_counts)))
common <- sort(common)

X <- X[common, , drop=FALSE]
Z <- Z[common, , drop=FALSE]
Y_counts <- Y_counts[, common, drop=FALSE]

stopifnot(
  identical(rownames(X), rownames(Z)),
  identical(rownames(X), colnames(Y_counts))
)

set.seed(42)
donors <- sort(unique(seu_covid$subject_true))
donors
# "NS1" "TP7" "TS5"

for (i in seq_along(donors)) {
  
  test_donor   <- donors[i]
  train_donors <- setdiff(donors, test_donor)
  
  train_cells <- colnames(seu_covid)[seu_covid$subject_true %in% train_donors]
  test_cells  <- colnames(seu_covid)[seu_covid$subject_true %in% test_donor]
  
  split_name <- paste0("LODO_", i, "_test_", test_donor)
  out_base_i <- file.path(out_base, split_name)
  out_ctp_i  <- file.path(out_ctp, split_name)
  
  dir.create(out_base_i, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_ctp_i, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(train_cells, file.path(out_base_i, "train_cells.rds"))
  saveRDS(test_cells,  file.path(out_base_i, "test_cells.rds"))
  saveRDS(X,           file.path(out_base_i, "X.rds"))
  saveRDS(Z,           file.path(out_base_i, "Z.rds"))
  saveRDS(Y_counts,    file.path(out_base_i, "adt_mat.rds"))
  saveRDS(seu_covid,   file.path(out_base_i, "seu_final.rds"))
  saveRDS(hvg,   file.path(out_base_i, "hvg_genes.rds"))
  rna_train <- X[train_cells, , drop = FALSE]
  rna_test  <- X[test_cells,  , drop = FALSE]
  
  adt_train <- t(Y_counts[, train_cells, drop = FALSE])
  adt_test  <- t(Y_counts[, test_cells,  drop = FALSE])
  
  stopifnot(identical(rownames(rna_train), rownames(adt_train)))
  stopifnot(identical(rownames(rna_test),  rownames(adt_test)))
  
  write.csv(rna_train, file.path(out_ctp_i, "rna_train.csv"))
  write.csv(rna_test,  file.path(out_ctp_i, "rna_test.csv"))
  write.csv(adt_train, file.path(out_ctp_i, "adt_train.csv"))
  write.csv(adt_test,  file.path(out_ctp_i, "adt_test.csv"))
  
  cat("\nDONE:", split_name, "\n")
  cat("Test donor :", test_donor, "\n")
  cat("Train donors:", paste(train_donors, collapse = ", "), "\n")
  cat("n_train:", length(train_cells), " n_test:", length(test_cells), "\n")
}







