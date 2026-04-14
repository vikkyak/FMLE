#' @keywords internal
#' @importFrom magrittr %>%
.remove_doublets_fmle <- function(object, samples = NULL, remove_cells = TRUE,
                            seed = 42, print_plots = TRUE, verbose = FALSE, ...) {
  set.seed(seed = seed)
  
  if (length(SeuratObject::Layers(object)) > 1) {
    object_joined <- SeuratObject::JoinLayers(object)
  } else {
    object_joined <- object
  }
  
  sce <- Seurat::as.SingleCellExperiment(object_joined)
  sce <- scDblFinder::scDblFinder(sce, samples = samples, verbose = verbose, ...)
  
  if (is.null(samples)) {
    metadata <- sce@colData@listData %>%
      as.data.frame() %>%
      dplyr::select("scDblFinder.class")
    
    p <- ggplot2::ggplot(
      metadata,
      ggplot2::aes(x = "", fill = scDblFinder.class, label = ggplot2::after_stat(count))
    ) +
      ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") +
      ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = 1.2) +
      ggplot2::labs(title = "Number of doublets by sample", fill = "Type") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = -0.5, hjust = 1)) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("# cells")
    
    if (print_plots) base::print(p)
  } else {
    metadata <- sce@colData@listData %>%
      as.data.frame() %>%
      dplyr::select("scDblFinder.sample", "scDblFinder.class")
    
    p <- ggplot2::ggplot(
      metadata,
      ggplot2::aes(x = scDblFinder.sample, fill = scDblFinder.class, label = ggplot2::after_stat(count))
    ) +
      ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") +
      ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = 1.2) +
      ggplot2::labs(title = "Number of doublets by sample", fill = "Type") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = -0.5, hjust = 1)) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("# cells")
    
    if (print_plots) base::print(p)
  }
  
  object[["scDblFinder.class"]] <- sce@colData@listData[["scDblFinder.class"]]
  
  if (remove_cells) {
    object <- base::subset(object, subset = scDblFinder.class == "singlet")
    object@meta.data[["scDblFinder.class"]] <- NULL
  }
  
  return(list(object = object, plots = p))
}

.mad_filtering_fmle <- function(object, samples = NULL, nmads = 3, type = "both",
                          mttype = "higher", remove_cells = TRUE, print_plots = TRUE,
                          seed = 42, min.features = NULL, verbose = FALSE, ...) {
  set.seed(seed)
  
  if (is.null(samples)) {
    batch <- NULL
  } else {
    batch <- object@meta.data %>% dplyr::pull(samples)
  }
  
  nCount_ol <- scater::isOutlier(object@meta.data$nCount_RNA,
                                 nmads = nmads, type = type, log = TRUE, batch = batch)
  nFeature_ol <- scater::isOutlier(object@meta.data$nFeature_RNA,
                                   nmads = nmads, type = type, log = TRUE, batch = batch)
  
  if (!is.null(min.features)) {
    TF <- (object@meta.data$nFeature_RNA < min.features)
    nFeature_ol <- (nFeature_ol | TF)
  }
  
  if (!is.null(object@meta.data$mito_percent)) {
    pMito_ol <- scater::isOutlier(object@meta.data$mito_percent,
                                  nmads = nmads, type = mttype, log = FALSE, batch = batch)
    object@meta.data$mad_filtered <- (nCount_ol | nFeature_ol | pMito_ol)
  } else {
    object@meta.data$mad_filtered <- (nCount_ol | nFeature_ol)
  }
  
  metadata <- object@meta.data %>%
    dplyr::select(tidyselect::any_of(c("nCount_RNA", "nFeature_RNA", "mito_percent", "mad_filtered"))) %>%
    dplyr::rename(Filtered = "mad_filtered")
  
  p <- ggplot2::ggplot(metadata, ggplot2::aes(x = nCount_RNA, y = nFeature_RNA, color = Filtered)) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c("darkgreen", "darkred")) +
    ggplot2::scale_x_continuous(trans = "log10") +
    ggplot2::scale_y_continuous(trans = "log10")
  
  if (print_plots) base::print(p)
  
  if (is.null(samples)) {
    p <- ggplot2::ggplot(metadata, ggplot2::aes(x = "", fill = Filtered, label = ggplot2::after_stat(count))) +
      ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") +
      ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = 1.2) +
      ggplot2::labs(title = "Number of filtered cells by sample", fill = "Filtered") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = -0.5, hjust = 1)) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("# cells")
    
    if (print_plots) base::print(p)
  } else {
    metadata <- cbind(metadata, data.frame(samples = object@meta.data %>% dplyr::pull(samples)))
    
    p <- ggplot2::ggplot(metadata, ggplot2::aes(x = samples, fill = Filtered, label = ggplot2::after_stat(count))) +
      ggplot2::theme_bw() +
      ggplot2::geom_bar(position = "identity", stat = "count") +
      ggplot2::scale_fill_manual(values = pals::kelly()[3:4]) +
      ggplot2::geom_text(stat = "count", vjust = 1.2) +
      ggplot2::labs(title = "Number of quality filtered cells by sample", fill = "Filtered") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = -0.5, hjust = 1)) +
      ggplot2::xlab("Sample") +
      ggplot2::ylab("# cells")
    
    if (print_plots) base::print(p)
  }
  
  if (remove_cells) {
    object <- base::subset(object, subset = mad_filtered == "FALSE")
    object@meta.data[["mad_filtered"]] <- NULL
  }
  
  return(list(object = object, plots = p))
}


.cluster_data_fmle <- function(object, resolution = 0.8, npcs_calculate = 100, npcs = NULL,
                         seed = 42, verbose = FALSE) {
  set.seed(seed)
  default_assay <- Seurat::DefaultAssay(object)
  
  if (default_assay == "integrated") {
    object <- object %>%
      Seurat::ScaleData(verbose = verbose) %>%
      Seurat::RunPCA(assay = default_assay, npcs = npcs_calculate, verbose = verbose)
  } else {
    object <- object %>%
      Seurat::NormalizeData(verbose = verbose) %>%
      Seurat::ScaleData(verbose = verbose) %>%
      Seurat::FindVariableFeatures(verbose = verbose) %>%
      Seurat::RunPCA(assay = default_assay, npcs = npcs_calculate, verbose = verbose)
  }
  
  if (is.null(npcs)) {
    ndims <- ceiling(
      intrinsicDimension::maxLikGlobalDimEst(
        object@reductions[[paste0("pca")]]@cell.embeddings,
        k = 20
      )[["dim.est"]]
    )
  } else {
    ndims <- npcs
  }
  
  print(paste0("Number of used dimensions for clustering: ", ndims))
  
  object <- object %>%
    Seurat::RunUMAP(dims = 1:ndims, verbose = verbose) %>%
    Seurat::FindNeighbors(dims = 1:ndims, reduction = "pca", verbose = verbose) %>%
    Seurat::FindClusters(resolution = resolution, verbose = verbose)
  
  return(object)
}

.empty_drops_fmle <- function(object, lower = 100, FDR = 0.01, samples = NULL, seed = 42, ...) {
  set.seed(seed)
  
  if (is.null(samples)) {
    e.out <- DropletUtils::emptyDrops(object@assays$RNA@counts, lower = lower, ...)
    is.cell <- e.out$FDR <= FDR
    object <- object[, which(is.cell)]
  } else {
    object_split <- Seurat::SplitObject(object, split.by = samples)
    for (i in names(object_split)) {
      e.out <- DropletUtils::emptyDrops(object_split[[i]]@assays$RNA@counts, lower = lower, ...)
      is.cell <- e.out$FDR <= FDR
      object_split[[i]] <- object_split[[i]][, which(is.cell)]
    }
    object <- merge(object_split[[1]], y = object_split[2:length(object_split)])
  }
  
  return(object)
}





