preprocess_fmle <- function (object,
                                 remove_doublets = TRUE,
                                 low_qc_cell_removal = TRUE,
                                 anno_level = 2,
                                 samples = NULL,
                                 integrate_data = FALSE,
                                 remove_empty_droplets = FALSE,
                                 lower = 100,
                                 FDR = 0.01,
                                 resolution = 0.8,
                                 seed = 42,
                                 return_plots = FALSE,
                                 print_plots = TRUE,
                                 species = "Hs",
                                 min.features = NULL,
                                 verbose = FALSE) {
  
  set.seed(seed)
  plot_list <- list()
  Seurat::DefaultAssay(object) <- "RNA"
  
  # empty droplets (optional)
  if (remove_empty_droplets) {
    object <- empty_drops(object = object, lower = lower,
                          FDR = FDR, samples = samples, seed = seed)
  }
  
  # mito_percent (not used here but fine)
  if (!("mito_percent" %in% names(object@meta.data))) {
    if (species == "Hs") {
      object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
    } else {
      object$mito_percent <- Seurat::PercentageFeatureSet(object, pattern = "^mt-")
    }
  }
  
  # doublet removal
  if (remove_doublets) {
    message("Start remove doublets")
    object <- object %>% remove_doublets(
      samples     = samples,
      print_plots = print_plots,
      seed        = seed,
      verbose     = verbose
    )
    plot_list[["doublets"]] <- object[[2]]
    object <- object[[1]]
  }
  
  # low-quality cell removal (MAD-based)
  if (low_qc_cell_removal) {
    message("Start low quality cell removal")
    object <- object %>% mad_filtering(
      samples     = samples,
      print_plots = print_plots,
      seed        = seed,
      min.features = min.features,
      verbose      = verbose
    )
    plot_list[["low_qc_cells"]] <- object[[2]]
    object <- object[[1]]
  }
  
  # optional integration
  if (integrate_data) {
    message("Start integrate data")
    object <- integrate_samples(object, samples = samples,
                                seed = seed, verbose = verbose)
  }
  
  message("Start clustering data")
  object <- cluster_data(object, resolution = resolution, seed = seed)
  Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]
  
  if (return_plots) {
    return(list(object = object, plots = plot_list))
  } else {
    return(object)
  }
}
