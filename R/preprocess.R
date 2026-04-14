#' Preprocess a Seurat object for FMLE workflows
#'
#' Performs optional empty-droplet removal, mitochondrial percentage calculation,
#' doublet removal, low-quality cell filtering, optional clustering.
#'
#' @param object A Seurat object.
#' @param remove_doublets Logical; whether to remove doublets.
#' @param low_qc_cell_removal Logical; whether to remove low-quality cells.
#' @param samples Optional metadata column used for sample-wise processing.
#' @param remove_empty_droplets Logical; whether to apply empty droplet filtering.
#' @param lower Lower count threshold for empty droplet filtering.
#' @param FDR FDR threshold for empty droplet filtering.
#' @param resolution Clustering resolution.
#' @param seed Random seed.
#' @param return_plots Logical; whether to return QC plots.
#' @param print_plots Logical; whether to print plots.
#' @param species Species string, e.g. "Hs".
#' @param min.features Optional minimum feature threshold in QC filtering.
#' @param verbose Logical.
#'
#' @return A Seurat object, or a list with object and plots if `return_plots=TRUE`.
#' @export

preprocess <- function(object,
                       remove_doublets = TRUE,
                       low_qc_cell_removal = TRUE,
                       anno_level = 2,
                       samples = NULL,
                       integrate_data = FALSE,
                       remove_empty_droplets = FALSE,
                       do_clustering = FALSE,
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
    object <- .empty_drops_fmle(object = object, lower = lower,
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
    object <- .remove_doublets_fmle(
      object      = object,
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
    object <- .mad_filtering_fmle(
      object      = object,
      samples     = samples,
      print_plots = print_plots,
      seed        = seed,
      min.features = min.features,
      verbose      = verbose
    )
    plot_list[["low_qc_cells"]] <- object[[2]]
    object <- object[[1]]
  }
  

  if (do_clustering) {
    message("Start clustering data")
    object <- .cluster_data_fmle(
      object = object,
      resolution = resolution,
      seed = seed,
      verbose = verbose
    )
    Seurat::Idents(object) <- object@meta.data[["seurat_clusters"]]
  }
  
  if (return_plots) {
    return(list(object = object, plots = plot_list))
  } else {
    return(object)
  }
}
