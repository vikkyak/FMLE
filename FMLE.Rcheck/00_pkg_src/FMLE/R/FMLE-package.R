#' FMLE: Fuzzy Mixture of Linear Experts for Single-Cell Protein Prediction
#'
#' The FMLE package implements fuzzy mixture-of-experts models tailored to
#' predicting surface protein abundance from single-cell RNA data, such as
#' CITE-seq experiments. A fuzzy c-means gating network in a low-dimensional
#' representation (e.g. principal components, latent embeddings, spatial
#' coordinates) softly assigns each cell to multiple linear experts that act
#' on high-dimensional gene expression features.
#'
#' The package provides:
#'
#' - \code{\link{fmle_train}} and \code{\link{fmle_predict}} for single-protein
#'   (single-task) regression.
#' - \code{\link{fmle_train_mt}} and \code{\link{fmle_predict_mt}} for
#'   multi-protein (multi-task) prediction with shared fuzzy gates and
#'   protein-specific experts.
#' - Optional L1 penalization for expert coefficients (lasso plus post-ridge),
#'   combined with a small ridge term for numerical stability.
#' - Prediction uncertainty with separate aleatoric and parameter components.
#'
#' While FMLE can be applied to general high-dimensional regression problems,
#' its primary motivation is accurate and interpretable prediction of protein
#' abundance at single-cell resolution, across donors, batches and antibody
#' panels.
#'
#' @keywords internal
"_PACKAGE"
