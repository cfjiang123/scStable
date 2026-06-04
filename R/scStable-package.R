#' scStable: Synthetic Single-Cell RNA-seq Replicate Generation
#'
#' The \pkg{scStable} package implements the scStable method for generating
#' synthetic single-cell RNA-seq (scRNA-seq) replicates anchored to bulk
#' RNA-seq data. The workflow proceeds in five steps:
#'
#' \enumerate{
#'   \item \code{\link{synthreplicate_prep}} -- harmonise bulk and single-cell
#'         matrices to a shared, informative gene set and (optionally) compute
#'         single-cell principal components.
#'   \item \code{\link{fit_bulk}} -- fit a per-gene log-normal/multivariate
#'         model to the bulk RNA-seq matrix.
#'   \item \code{\link{scDesign3_fit}} -- fit a single-cell generative model
#'         with \pkg{scDesign3}.
#'   \item \code{\link{synthreplicate_gen_bulk}} -- sample synthetic bulk
#'         replicates from the fitted bulk model.
#'   \item \code{\link{synthreplicate_gen_sc}} -- propagate the synthetic bulk
#'         variation into newly simulated single-cell count matrices.
#' }
#'
#' The convenience wrapper \code{\link{synthreplicate_from_tissue}} runs the
#' full scStable pipeline end-to-end for a single tissue.
#'
#' @keywords internal
"_PACKAGE"
