#' scStable: stability-aware analysis with synthetic scRNA-seq samples
#'
#' \pkg{scStable} generates multiple synthetic scRNA-seq samples from a single
#' scRNA-seq sample (the \emph{scRNA-seq reference}). The samples keep the same
#' cells as the reference and differ in gene expression by between-sample
#' variation learned from multi-sample bulk RNA-seq data.
#'
#' Two modes:
#' \itemize{
#'   \item \strong{Reference mode}: a tissue- and condition-matched
#'         multi-sample \emph{bulk RNA-seq reference} is supplied and its
#'         between-sample variation is estimated with \code{\link{fit_bulk}}.
#'   \item \strong{Bulk-reference-free mode}: pre-estimated tissue-specific
#'         parameters (e.g. learned from GTEx) are used instead; see
#'         \code{\link{synthreplicate_from_tissue}}.
#' }
#'
#' Workflow (preprocessing with \code{\link{synthreplicate_prep}}, then):
#' \enumerate{
#'   \item \code{\link{scDesign3_fit}} -- fit a cell-label-free scRNA-seq
#'         generative model (PCs as covariates).
#'   \item \code{\link{fit_bulk}} -- estimate bulk-derived between-sample
#'         variation (or use pre-estimated parameters).
#'   \item \code{\link{synthreplicate_gen_bulk}} and
#'         \code{\link{synthreplicate_gen_sc}} -- map the variation onto the
#'         scRNA-seq model, with scale factor \eqn{\lambda}
#'         (\code{scaling_factor}).
#'   \item \code{\link{synthreplicate_gen_sc}} -- generate the synthetic
#'         scRNA-seq samples.
#' }
#'
#' @keywords internal
"_PACKAGE"
