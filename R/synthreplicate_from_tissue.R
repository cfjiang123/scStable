###############################################################################
# scStable - bulk-reference-free mode (end-to-end wrapper)
###############################################################################

#' Run scStable in the bulk-reference-free mode
#'
#' When no tissue- and condition-matched bulk RNA-seq reference is available,
#' scStable uses \strong{pre-estimated, tissue-specific parameters}
#' \eqn{(\hat\mu, \hat\Sigma)} of the Gaussian bulk model, learned from a large
#' external resource such as GTEx. \code{synthreplicate_from_tissue} runs the
#' complete workflow (Steps 1-4) in this mode and writes the synthetic
#' scRNA-seq samples to \code{replicate_dir}.
#'
#' @details
#' The between-sample variation can be supplied in either of two ways:
#' \describe{
#'   \item{\code{bulk_params} (recommended)}{A pre-estimated parameter set: the
#'     list returned by \code{\link{fit_bulk}} (elements \code{mu}, \code{cov},
#'     \code{optimal_c}, \code{bulk_mean}, \code{n_samples}), or the path to an
#'     \code{.rds} file containing it. \code{cov} may be a full matrix or a
#'     vector of gene variances (diagonal model). Parameters are subset to the
#'     highly variable genes of \code{scRNA_matrix}.}
#'   \item{\code{tissue_name} + \code{gtex_data_dir}}{A per-tissue
#'     \code{SummarizedExperiment} at \code{<gtex_data_dir>/<tissue_name>.RDS}
#'     (raw counts in \code{assay()}, gene symbols in
#'     \code{rowData()$Description}); the parameters are estimated from it with
#'     \code{\link{fit_bulk}}.}
#' }
#' In reference mode (a matched bulk RNA-seq reference is available), call the
#' step functions directly; see the package vignette.
#'
#' @param tissue_name Tissue identifier used to locate
#'   \code{<gtex_data_dir>/<tissue_name>.RDS}. Ignored if \code{bulk_params} is
#'   given.
#' @param scRNA_matrix Single-sample scRNA-seq count matrix (genes x cells)
#'   with gene names as \code{rownames}.
#' @param gtex_data_dir Directory with per-tissue \code{.RDS} files. Ignored if
#'   \code{bulk_params} is given.
#' @param save_dir Directory for the fitted scDesign3 objects.
#' @param replicate_dir Directory for the synthetic scRNA-seq samples.
#' @param Cell_label Optional cell-level covariate \code{DataFrame}.
#' @param number.pc,number.gene Passed to \code{\link{synthreplicate_prep}}.
#' @param n_cores_bulk Passed to \code{\link{fit_bulk}}.
#' @param use.option,celltype_col,mu_formula,n_cores_marginal,family_copula,parallel_para
#'   Passed to \code{\link{scDesign3_fit}}.
#' @param use.cor,min.eig,number.replicate Passed to
#'   \code{\link{synthreplicate_gen_bulk}}. \code{use.cor = NULL} uses the full
#'   covariance for >= 30 bulk samples and the diagonal otherwise.
#' @param match.option,scaling_factor Passed to
#'   \code{\link{synthreplicate_gen_sc}}; \code{scaling_factor} is the scale
#'   factor \eqn{\lambda}.
#' @param bulk_params Pre-estimated bulk parameters (list or \code{.rds} path);
#'   see Details.
#'
#' @return (Invisibly) a list with the outputs of each step (\code{ouput1}:
#'   preprocessing, \code{ouput2}: bulk parameters used, \code{ouput3}:
#'   scDesign3 fit, \code{ouput4}: synthetic bulk, \code{ouput5}: paths of the
#'   synthetic samples), \code{tissue_name} and \code{bulkRNA_matrix} (the
#'   loaded GTEx matrix, or \code{NULL} when \code{bulk_params} is used).
#'
#' @importFrom SummarizedExperiment assay rowData
#' @export
synthreplicate_from_tissue <- function(
    tissue_name        = NULL,
    scRNA_matrix,
    gtex_data_dir      = NULL,
    save_dir,
    replicate_dir,
    Cell_label         = NULL,

    # synthreplicate_prep parameters
    number.pc          = 10,
    number.gene        = 1500,

    # fit_bulk parameters
    n_cores_bulk       = 20,

    # scDesign3_fit parameters
    use.option         = 2,
    celltype_col       = "cell_type",
    mu_formula         = "cell_type",
    n_cores_marginal   = 10,
    family_copula      = "nb",
    parallel_para      = "pbmcmapply",

    # synthreplicate_gen_bulk parameters
    use.cor            = NULL,
    min.eig            = 1,
    number.replicate   = 100,

    # synthreplicate_gen_sc parameters
    match.option       = 2,
    scaling_factor     = 1,

    # pre-estimated bulk parameters
    bulk_params        = NULL
) {
  if (is.null(scRNA_matrix) || is.null(rownames(scRNA_matrix))) {
    stop("scRNA_matrix must be provided with gene names as rownames")
  }
  bulkRNA_matrix <- NULL

  if (!is.null(bulk_params)) {
    # -- A. pre-estimated parameters --
    if (is.character(bulk_params)) bulk_params <- readRDS(bulk_params)
    need <- c("mu", "cov", "optimal_c", "bulk_mean", "n_samples")
    miss <- setdiff(need, names(bulk_params))
    if (length(miss)) stop("bulk_params is missing: ", paste(miss, collapse = ", "))

    genes <- intersect(rownames(scRNA_matrix), names(bulk_params$mu))
    if (length(genes) == 0) stop("No genes shared between scRNA_matrix and bulk_params")
    message(sprintf("Using pre-estimated bulk parameters (%d shared genes)", length(genes)))

    ouput1 <- synthreplicate_prep(
      bulkRNA_matrix = NULL,
      scRNA_matrix   = scRNA_matrix[genes, , drop = FALSE],
      number.pc      = number.pc,
      number.gene    = number.gene
    )
    hvg <- rownames(ouput1$sc)
    cov_h <- if (is.matrix(bulk_params$cov)) bulk_params$cov[hvg, hvg, drop = FALSE] else bulk_params$cov[hvg]
    ouput2 <- list(mu        = bulk_params$mu[hvg],
                   cov       = cov_h,
                   optimal_c = bulk_params$optimal_c[hvg],
                   bulk_mean = bulk_params$bulk_mean[hvg],
                   n_samples = bulk_params$n_samples)
    if (is.null(use.cor) && !is.matrix(cov_h)) use.cor <- 2
  } else {
    # -- B. estimate parameters from a public tissue-matched bulk reference --
    if (is.null(tissue_name) || is.null(gtex_data_dir)) {
      stop("Provide bulk_params, or both tissue_name and gtex_data_dir")
    }
    bulk_file_path <- file.path(gtex_data_dir, paste0(tissue_name, ".RDS"))
    if (!file.exists(bulk_file_path)) {
      bulk_file_path <- file.path(gtex_data_dir, paste0(tolower(tissue_name), ".RDS"))
    }
    if (!file.exists(bulk_file_path)) {
      available <- gsub("\\.RDS$", "", list.files(gtex_data_dir, pattern = "\\.RDS$",
                                                  ignore.case = TRUE), ignore.case = TRUE)
      stop(sprintf("Bulk file not found for '%s' in %s. Available tissues: %s",
                   tissue_name, gtex_data_dir, paste(available, collapse = ", ")))
    }
    message(sprintf("Estimating bulk parameters from %s", bulk_file_path))
    bulkRNA <- readRDS(bulk_file_path)
    bulkRNA_matrix <- assay(bulkRNA)
    rownames(bulkRNA_matrix) <- rowData(bulkRNA)$Description

    ouput1 <- synthreplicate_prep(
      bulkRNA_matrix = bulkRNA_matrix,
      scRNA_matrix   = scRNA_matrix,
      number.pc      = number.pc,
      number.gene    = number.gene
    )
    ouput2 <- fit_bulk(bulkRNA_matrix = ouput1$bulk, n_cores = n_cores_bulk)
  }
  lib_size <- sum(ouput2$bulk_mean)

  message("Step 1: fitting the scRNA-seq generative model")
  ouput3 <- scDesign3_fit(
    scRNA_matrix     = ouput1$sc,
    top_pcs          = ouput1$pca,
    Cell_label       = Cell_label,
    save_dir         = save_dir,
    use.option       = use.option,
    celltype_col     = celltype_col,
    mu_formula       = mu_formula,
    n_cores_marginal = n_cores_marginal,
    family_copula    = family_copula,
    parallel_para    = parallel_para
  )

  message("Step 3: drawing synthetic bulk samples")
  ouput4 <- synthreplicate_gen_bulk(
    mu               = ouput2$mu,
    cov              = ouput2$cov,
    optimal_c        = ouput2$optimal_c,
    scRNA_matrix     = ouput1$sc,
    use.cor          = use.cor,
    number.replicate = number.replicate,
    min.eig          = min.eig,
    lib_size         = lib_size,
    n_samples        = ouput2$n_samples
  )

  message("Steps 3-4: generating synthetic scRNA-seq samples")
  ouput5 <- synthreplicate_gen_sc(
    bulk_synth     = ouput4$sampling_bulk,
    mu             = ouput2$mu,
    d              = ouput4$d,
    optimal_c      = ouput2$optimal_c,
    scRNA_matrix   = ouput1$sc,
    scRNA_list     = ouput3,
    match.option   = match.option,
    scaling_factor = scaling_factor,
    save.dir       = replicate_dir,
    lib_size       = lib_size
  )
  message(sprintf("%d synthetic samples written to %s", number.replicate, replicate_dir))

  invisible(list(
    ouput1 = ouput1,
    ouput2 = ouput2,
    ouput3 = ouput3,
    ouput4 = ouput4,
    ouput5 = ouput5,
    tissue_name = tissue_name,
    bulkRNA_matrix = bulkRNA_matrix
  ))
}
