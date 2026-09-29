###############################################################################
# scStable - Step 1: model the single-sample scRNA-seq reference
###############################################################################

#' Step 1: fit the scRNA-seq generative model
#'
#' \code{scDesign3_fit} wraps the \pkg{scDesign3} pipeline
#' (\code{construct_data}, \code{fit_marginal}, \code{fit_copula},
#' \code{extract_para}) to fit a negative-binomial marginal model with a
#' Gaussian copula to the scRNA-seq reference. By default
#' (\code{use.option = 2}) the model is cell-label-free: the top PCs from
#' \code{\link{synthreplicate_prep}} are the cell covariates
#' (\code{mu ~ pc1 + ... + pck}). All fitted objects are saved to
#' \code{save_dir} and also returned.
#'
#' @param scRNA_matrix Single-cell count matrix (genes x cells).
#' @param top_pcs Matrix of single-cell PC embeddings (cells x PCs).
#' @param Cell_label Optional \code{DataFrame} of cell-level covariates (used
#'   when \code{use.option = 1}).
#' @param save_dir Directory in which to write the fitted scDesign3 objects.
#' @param use.option \code{2} (default): cell-label-free model using the PCs as
#'   covariates. \code{1}: user-defined covariates from \code{Cell_label} and
#'   \code{mu_formula} (PCs are added as columns \code{pc1}, \code{pc2}, ...).
#' @param assay_use,celltype_col,pseudotime_col,spatial_col,other_covariates,corr_by
#'   Arguments forwarded to \code{scDesign3::construct_data}.
#' @param predictor,mu_formula,sigma_formula,family_marginal,n_cores_marginal,usebam,parallel_marginal
#'   Arguments forwarded to \code{scDesign3::fit_marginal}.
#' @param family_copula,copula,n_cores_copula,parallel_copula
#'   Arguments forwarded to \code{scDesign3::fit_copula}.
#' @param n_cores_para,family_para,parallel_para
#'   Arguments forwarded to \code{scDesign3::extract_para}.
#'
#' @return (Invisibly) a list with the constructed \code{SingleCellExperiment}
#'   (\code{scRNA_sce_pc}), the constructed data (\code{scRNA_data_pc}), the
#'   fitted marginals (\code{scRNA_marginal_pc}), the fitted copula
#'   (\code{scRNA_copula_pc}) and the extracted parameters
#'   (\code{scRNA_para_pc}).
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors DataFrame
#' @importFrom scDesign3 construct_data fit_marginal fit_copula extract_para
#' @export
scDesign3_fit <- function(
    scRNA_matrix,
    top_pcs,
    Cell_label           = NULL,
    save_dir,
    use.option = 2,

    ## construct_data args
    assay_use            = "counts",
    celltype_col         = "cell_type",
    pseudotime_col       = NULL,
    spatial_col          = NULL,
    other_covariates     = NULL,
    corr_by              = "1",

    ## fit_marginal args
    predictor            = "gene",
    mu_formula           = "cell_type",
    sigma_formula        = 1,
    family_marginal      = "nb",
    n_cores_marginal     = 10,
    usebam               = FALSE,
    parallel_marginal    = "pbmcmapply",

    ## fit_copula args
    family_copula        = "nb",
    copula               = "gaussian",
    n_cores_copula       = 10,
    parallel_copula      = "pbmcmapply",

    ## extract_para args
    n_cores_para         = 10,
    family_para          = "nb",
    parallel_para        = "pbmcmapply"
) {
  # 1. ensure save_dir exists
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # 2. build the SingleCellExperiment with PCs as cell covariates
  if (!use.option %in% c(1, 2)) stop("use.option must be 1 or 2")
  if (use.option == 1) {
    if (is.null(Cell_label)) stop("use.option = 1 requires Cell_label")
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays  = list(counts = scRNA_matrix),
      colData = Cell_label
    )
    for (i in seq_len(ncol(top_pcs))) {
      sce[[paste0("pc", i)]] <- top_pcs[, i]
    }
  }
  if (use.option == 2) {
    Cell_label <- S4Vectors::DataFrame(
      cell_type = rep(1, ncol(scRNA_matrix)),
      row.names = colnames(scRNA_matrix)
    )

    other_covariates <- paste0("pc", seq_len(ncol(top_pcs)))
    mu_formula <- paste(other_covariates, collapse = " + ")

    sce <- SingleCellExperiment::SingleCellExperiment(
      assays  = list(counts = scRNA_matrix),
      colData = Cell_label
    )
    for (i in seq_len(ncol(top_pcs))) {
      sce[[paste0("pc", i)]] <- top_pcs[, i]
    }
  }


  # 4. construct data
  data_pc <- scDesign3::construct_data(
    sce               = sce,
    assay_use         = assay_use,
    celltype          = celltype_col,
    pseudotime        = pseudotime_col,
    spatial           = spatial_col,
    other_covariates  = other_covariates,
    corr_by           = corr_by
  )

  # 5. fit marginal
  message("Fitting marginal...")
  marginal_pc <- scDesign3::fit_marginal(
    data            = data_pc,
    predictor       = predictor,
    mu_formula      = mu_formula,
    sigma_formula   = sigma_formula,
    family_use      = family_marginal,
    n_cores         = n_cores_marginal,
    usebam          = usebam,
    parallelization = parallel_marginal
  )

  # 6. fit copula
  message("Fitting copula...")
  copula_pc <- scDesign3::fit_copula(
    sce             = sce,
    assay_use       = assay_use,
    marginal_list   = marginal_pc,
    family_use      = family_copula,
    copula          = copula,
    n_cores         = n_cores_copula,
    input_data      = data_pc$dat,
    parallelization = parallel_copula
  )

  # 7. extract parameters
  message("Extracting parameters...")
  para_pc <- scDesign3::extract_para(
    sce             = sce,
    marginal_list   = marginal_pc,
    n_cores         = n_cores_para,
    family_use      = family_para,
    new_covariate   = data_pc$newCovariate,
    data            = data_pc$dat,
    parallelization = parallel_para
  )

  # 8. save objects
  saveRDS(sce,         file = file.path(save_dir, "scRNA_sce_pc.rds"))
  saveRDS(data_pc,     file = file.path(save_dir, "scRNA_data_pc.rds"))
  saveRDS(marginal_pc, file = file.path(save_dir, "scRNA_marginal_pc.rds"))
  saveRDS(copula_pc,   file = file.path(save_dir, "scRNA_copula_pc.rds"))
  saveRDS(para_pc,     file = file.path(save_dir, "scRNA_para_pc.rds"))

  # 9. return list
  invisible(list(
    scRNA_sce_pc      = sce,
    scRNA_data_pc     = data_pc,
    scRNA_marginal_pc = marginal_pc,
    scRNA_copula_pc   = copula_pc,
    scRNA_para_pc     = para_pc
  ))
}
