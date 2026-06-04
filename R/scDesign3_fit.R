###############################################################################
# scStable - step 3: fit scDesign3 for scRNA-seq
# input : filtered scRNA-seq matrix (with names), top PCs (with cell names),
#         Cell_label (result from step 1)
# output: scDesign3 output
###############################################################################

#' Fit a single-cell generative model for scStable
#'
#' \code{scDesign3_fit} wraps the \pkg{scDesign3} modelling pipeline
#' (\code{construct_data}, \code{fit_marginal}, \code{fit_copula},
#' \code{extract_para}) used by \pkg{scStable} to learn a generative model of
#' the reference single-cell data. Principal components from
#' \code{\link{synthreplicate_prep}} can be used as continuous covariates. All
#' intermediate objects are saved to \code{save_dir} and also returned.
#'
#' @param scRNA_matrix Single-cell count matrix (genes x cells).
#' @param top_pcs Matrix of single-cell PC embeddings (cells x PCs).
#' @param Cell_label Optional \code{DataFrame} of cell-level covariates (used
#'   when \code{use.option = 1}).
#' @param save_dir Directory in which to write the fitted scDesign3 objects.
#' @param use.option Integer (1, 2 or 3) selecting how covariates / formulas are
#'   constructed; see Details.
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
    use.option = 1,

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

  # 2. default Cell_label
  if (use.option == 1) {
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

  if (use.option == 3) {
    Cell_label <- S4Vectors::DataFrame(
      cell_type = rep(1, ncol(scRNA_matrix)),
      row.names = colnames(scRNA_matrix)
    )

    other_covariates <- paste0("pc", seq_len(ncol(top_pcs)))
    mu_formula <- paste(other_covariates, collapse = " + ")
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
