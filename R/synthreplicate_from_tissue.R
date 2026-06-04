###############################################################################
# scStable - end-to-end wrapper: run the full pipeline for one tissue
###############################################################################

#' Run the full scStable pipeline for a single tissue
#'
#' \code{synthreplicate_from_tissue} is a convenience wrapper that executes the
#' complete \pkg{scStable} workflow end-to-end: it loads a GTEx bulk RNA-seq
#' object for the requested tissue, preprocesses it together with a reference
#' single-cell matrix (\code{\link{synthreplicate_prep}}), fits the bulk
#' distribution (\code{\link{fit_bulk}}), fits the single-cell model
#' (\code{\link{scDesign3_fit}}), samples synthetic bulk replicates
#' (\code{\link{synthreplicate_gen_bulk}}) and finally generates synthetic
#' single-cell replicates (\code{\link{synthreplicate_gen_sc}}).
#'
#' @param tissue_name Character; tissue identifier used to locate
#'   \code{<gtex_data_dir>/<tissue_name>.RDS}.
#' @param scRNA_matrix Reference single-cell count matrix (genes x cells) with
#'   gene names as \code{rownames}.
#' @param gtex_data_dir Directory containing per-tissue GTEx \code{.RDS} files.
#' @param save_dir Directory for fitted scDesign3 objects.
#' @param replicate_dir Directory for the generated synthetic replicate files.
#' @param Cell_label Optional cell-level covariate \code{DataFrame}.
#' @param number.pc,number.gene Passed to \code{\link{synthreplicate_prep}}.
#' @param n_cores_bulk Passed to \code{\link{fit_bulk}}.
#' @param use.option,celltype_col,mu_formula,n_cores_marginal,family_copula,parallel_para
#'   Passed to \code{\link{scDesign3_fit}}.
#' @param use.cor,min.eig,number.replicate Passed to
#'   \code{\link{synthreplicate_gen_bulk}}.
#' @param match.option,scaling_factor Passed to
#'   \code{\link{synthreplicate_gen_sc}}.
#'
#' @return (Invisibly) a list collecting the outputs of each pipeline step
#'   (\code{ouput1}..\code{ouput5}), the \code{tissue_name} and the loaded
#'   \code{bulkRNA_matrix}.
#'
#' @importFrom SummarizedExperiment assay rowData
#' @export
synthreplicate_from_tissue <- function(
    tissue_name,
    scRNA_matrix,
    gtex_data_dir      = "/home/chengfeng/scRobust/data/GTEx/tissue_data/",
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
    use.cor            = 3,
    min.eig            = 1,
    number.replicate   = 100,

    # synthreplicate_gen_sc parameters
    match.option       = 2,
    scaling_factor     = 1
) {

  # -- 0. Validate inputs --
  if (is.null(tissue_name) || !is.character(tissue_name) || nchar(tissue_name) == 0) {
    stop("tissue_name must be a non-empty character string")
  }

  if (is.null(scRNA_matrix)) {
    stop("scRNA_matrix must be provided")
  }

  if (is.null(rownames(scRNA_matrix))) {
    stop("scRNA_matrix must have gene names as rownames")
  }

  # -- 1. Construct file path and load bulk RNA-seq data --
  message("===============================================================")
  message(sprintf("Loading bulk RNA-seq data for tissue: %s", tissue_name))
  message("===============================================================")

  # Ensure directory path ends with "/"
  if (!endsWith(gtex_data_dir, "/")) {
    gtex_data_dir <- paste0(gtex_data_dir, "/")
  }

  # Construct full file path
  bulk_file_path <- paste0(gtex_data_dir, tissue_name, ".RDS")

  # Check if file exists
  if (!file.exists(bulk_file_path)) {
    # Try lowercase
    bulk_file_path_lower <- paste0(gtex_data_dir, tolower(tissue_name), ".RDS")
    if (file.exists(bulk_file_path_lower)) {
      bulk_file_path <- bulk_file_path_lower
    } else {
      # List available tissues
      available_files <- list.files(gtex_data_dir, pattern = "\\.RDS$", ignore.case = TRUE)
      available_tissues <- gsub("\\.RDS$", "", available_files, ignore.case = TRUE)
      stop(sprintf(
        "Bulk RNA-seq file not found: %s\nAvailable tissues: %s",
        bulk_file_path,
        paste(available_tissues, collapse = ", ")
      ))
    }
  }

  # Load bulk RNA-seq data
  bulkRNA <- readRDS(bulk_file_path)
  bulkRNA_matrix <- assay(bulkRNA)
  rownames(bulkRNA_matrix) <- rowData(bulkRNA)$Description


  message(sprintf("Loaded bulk RNA-seq matrix: %d genes x %d samples",
                  nrow(bulkRNA_matrix), ncol(bulkRNA_matrix)))

  # -- 2. Run synthreplicate_prep --
  message("\n===============================================================")
  message("Step 1: Preprocessing bulk and single-cell data")
  message("===============================================================")

  ouput1 <- synthreplicate_prep(
    bulkRNA_matrix = bulkRNA_matrix,
    scRNA_matrix   = scRNA_matrix,
    number.pc      = number.pc,
    number.gene    = number.gene
  )

  message(sprintf("Selected %d highly variable genes", length(ouput1$hvg)))
  message(sprintf("Filtered bulk matrix: %d genes x %d samples",
                  nrow(ouput1$bulk), ncol(ouput1$bulk)))
  message(sprintf("Filtered sc matrix: %d genes x %d cells",
                  nrow(ouput1$sc), ncol(ouput1$sc)))

  # -- 3. Fit bulk distribution --
  message("\n===============================================================")
  message("Step 2: Fitting bulk RNA-seq distribution")
  message("===============================================================")

  ouput2 <- fit_bulk(
    bulkRNA_matrix = ouput1$bulk,
    n_cores        = n_cores_bulk
  )

  message("Bulk distribution fitting completed")

  # -- 4. Fit scDesign3 --
  message("\n===============================================================")
  message("Step 3: Fitting scDesign3 model for single-cell data")
  message("===============================================================")

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

  message("scDesign3 fitting completed")

  # -- 5. Generate synthetic bulk samples --
  message("\n===============================================================")
  message("Step 4: Generating synthetic bulk samples")
  message("===============================================================")

  ouput4 <- synthreplicate_gen_bulk(
    bulkRNA_matrix   = ouput1$bulk,
    mu               = ouput2$mu,
    cov              = ouput2$cov,
    optimal_c        = ouput2$optimal_c,
    scRNA_matrix     = ouput1$sc,
    use.cor          = use.cor,
    number.replicate = number.replicate,
    min.eig          = min.eig
  )

  message(sprintf("Generated %d synthetic bulk replicates", number.replicate))

  # -- 6. Generate synthetic single-cell data --
  message("\n===============================================================")
  message("Step 5: Generating synthetic single-cell replicates")
  message("===============================================================")

  ouput5 <- synthreplicate_gen_sc(
    bulkRNA_matrix = ouput1$bulk,
    bulk_synth     = ouput4$sampling_bulk,
    mu             = ouput2$mu,
    d              = ouput4$d,
    optimal_c      = ouput2$optimal_c,
    scRNA_matrix   = ouput1$sc,
    scRNA_list     = ouput3,
    match.option   = match.option,
    scaling_factor = scaling_factor,
    save.dir       = replicate_dir
  )

  message(sprintf("Synthetic single-cell replicates saved to: %s", replicate_dir))

  # -- 7. Summary and return --
  message("\n===============================================================")
  message("scStable pipeline completed successfully!")
  message("===============================================================")
  message(sprintf("Tissue: %s", tissue_name))
  message(sprintf("Number of replicates: %d", number.replicate))
  message(sprintf("Synthetic sc files: %s/replicate*.csv", replicate_dir))

  # Return results
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
