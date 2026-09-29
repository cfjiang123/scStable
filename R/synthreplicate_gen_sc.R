###############################################################################
# scStable - Step 3 (part 2) + Step 4: map variation and generate samples
###############################################################################

#' Steps 3-4: map between-sample variation and generate synthetic scRNA-seq samples
#'
#' For each synthetic bulk profile from \code{\link{synthreplicate_gen_bulk}},
#' \code{synthreplicate_gen_sc} converts it back to the pseudo-bulk scale
#' (removing the offset \eqn{\xi_g} and inverting the log transform), computes
#' gene-wise multiplicative adjustments \eqn{\omega_{gr}^\lambda}, re-weights the
#' scDesign3 mean matrix and simulates a new count matrix with
#' \code{scDesign3::simu_new}. All synthetic samples contain the same cells as
#' the scRNA-seq reference. Each sample is written to
#' \code{<save.dir>/replicate<r>.csv} (tab-separated, genes x cells).
#'
#' @param bulkRNA_matrix Bulk matrix from \code{\link{synthreplicate_prep}}
#'   (reference mode); only its mean library size is used. May be \code{NULL}
#'   if \code{lib_size} is given (bulk-reference-free mode).
#' @param bulk_synth Synthetic bulk matrix (genes x samples),
#'   \code{synthreplicate_gen_bulk()$sampling_bulk}.
#' @param mu Gene-wise bulk mean from \code{\link{fit_bulk}}.
#' @param d Per-gene offset, \code{synthreplicate_gen_bulk()$d}.
#' @param optimal_c Per-gene pseudo-counts from \code{\link{fit_bulk}}.
#' @param scRNA_matrix Single-cell matrix from \code{\link{synthreplicate_prep}}.
#' @param scRNA_list The list returned by \code{\link{scDesign3_fit}}.
#' @param match.option \code{2} (default, used in the manuscript) maps each
#'   synthetic bulk profile to a predicted pseudo-bulk and divides by the
#'   observed pseudo-bulk; \code{1} uses the ratio to the bulk mean directly.
#' @param scaling_factor The scale factor \eqn{\lambda}: \code{0} no injected
#'   variation, \code{1} realistic bulk-derived variation, \code{>1} amplified
#'   variation for stress-testing.
#' @param use.pc Logical / integer; if true, genes whose simulated total is
#'   extremely inflated relative to the reference are rescaled after simulation.
#' @param n_cores Cores passed to \code{scDesign3::simu_new}.
#' @param sc_quantile Reserved; not used.
#' @param save.dir Directory in which synthetic samples are written.
#' @param lib_size Mean bulk library size; defaults to
#'   \code{sum(bulkRNA_matrix) / ncol(bulkRNA_matrix)}.
#'
#' @return Invisibly, the paths of the written \code{replicate<r>.csv} files.
#'
#' @importFrom scDesign3 simu_new
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stats quantile IQR median
#' @importFrom utils write.table
#' @export
synthreplicate_gen_sc <-function(bulkRNA_matrix = NULL, bulk_synth, mu, d, optimal_c, scRNA_matrix, scRNA_list,
                                 match.option = 2, scaling_factor = 1, use.pc = 1, n_cores = 1, sc_quantile = 0.995, save.dir,
                                 lib_size = NULL){
  if (!match.option %in% c(1, 2)) stop("match.option must be 1 or 2")
  if (is.null(lib_size)) {
    if (is.null(bulkRNA_matrix)) stop("Provide either bulkRNA_matrix or lib_size")
    lib_size <- sum(bulkRNA_matrix) / ncol(bulkRNA_matrix)
  }
  if (!identical(rownames(bulk_synth), rownames(scRNA_matrix))) {
    stop("rownames(bulk_synth) must match rownames(scRNA_matrix) in the same order")
  }
  if (!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
  pseudo_bulk = rowSums(scRNA_matrix)
  e_v = quantile(pseudo_bulk, 0.75, na.rm = TRUE) + 1.5 * IQR(pseudo_bulk, na.rm = TRUE)
  names(pseudo_bulk) <- rownames(scRNA_matrix)
  alpha = lib_size / sum(pseudo_bulk) # scaling
  # create mapping vector
  per_gene_mapping = list()
  if(match.option == 1){ ## precise matching
    for (i in 1:ncol(bulk_synth)) {
      per_gene_mapping[[i]] = (exp(bulk_synth[,i] - d) - optimal_c) / (exp(mu - d) - optimal_c)
      e_v_m = quantile(per_gene_mapping[[i]], 0.75, na.rm = TRUE) + 1.5 * IQR(per_gene_mapping[[i]], na.rm = TRUE)
      per_gene_mapping[[i]][per_gene_mapping[[i]]>e_v_m] = e_v_m
      per_gene_mapping[[i]] = per_gene_mapping[[i]] ^ scaling_factor
    }
  }

  if(match.option == 2){
    for (i in 1:ncol(bulk_synth)) {
      v <- bulk_synth[,i]
      pseudo_bulk_hat = (exp(v - d) - optimal_c) / alpha
      pseudo_bulk_hat[pseudo_bulk_hat < 0] = 0
      pseudo_bulk_hat[pseudo_bulk_hat > e_v] = e_v
      per_gene_mapping[[i]] = pseudo_bulk_hat / pseudo_bulk
      e_v_m = quantile(per_gene_mapping[[i]], 0.75, na.rm = TRUE) + 1.5 * IQR(per_gene_mapping[[i]], na.rm = TRUE)
      per_gene_mapping[[i]][per_gene_mapping[[i]]>e_v_m] = e_v_m
      per_gene_mapping[[i]] = per_gene_mapping[[i]] ^ scaling_factor
    }
  }


  # generate new replicate
  Generate_multiple_counts = lapply(1:ncol(bulk_synth), function(i){
    set.seed(i)
    mean_mat = scRNA_list$scRNA_para_pc$mean_mat

    per_gene_mapping[[i]][per_gene_mapping[[i]] == 0] = 1e-6
    mean_mat_weighted <- sweep(mean_mat, 2, per_gene_mapping[[i]], `*`)

    colnames(mean_mat_weighted) = colnames(scRNA_list$scRNA_para_pc$mean_mat)


    newcount <- simu_new(
      sce = scRNA_list$scRNA_sce_pc,
      mean_mat = mean_mat_weighted,
      sigma_mat = scRNA_list$scRNA_para_pc$sigma_mat,
      zero_mat = scRNA_list$scRNA_para_pc$zero_mat,
      quantile_mat = NULL,
      copula_list = scRNA_list$scRNA_copula_pc$copula_list,
      n_cores = n_cores,
      family_use = "nb",
      input_data = scRNA_list$scRNA_data_pc$dat,
      new_covariate = scRNA_list$scRNA_data_pc$newCovariate,
      parallelization = "bpmcmapply",
      BPPARAM = BiocParallel::MulticoreParam(),
      important_feature = scRNA_list$scRNA_copula_pc$important_feature,
      filtered_gene = scRNA_list$scRNA_data_pc$filtered_gene
    )
    if(use.pc){
      rs_new <- rowSums(newcount, na.rm = TRUE)
      rs_sc  <- rowSums(scRNA_matrix,  na.rm = TRUE)

      fc <- rs_new / (rs_sc + 1e-6)
      fc_reverse <- rs_sc / (rs_new + 1e-6)
      fc_reverse = fc_reverse[rs_new != 0]

      genes_fc2_high <- names(fc)[fc > 3 * max(fc_reverse) & rs_new > median(rs_new)]
      newcount[genes_fc2_high, ] = newcount[genes_fc2_high, ] * fc_reverse[genes_fc2_high]
    }
    filename <- paste0(save.dir,'/replicate',i,'.csv')
    write.table(newcount, filename, sep = "\t", row.names = TRUE, col.names = TRUE)

    rm(newcount)
    rm(mean_mat)
    filename
  })
  invisible(unlist(Generate_multiple_counts))
}
