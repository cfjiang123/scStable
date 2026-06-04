###############################################################################
# scStable - step 5: generate new sc
# input : result from step 2, result from step 3, result from step 4
# output: new scRNA-seq
###############################################################################

#' Generate synthetic single-cell RNA-seq replicates for scStable
#'
#' \code{synthreplicate_gen_sc} is the final step of the \pkg{scStable}
#' workflow. For each synthetic bulk replicate it derives a per-gene mapping
#' (fold-change) factor, re-weights the scDesign3 mean matrix accordingly, and
#' simulates a new single-cell count matrix with \code{scDesign3::simu_new}. Each
#' synthetic replicate is written to disk as a tab-separated file.
#'
#' @param bulkRNA_matrix Filtered bulk matrix from \code{\link{synthreplicate_prep}}.
#' @param bulk_synth Synthetic bulk matrix (genes x replicates) from
#'   \code{\link{synthreplicate_gen_bulk}}.
#' @param mu Gene-wise mean vector from \code{\link{fit_bulk}}.
#' @param d Per-gene shift vector from \code{\link{synthreplicate_gen_bulk}}.
#' @param optimal_c Per-gene log offsets from \code{\link{fit_bulk}}.
#' @param scRNA_matrix Filtered single-cell matrix from
#'   \code{\link{synthreplicate_prep}}.
#' @param scRNA_list The list returned by \code{\link{scDesign3_fit}}.
#' @param match.option Integer; \code{1} precise per-gene matching, \code{2}
#'   pseudo-bulk ratio matching.
#' @param scaling_factor Numeric exponent applied to the per-gene mapping.
#' @param use.pc Logical / integer flag controlling the high fold-change
#'   correction step.
#' @param n_cores Integer; cores passed to \code{scDesign3::simu_new}.
#' @param sc_quantile Numeric quantile used for capping (reserved).
#' @param save.dir Directory in which synthetic replicate files are written.
#'
#' @return Called for its side effects: writes one
#'   \code{replicate<i>.csv} file per synthetic replicate to \code{save.dir}.
#'   Returns \code{NULL} invisibly.
#'
#' @importFrom scDesign3 simu_new
#' @importFrom BiocParallel MulticoreParam
#' @importFrom stats quantile IQR median
#' @importFrom utils write.table
#' @export
synthreplicate_gen_sc <-function(bulkRNA_matrix, bulk_synth, mu, d, optimal_c, scRNA_matrix, scRNA_list,
                                 match.option = 1, scaling_factor = 1, use.pc = 1, n_cores = 1, sc_quantile = 0.995, save.dir){
  if (!dir.exists(save.dir)) dir.create(save.dir, recursive = TRUE)
  pseudo_bulk = rowSums(scRNA_matrix)
  e_v = quantile(pseudo_bulk, 0.75, na.rm = TRUE) + 1.5 * IQR(pseudo_bulk, na.rm = TRUE)
  names(pseudo_bulk) <- rownames(scRNA_matrix)
  alpha = sum(bulkRNA_matrix) / ncol(bulkRNA_matrix) / sum(pseudo_bulk) # scaling
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

  })
}
