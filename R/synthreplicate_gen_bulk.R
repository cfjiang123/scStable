###############################################################################
# scStable - step 4: generate new bulk
# input : result from step 2 (fit_bulk), result from step 3 (scDesign3_fit)
# output: new bulk
###############################################################################

#' Generate synthetic bulk RNA-seq replicates for scStable
#'
#' \code{synthreplicate_gen_bulk} draws synthetic bulk RNA-seq profiles from the
#' bulk model fitted by \code{\link{fit_bulk}}. A pseudo-bulk profile derived
#' from the single-cell data is used to set per-gene truncation bounds so that
#' sampled (log-)expression values remain biologically valid. Three sampling
#' modes are available via \code{use.cor}.
#'
#' @param bulkRNA_matrix Filtered bulk matrix (genes x samples) from
#'   \code{\link{synthreplicate_prep}}.
#' @param mu Gene-wise mean vector from \code{\link{fit_bulk}}.
#' @param cov Gene-gene covariance matrix from \code{\link{fit_bulk}}.
#' @param optimal_c Per-gene log offsets from \code{\link{fit_bulk}}.
#' @param scRNA_matrix Filtered single-cell matrix (genes x cells) from
#'   \code{\link{synthreplicate_prep}}.
#' @param use.cor Integer; \code{1} truncated multivariate normal (correlated),
#'   \code{2} independent truncated normals, \code{3} plain multivariate normal.
#' @param min.eig Numeric tolerance passed to
#'   \code{corpcor::make.positive.definite} when \code{use.cor = 1}.
#' @param number.replicate Integer; number of synthetic bulk replicates to draw.
#'
#' @return A list with elements \code{d} (per-gene shift between pseudo-bulk and
#'   real bulk) and \code{sampling_bulk} (genes x replicates matrix of synthetic
#'   bulk profiles).
#'
#' @importFrom corpcor make.positive.definite
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS mvrnorm
#' @export
synthreplicate_gen_bulk <-function(bulkRNA_matrix, mu, cov, optimal_c, scRNA_matrix, use.cor = 1,
                              min.eig = 1e-2, number.replicate = 10){
  pseudo_bulk = rowSums(scRNA_matrix)
  names(pseudo_bulk) <- rownames(scRNA_matrix)
  pseudo_bulk = pseudo_bulk * sum(bulkRNA_matrix) / ncol(bulkRNA_matrix) / sum(pseudo_bulk) # scaling
  log_pseudo_bulk = log(pseudo_bulk + optimal_c)
  d = mu - log_pseudo_bulk # shift factor between pseudo bulk and real bulk
  lower_bound = d + log(optimal_c)
  upper_bound = rep(Inf, length(lower_bound))

    if(use.cor == 1){
      cov_p <- make.positive.definite(cov, tol=min.eig)
      sampling_bulk = t(tmvtnorm::rtmvnorm(n = number.replicate, mean = mu, sigma = cov_p,
                                           lower = lower_bound, upper = upper_bound,
                                           algorithm = "gibbs",start = mu))
    }else if(use.cor == 2){
      cov_p = diag(cov)
      samples_list <- lapply(seq_along(mu), function(i) {
        rtruncnorm(
          n     = number.replicate,
          a     = lower_bound[i],
          b     = upper_bound[i],
          mean  = mu[i],
          sd    = sqrt(cov_p[i])
        )
      })
      samples_mat <- do.call(cbind, samples_list)
      sampling_bulk <- t(samples_mat)
    } else if(use.cor == 3){
      sampling_bulk <- replicate(
        number.replicate,
        MASS::mvrnorm(n = 1, mu = mu, Sigma = cov)
      )
    }
    rownames(sampling_bulk) <- names(mu)
    colnames(sampling_bulk) <- paste0("rep", seq_len(number.replicate))
  return(list(d = d, sampling_bulk = sampling_bulk))
}
