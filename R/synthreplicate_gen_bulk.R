###############################################################################
# scStable - Step 3 (part 1): draw synthetic bulk samples
###############################################################################

#' Step 3 (part 1): draw synthetic bulk RNA-seq samples
#'
#' \code{synthreplicate_gen_bulk} draws \code{number.replicate} synthetic bulk
#' profiles on the log scale from the Gaussian bulk model
#' \eqn{N(\hat\mu, \hat\Sigma)} and computes the gene-specific offset
#' \eqn{\xi_g} (returned as \code{d}) that aligns the scRNA-seq pseudo-bulk with
#' the bulk RNA-seq reference.
#'
#' @param bulkRNA_matrix Bulk matrix (genes x samples) from
#'   \code{\link{synthreplicate_prep}} (reference mode). Only its mean library
#'   size and number of samples are used. May be \code{NULL} if \code{lib_size}
#'   is given (bulk-reference-free mode).
#' @param mu,cov,optimal_c Gene-wise mean, gene-gene covariance and
#'   pseudo-counts from \code{\link{fit_bulk}} (or a pre-estimated parameter
#'   set), in the same gene order as \code{scRNA_matrix}.
#' @param scRNA_matrix Single-cell count matrix (genes x cells) from
#'   \code{\link{synthreplicate_prep}}.
#' @param use.cor Covariance used for sampling: \code{3} full covariance
#'   (multivariate normal), \code{2} diagonal covariance (independent truncated
#'   normals), \code{1} truncated multivariate normal (Gibbs sampler). The
#'   default \code{NULL} follows the manuscript: \code{3} if the bulk reference
#'   has at least 30 samples, otherwise \code{2}.
#' @param min.eig Tolerance passed to \code{corpcor::make.positive.definite}
#'   when \code{use.cor = 1}.
#' @param number.replicate Number of synthetic samples \eqn{R} to draw.
#' @param lib_size Mean bulk library size. Defaults to
#'   \code{sum(bulkRNA_matrix) / ncol(bulkRNA_matrix)}; supply it (e.g.
#'   \code{sum(fit_bulk()$bulk_mean)}) when \code{bulkRNA_matrix} is \code{NULL}.
#' @param n_samples Number of bulk samples behind \code{cov}, used only to
#'   choose \code{use.cor} when it is \code{NULL}. Defaults to
#'   \code{ncol(bulkRNA_matrix)}.
#'
#' @return A list with \code{d} (per-gene offset \eqn{\xi_g}) and
#'   \code{sampling_bulk} (genes x \code{number.replicate} matrix of synthetic
#'   bulk profiles on the log scale).
#'
#' @importFrom corpcor make.positive.definite
#' @importFrom tmvtnorm rtmvnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom MASS mvrnorm
#' @export
synthreplicate_gen_bulk <- function(bulkRNA_matrix = NULL, mu, cov, optimal_c, scRNA_matrix,
                                    use.cor = NULL, min.eig = 1e-2, number.replicate = 100,
                                    lib_size = NULL, n_samples = NULL) {
  use_bulk <- is.null(lib_size)
  if (use_bulk && is.null(bulkRNA_matrix)) stop("Provide either bulkRNA_matrix or lib_size")
  if (is.null(n_samples) && !is.null(bulkRNA_matrix)) n_samples <- ncol(bulkRNA_matrix)
  if (is.null(use.cor)) {
    if (is.null(n_samples)) stop("Set use.cor, or provide bulkRNA_matrix / n_samples")
    use.cor <- if (n_samples >= 30) 3 else 2
  }
  if (!identical(names(mu), rownames(scRNA_matrix))) {
    stop("names(mu) must match rownames(scRNA_matrix) in the same order")
  }

  pseudo_bulk = rowSums(scRNA_matrix)
  names(pseudo_bulk) <- rownames(scRNA_matrix)
  pseudo_bulk = if (use_bulk) {
    pseudo_bulk * sum(bulkRNA_matrix) / ncol(bulkRNA_matrix) / sum(pseudo_bulk) # scaling
  } else {
    pseudo_bulk * lib_size / sum(pseudo_bulk)
  }
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
      cov_p = if (is.matrix(cov)) diag(cov) else cov
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
    } else {
      stop("use.cor must be 1, 2 or 3")
    }
    rownames(sampling_bulk) <- names(mu)
    colnames(sampling_bulk) <- paste0("rep", seq_len(number.replicate))
  return(list(d = d, sampling_bulk = sampling_bulk))
}
