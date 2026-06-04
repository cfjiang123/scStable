###############################################################################
# scStable - step 2: fit bulk distribution
# input : filtered bulk RNA-seq matrix (with names) (result from step 1)
# output: bulk RNA-seq matrix distribution, mu, cov, gene-specific constant
###############################################################################

#' Fit the bulk RNA-seq distribution for scStable
#'
#' \code{fit_bulk} estimates, for each gene, a log-offset constant \code{c} that
#' best normalises the bulk RNA-seq counts (\code{log(count + c)}), then computes
#' the gene-wise mean vector and the gene-gene covariance matrix used by
#' \pkg{scStable} to sample synthetic bulk replicates. The search for the
#' optimal \code{c} is parallelised across genes with \pkg{foreach} /
#' \pkg{doParallel}.
#'
#' @param bulkRNA_matrix Numeric matrix of bulk RNA-seq counts (genes x samples),
#'   typically the \code{bulk} element returned by \code{\link{synthreplicate_prep}}.
#' @param n_cores Integer; number of worker cores for the parallel grid search.
#' @param p_val_threshold Numeric; normality-test p-value at which the search for
#'   \code{c} stops early.
#' @param step Numeric; increment of the \code{c} grid.
#' @param c_value_max Numeric; maximum value of \code{c} to evaluate.
#' @param normal_test Character; \code{"auto"} chooses Shapiro-Wilk for small
#'   samples and Anderson-Darling otherwise; may also be \code{"shapiro"} or
#'   \code{"ad"}.
#' @param sw_limit Integer; sample-size ceiling under which Shapiro-Wilk is used
#'   when \code{normal_test = "auto"}.
#'
#' @return A list with elements \code{bulk_data_counts} (log-transformed counts),
#'   \code{mu} (gene means), \code{cov} (gene-gene covariance) and
#'   \code{optimal_c} (per-gene log offsets).
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @importFrom nortest ad.test
#' @importFrom stats shapiro.test cov setNames
#' @export
fit_bulk <- function(
    bulkRNA_matrix,
    n_cores         = 20,
    p_val_threshold = 0.5,
    step            = 1,
    c_value_max     = 200,
    normal_test     = c("auto"), # , "shapiro", "ad"
    sw_limit        = 200
) {
  # -- 0. Define internal helper: choose c to maximize normality p-value --
  find_optimal_c <- function(
    gene_counts,
    p_val_threshold,
    step,
    c_value_max,
    normal_test,
    sw_limit
  ) {
    optimal_c       <- 1
    optimal_p_value <- 0

    for (c_val in seq(1, c_value_max, by = step)) {
      log_counts <- log(gene_counts + c_val)
      n_samples  <- length(log_counts)

      # pick normality test
      if (normal_test == "shapiro" ||
          (normal_test == "auto" && n_samples <= sw_limit)) {
        p_val <- tryCatch(shapiro.test(log_counts)$p.value,
                          error = function(e) NA_real_)
      } else {
        p_val <- tryCatch(nortest::ad.test(log_counts)$p.value,
                          error = function(e) NA_real_)
      }

      # update if improved
      if (!is.na(p_val) && p_val > optimal_p_value) {
        optimal_c       <- c_val
        optimal_p_value <- p_val
        if (optimal_p_value >= p_val_threshold) break
      }
    }

    list(optimal_c       = optimal_c,
         optimal_p_value = optimal_p_value)
  }

  # -- 1. Spin up cluster and register for foreach --
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  # -- 2. Parallel loop over genes --
  results <- foreach::foreach(
    i         = seq_len(nrow(bulkRNA_matrix)),
    .combine  = rbind,
    .packages = "nortest"
  ) %dopar% {
    find_optimal_c(
      gene_counts    = bulkRNA_matrix[i, ],
      p_val_threshold = p_val_threshold,
      step            = step,
      c_value_max     = c_value_max,
      normal_test     = normal_test,
      sw_limit        = sw_limit
    )
  }

  # -- 3. Tear down cluster --
  parallel::stopCluster(cl)

  # -- 4. Assign row & column names --
  rownames(results) <- rownames(bulkRNA_matrix)
  colnames(results) <- c("optimal_c", "optimal_p_value")

  # -- 5. Extract c-values and log-transform counts --
  optimal_c        <- setNames(as.numeric(results[, "optimal_c"]),
                               rownames(results))
  bulk_data_counts <- log(bulkRNA_matrix + optimal_c)

  # -- 6. Compute summary statistics --
  mu    <- rowMeans(bulk_data_counts)
  cov_m <- stats::cov(t(bulk_data_counts))

  # -- 7. Return as list --
  list(
    bulk_data_counts = bulk_data_counts,
    mu               = mu,
    cov              = cov_m,
    optimal_c        = optimal_c
  )
}
