###############################################################################
# scStable - sample stability metrics for synthetic scRNA-seq samples
#
# Load with:
#   source(system.file("scripts", "scStable_stability.R", package = "scStable"))
#
# Notation follows the manuscript (Methods, "Sample stability metrics"):
#   Y_0       single-sample scRNA-seq reference (genes x cells counts)
#   Y_1..Y_R  synthetic scRNA-seq samples written by synthreplicate_gen_sc()
#
# Functions
#   list_synthetic_samples()      paths of replicate<r>.csv files, in order of r
#   read_synthetic_sample()       read one synthetic sample as a count matrix
#   quantile_interval()           (1 - alpha) quantile interval and length L_alpha
#   de_padj()                     Seurat Bonferroni-adjusted p-values for one sample
#   stable_de()                   stable DE genes (p_adj < alpha in Y_0 and L_alpha < L_cut)
#   entropy_stability()           per-cell entropy-based clustering stability S_i
#   cluster_agreement()           ARI, NMI, Jaccard and FMI between two labelings
#   seurat_louvain_grid()         Seurat Louvain labels for a hyperparameter grid
#   stability_select()            Stab(A_l): agreement of synthetic vs reference output
#   stability_path()              stability_select() along a path of scale factors lambda
#
# Arguments named `samples` accept either a list of count matrices or a
# character vector of replicate file paths; files are read one at a time.
###############################################################################


# ---------------------------------------------------------------------------
# Reading synthetic samples
# ---------------------------------------------------------------------------

#' Paths of the synthetic samples in a directory, ordered by replicate index r.
#' Input : dir     directory passed as save.dir to synthreplicate_gen_sc()
#'         pattern file name pattern
#' Output: character vector of file paths
list_synthetic_samples <- function(dir, pattern = "^replicate\\d+\\.csv$") {
  f <- list.files(dir, pattern, full.names = TRUE)
  if (!length(f)) stop("No synthetic samples matching '", pattern, "' in ", dir)
  f[order(as.integer(gsub("\\D", "", basename(f))))]
}

#' Read one synthetic sample.
#' Input : file  path of a replicate<r>.csv file (tab-separated, genes x cells)
#' Output: numeric matrix, genes x cells, with gene and cell names
read_synthetic_sample <- function(file) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    # write.table() header has one field fewer than the data rows (no row-name column)
    cells <- scan(file, what = "", sep = "\t", nlines = 1, quiet = TRUE)
    d <- data.table::fread(file, sep = "\t", skip = 1, header = FALSE, data.table = FALSE)
    m <- as.matrix(d[, -1, drop = FALSE])
    dimnames(m) <- list(d[[1]], cells)
    m
  } else {
    as.matrix(read.table(file, sep = "\t", header = TRUE, check.names = FALSE))
  }
}

.get_sample <- function(samples, r) {
  if (is.character(samples)) read_synthetic_sample(samples[[r]]) else samples[[r]]
}


# ---------------------------------------------------------------------------
# Quantile-interval stability (Eq. rep_ci)
# ---------------------------------------------------------------------------

#' (1 - alpha) sample-based quantile interval of a scalar summary T.
#' Input : T_syn  features x R numeric matrix; column r holds T(A(Y_r)) for every
#'                feature (gene, edge, ...). A list of R named vectors is also accepted.
#'         T_ref  optional named vector T(A(Y_0)) on the reference
#'         alpha  interval level; the interval is [T_{alpha/2}, T_{1-alpha/2}]
#'         na_as  value substituted for NA in T_syn before taking quantiles
#'                (NULL drops NA values)
#' Output: data.frame with one row per feature:
#'         feature, ref (if T_ref given), lower, upper, L (= upper - lower), n (non-NA samples)
quantile_interval <- function(T_syn, T_ref = NULL, alpha = 0.05, na_as = NULL) {
  if (is.list(T_syn) && !is.data.frame(T_syn)) {
    feats <- unique(unlist(lapply(T_syn, names)))
    T_syn <- sapply(T_syn, function(v) v[feats])
    rownames(T_syn) <- feats
  }
  T_syn <- as.matrix(T_syn)
  if (!is.null(na_as)) T_syn[is.na(T_syn)] <- na_as
  q <- t(apply(T_syn, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) >= 2) quantile(x, c(alpha / 2, 1 - alpha / 2), names = FALSE) else c(NA, NA)
  }))
  out <- data.frame(feature = rownames(T_syn), lower = q[, 1], upper = q[, 2],
                    L = q[, 2] - q[, 1], n = rowSums(!is.na(T_syn)), row.names = rownames(T_syn))
  if (!is.null(T_ref)) out <- cbind(out[, 1, drop = FALSE], ref = unname(T_ref[out$feature]), out[, -1])
  out
}


# ---------------------------------------------------------------------------
# Stable DE genes
# ---------------------------------------------------------------------------

#' Seurat Bonferroni-adjusted p-values of one DE test on one sample.
#' Input : counts    genes x cells count matrix
#'         group     cell labels, named by cell or in the column order of counts
#'         ident.1, ident.2  the two groups compared
#'         test.use  any Seurat FindMarkers test ("wilcox", "negbinom", "bimod", "t", "LR", "poisson", ...)
#'         na_as     p_adj given to genes that FindMarkers does not test
#'                   (filtered by min.pct / logfc.threshold)
#'         ...       further arguments passed to Seurat::FindMarkers
#' Output: named numeric vector of p_adj for every gene in rownames(counts)
de_padj <- function(counts, group, ident.1, ident.2, test.use = "wilcox", na_as = 1,
                    logfc.threshold = 0.1, min.pct = 0.1, ...) {
  if (!is.null(names(group))) group <- group[colnames(counts)]
  genes <- rownames(counts)
  seu <- Seurat::CreateSeuratObject(Matrix::Matrix(as.matrix(counts), sparse = TRUE),
                                    min.cells = 0, min.features = 0)
  lookup <- setNames(genes, rownames(seu))        # Seurat may rename genes (e.g. "_" -> "-")
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  Seurat::Idents(seu) <- factor(unname(group))
  res <- Seurat::FindMarkers(seu, ident.1 = ident.1, ident.2 = ident.2, test.use = test.use,
                             logfc.threshold = logfc.threshold, min.pct = min.pct,
                             verbose = FALSE, ...)
  p <- setNames(rep(na_as, length(genes)), genes)
  p[lookup[rownames(res)]] <- res$p_val_adj
  p
}

#' Stable DE genes across synthetic samples.
#' Between cell types: one reference, cells labelled by `group`.
#' Between conditions: `ref`/`samples` hold condition ident.1 and `ref2`/`samples2`
#' hold condition ident.2 (generated separately); sample r of each set is paired.
#' Input : ref       genes x cells reference counts Y_0
#'         samples   synthetic samples Y_1..Y_R (list of matrices or file paths)
#'         group     cell labels of ref (between cell types); ignored if ref2 is given
#'         ident.1, ident.2  groups (or conditions) compared
#'         ref2, samples2    second-condition reference and synthetic samples (optional)
#'         test.use  Seurat DE test
#'         alpha     significance level for reference p_adj and interval level (default 0.05)
#'         L_cut     stability cut-off on the interval length L_alpha (default 0.05)
#'         na_as     p_adj of untested genes (see de_padj)
#'         ...       passed to Seurat::FindMarkers
#' Output: data.frame with one row per gene:
#'         gene, ref_padj, lower, upper, L, n, de (ref_padj < alpha),
#'         stable (de & L < L_cut), sorted by ref_padj.
#'         The per-sample p_adj matrix (genes x R) is attached as attr(, "padj").
stable_de <- function(ref, samples, group = NULL, ident.1, ident.2, ref2 = NULL, samples2 = NULL,
                      test.use = "wilcox", alpha = 0.05, L_cut = 0.05, na_as = 1, ...) {
  pair <- !is.null(ref2)
  if (pair && length(samples2) != length(samples)) stop("samples and samples2 must have the same length")
  run <- function(a, b) {
    if (!pair) return(de_padj(a, group, ident.1, ident.2, test.use, na_as, ...))
    genes <- intersect(rownames(a), rownames(b))
    m <- cbind(a[genes, , drop = FALSE], b[genes, , drop = FALSE])
    colnames(m) <- make.unique(c(paste0("c1_", colnames(a)), paste0("c2_", colnames(b))))
    g <- rep(c(ident.1, ident.2), c(ncol(a), ncol(b)))
    de_padj(m, g, ident.1, ident.2, test.use, na_as, ...)
  }
  p_ref <- run(ref, ref2)
  p_syn <- sapply(seq_along(samples), function(r)
    run(.get_sample(samples, r), if (pair) .get_sample(samples2, r))[names(p_ref)])
  rownames(p_syn) <- names(p_ref)
  qi <- quantile_interval(p_syn, T_ref = p_ref, alpha = alpha)
  out <- data.frame(gene = qi$feature, ref_padj = qi$ref, lower = qi$lower, upper = qi$upper,
                    L = qi$L, n = qi$n, row.names = qi$feature)
  out$de     <- out$ref_padj < alpha
  out$stable <- out$de & !is.na(out$L) & out$L < L_cut
  out <- out[order(out$ref_padj), ]
  attr(out, "padj") <- p_syn
  out
}


# ---------------------------------------------------------------------------
# Entropy-based per-cell clustering stability (Eqs. cocluster_matrix - entropy_stability)
# ---------------------------------------------------------------------------

#' Per-cell entropy stability S_i from clusterings of the R synthetic samples.
#' Input : labels      cells x R matrix / data.frame of cluster labels, or a list of R
#'                     label vectors (same cell order in every sample)
#'         chunk_size  rows of the consensus matrix computed at a time (memory control)
#' Output: numeric vector S_i in [0, 1] (names = cell names if available);
#'         1 = decisive, consistent assignments across synthetic samples
entropy_stability <- function(labels, chunk_size = 1000) {
  if (is.list(labels) && !is.data.frame(labels)) labels <- do.call(cbind, lapply(labels, as.character))
  labels <- apply(as.matrix(labels), 2, function(x) match(x, unique(x)))
  n <- nrow(labels); R <- ncol(labels)
  S <- numeric(n)
  for (start in seq(1, n, by = chunk_size)) {
    idx <- start:min(n, start + chunk_size - 1)
    C <- matrix(0, length(idx), n)
    for (r in seq_len(R)) C <- C + outer(labels[idx, r], labels[, r], "==")
    C <- C / R
    H <- -C * log2(C) - (1 - C) * log2(1 - C)
    H[is.nan(H)] <- 0                              # H = 0 when C_ij is 0 or 1 (includes j = i)
    S[idx] <- 1 - rowSums(H) / (n - 1)
  }
  names(S) <- rownames(labels)
  S
}


# ---------------------------------------------------------------------------
# Stability-driven method and hyperparameter selection (Eqs. stability_score, variation_path)
# ---------------------------------------------------------------------------

#' Agreement between two labelings of the same cells.
#' Input : pred, ref  label vectors of equal length
#' Output: named vector c(ARI, NMI, Jaccard, FMI); NMI is normalised by max entropy
cluster_agreement <- function(pred, ref) {
  tab <- table(pred, ref)
  n   <- sum(tab)
  ch2 <- function(x) sum(x * (x - 1) / 2)
  a   <- ch2(tab); ra <- ch2(rowSums(tab)); ca <- ch2(colSums(tab)); N <- n * (n - 1) / 2
  expct <- ra * ca / N; mx <- (ra + ca) / 2
  ari <- if (mx == expct) as.numeric(a == mx) else (a - expct) / (mx - expct)
  ent <- function(p) { p <- p[p > 0] / n; -sum(p * log(p)) }
  pij <- tab / n; pi <- rowSums(pij); pj <- colSums(pij); nz <- pij > 0
  mi  <- sum(pij[nz] * log(pij[nz] / outer(pi, pj)[nz]))
  hmax <- max(ent(rowSums(tab)), ent(colSums(tab)))
  c(ARI = ari, NMI = if (hmax == 0) 1 else mi / hmax,
    Jaccard = a / (ra + ca - a), FMI = if (ra * ca == 0) 0 else a / sqrt(ra * ca))
}

#' Seurat Louvain clustering over a grid of hyperparameters (PCA computed once).
#' Usable directly as `run_fun` in stability_select().
#' Input : counts      genes x cells count matrix
#'         n_pc, k_nn, resolution  candidate values; all combinations are run
#'         seed        random seed
#' Output: integer matrix, cells x specifications, column names "npc<n>_knn<k>_res<r>"
seurat_louvain_grid <- function(counts, n_pc = c(10, 20, 50), k_nn = c(10, 20, 50),
                                resolution = c(0.4, 0.8, 1.2), seed = 1) {
  set.seed(seed)
  seu <- Seurat::CreateSeuratObject(Matrix::Matrix(as.matrix(counts), sparse = TRUE))
  seu <- Seurat::NormalizeData(seu, verbose = FALSE)
  seu <- Seurat::FindVariableFeatures(seu, verbose = FALSE)
  seu <- Seurat::ScaleData(seu, verbose = FALSE)
  seu <- Seurat::RunPCA(seu, npcs = max(n_pc), verbose = FALSE)
  grid <- expand.grid(n_pc = n_pc, k_nn = k_nn, res = resolution)
  out <- sapply(seq_len(nrow(grid)), function(i) {
    g  <- grid[i, ]
    nb <- Seurat::FindNeighbors(seu, dims = 1:g$n_pc, k.param = g$k_nn, verbose = FALSE)
    as.integer(Seurat::FindClusters(nb, resolution = g$res, random.seed = seed,
                                    verbose = FALSE)$seurat_clusters)
  })
  dimnames(out) <- list(colnames(counts), with(grid, sprintf("npc%d_knn%d_res%s", n_pc, k_nn, res)))
  out
}

#' Stability of candidate analysis specifications A_1..A_L:
#' Stab(A_l) = Agg( d(A_l(Y_r), A_l(Y_0)) : r = 1..R ).
#' Input : ref      genes x cells reference counts Y_0
#'         samples  synthetic samples Y_1..Y_R (list of matrices or file paths)
#'         run_fun  function(counts) returning a cells x L matrix of labels, one
#'                  column per specification (column names = specification names),
#'                  e.g. seurat_louvain_grid or a function running several methods
#'         agg      aggregation over samples (default mean)
#'         rank_by  metric used to order the specifications
#' Output: list with
#'         summary     data.frame, one row per specification: spec, n_clusters_ref,
#'                     aggregated ARI, NMI, Jaccard, FMI; ordered by rank_by (decreasing),
#'                     specifications with fewer than 2 reference clusters last
#'         per_sample  data.frame: spec, sample, ARI, NMI, Jaccard, FMI
#'         ref_labels  cells x L reference labels
stability_select <- function(ref, samples, run_fun, agg = mean, rank_by = "ARI") {
  ref_lab <- as.matrix(run_fun(ref))
  specs <- colnames(ref_lab)
  per <- do.call(rbind, lapply(seq_along(samples), function(r) {
    lab <- as.matrix(run_fun(.get_sample(samples, r)))
    do.call(rbind, lapply(specs, function(s)
      data.frame(spec = s, sample = r, t(cluster_agreement(lab[, s], ref_lab[, s])))))
  }))
  metrics <- c("ARI", "NMI", "Jaccard", "FMI")
  summ <- do.call(rbind, lapply(specs, function(s) {
    x <- per[per$spec == s, metrics]
    data.frame(spec = s, n_clusters_ref = length(unique(ref_lab[, s])),
               t(sapply(x, agg)), row.names = NULL)
  }))
  summ <- summ[order(summ$n_clusters_ref < 2, -summ[[rank_by]]), ]
  rownames(summ) <- NULL
  list(summary = summ, per_sample = per, ref_labels = ref_lab)
}

#' Stability of each specification along a path of scale factors lambda.
#' Input : ref               genes x cells reference counts Y_0
#'         samples_by_lambda named list, one element per lambda (names = lambda values),
#'                           each a list of matrices or file paths
#'         run_fun, agg      as in stability_select()
#' Output: data.frame: lambda, spec, n_clusters_ref, ARI, NMI, Jaccard, FMI
stability_path <- function(ref, samples_by_lambda, run_fun, agg = mean) {
  do.call(rbind, lapply(names(samples_by_lambda), function(l) {
    s <- stability_select(ref, samples_by_lambda[[l]], run_fun, agg)$summary
    cbind(lambda = as.numeric(l), s)
  }))
}
