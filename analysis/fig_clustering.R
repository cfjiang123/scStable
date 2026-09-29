# Clustering figure: stability of clustering methods and hyperparameters across
# synthetic samples, without using cell-type labels for the selection.
#   A  per-cell entropy stability (Seurat Louvain, lambda = 1) on a reference UMAP
#   B  agreement of six methods with the true cell types vs scale factor lambda
#   C  Louvain hyperparameter grid: stability (synthetic vs reference clustering)
#      selects a setting; reference UMAP coloured by default vs selected clusters
suppressPackageStartupMessages({
  library(Seurat)
  library(mclust)   # Mclust() needs mclustBIC() on the search path
  library(ggplot2)
  library(patchwork)
})
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# Settings scaled down for the 100-cell example. Manuscript (3,555 cells):
# DIMS = 1:30, K_NN = 20, resolution 0.5 / 0.8 / 0.8,
# grid n_pc {10, 20, 50} x k_nn {10, 20, 50} x res {0.4, 0.8, 1.2}, default npc50_knn20_res0.8.
METHODS    <- c("seurat_louvain", "seurat_louvain_refine", "seurat_slm",
                "kmeans", "hclust_ward", "gmm_mclust")
RESOLUTION <- c(seurat_louvain = 0.8, seurat_louvain_refine = 1.2, seurat_slm = 1.2)
K          <- 4     # number of clusters for kmeans / Ward / GMM
DIMS       <- 1:10
K_NN       <- 10
GRID       <- expand.grid(n_pc = c(5, 10, 20), k_nn = c(5, 10, 20), res = c(0.4, 0.8, 1.2))
DEFAULT    <- "npc20_knn20_res0.8"

ref       <- readRDS("data/synthetic/fit.rds")$prep$sc
cell_type <- readRDS("data/cell_types.rds")[colnames(ref)]
lambda_dirs <- list.files("data/synthetic", "^lambda_", full.names = TRUE)
read_samples <- function(dir) lapply(list.files(dir, "^replicate\\d+\\.csv$", full.names = TRUE),
  function(f) as.matrix(read.table(f, sep = "\t", header = TRUE, check.names = FALSE)))

preprocess <- function(counts, npcs) {
  set.seed(1)
  seu <- CreateSeuratObject(as(counts, "CsparseMatrix"))
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  RunPCA(seu, npcs = npcs, verbose = FALSE)
}

cluster_methods <- function(counts) {
  seu <- FindNeighbors(preprocess(counts, max(DIMS)), dims = DIMS, k.param = K_NN, verbose = FALSE)
  X   <- Embeddings(seu, "pca")[, DIMS]
  algo <- c(seurat_louvain = 1, seurat_louvain_refine = 2, seurat_slm = 3)
  set.seed(1)
  sapply(METHODS, function(m) as.integer(switch(m,
    kmeans      = kmeans(X, K, nstart = 10)$cluster,
    hclust_ward = cutree(hclust(dist(X), "ward.D2"), K),
    gmm_mclust  = Mclust(X, G = K, verbose = FALSE)$classification,
    FindClusters(seu, algorithm = algo[[m]], resolution = RESOLUTION[[m]],
                 random.seed = 1, verbose = FALSE)$seurat_clusters)))
}

cluster_grid <- function(counts) {
  seu <- preprocess(counts, max(GRID$n_pc))
  sapply(seq_len(nrow(GRID)), function(i) {
    g  <- GRID[i, ]
    nb <- FindNeighbors(seu, dims = 1:g$n_pc, k.param = g$k_nn, verbose = FALSE)
    as.integer(FindClusters(nb, resolution = g$res, random.seed = 1, verbose = FALSE)$seurat_clusters)
  })
}
grid_keys <- with(GRID, sprintf("npc%d_knn%d_res%s", n_pc, k_nn, res))

agreement <- function(pred, truth) {
  tab <- table(pred, truth); ch2 <- function(n) sum(n * (n - 1) / 2)
  a <- ch2(tab); b <- ch2(rowSums(tab)) - a; c <- ch2(colSums(tab)) - a
  c(ARI = adjustedRandIndex(pred, truth), NMI = aricode::NMI(pred, truth),
    Jaccard = a / (a + b + c), FMI = a / sqrt((a + b) * (a + c)))
}

set.seed(1)
ref_pca  <- irlba::prcomp_irlba(t(log1p(ref)), n = 20, center = TRUE, scale. = TRUE)
ref_umap <- umap::umap(ref_pca$x)$layout
plot_umap <- function(col, title, legend) {
  df <- data.frame(UMAP1 = ref_umap[, 1], UMAP2 = ref_umap[, 2], col = col)
  ggplot(df, aes(UMAP1, UMAP2, color = col)) + geom_point(size = 1.5, alpha = 0.8) +
    labs(title = title, color = legend) + theme_classic() + theme(aspect.ratio = 1)
}

# ---- B: six methods x lambda, agreement with the true cell types ----
labels <- lapply(lambda_dirs, function(d) lapply(read_samples(d), cluster_methods))
names(labels) <- sub("lambda_", "", basename(lambda_dirs))
metrics <- do.call(rbind, lapply(names(labels), function(l) do.call(rbind,
  lapply(seq_along(labels[[l]]), function(r) do.call(rbind, lapply(METHODS, function(m)
    data.frame(lambda = as.numeric(l), sample = r, method = m,
               t(agreement(labels[[l]][[r]][, m], cell_type)))))))))
write.csv(metrics, "results/clustering_B_method_metrics.csv", row.names = FALSE)
mean_metrics <- aggregate(cbind(ARI, NMI, Jaccard, FMI) ~ method + lambda, metrics, mean)
long <- do.call(rbind, lapply(c("ARI", "NMI", "Jaccard", "FMI"), function(k)
  data.frame(mean_metrics[, c("method", "lambda")], metric = k, value = mean_metrics[[k]])))
p_b <- ggplot(long, aes(lambda, value, color = method)) + geom_line(linewidth = 1) +
  facet_wrap(~metric, nrow = 2) + labs(x = "Scale factor", y = NULL, color = NULL) +
  theme_classic() + theme(legend.position = "bottom")
ggsave("figures/clustering_B_method_metrics.pdf", p_b, width = 9, height = 7)

# ---- A: entropy stability of each cell (Seurat Louvain, lambda = 1) ----
lab1 <- sapply(labels[["1"]], function(x) x[, "seurat_louvain"])
consensus <- Reduce(`+`, lapply(seq_len(ncol(lab1)), function(r) outer(lab1[, r], lab1[, r], "=="))) / ncol(lab1)
p <- pmin(pmax(consensus, 1e-10), 1 - 1e-10)
H <- -p * log2(p) - (1 - p) * log2(1 - p)
diag(H) <- NA
stability <- 1 - rowMeans(H, na.rm = TRUE)
write.csv(data.frame(cell = colnames(ref), cell_type, entropy_stability = stability),
          "results/clustering_A_entropy_stability.csv", row.names = FALSE)
p_a <- plot_umap(cell_type, "Cell types", "Cell type") |
  (plot_umap(stability, "Entropy stability", "Stability") + scale_color_viridis_c(limits = c(0, 1)))
ggsave("figures/clustering_A_entropy_stability.pdf", p_a, width = 10, height = 4.5)

# ---- C: Louvain grid; stability = agreement of synthetic with reference clustering ----
grid_ref <- cluster_grid(ref)
grid_syn <- lapply(read_samples("data/synthetic/lambda_1"), cluster_grid)
grid_res <- data.frame(param = grid_keys, GRID, t(sapply(seq_len(nrow(GRID)), function(i) c(
  stability_ARI = mean(sapply(grid_syn, function(s) adjustedRandIndex(s[, i], grid_ref[, i]))),
  truth_ARI     = adjustedRandIndex(grid_ref[, i], cell_type),
  n_clusters    = length(unique(grid_ref[, i]))))))
# a single cluster is trivially stable, so it is not eligible for selection
grid_res <- grid_res[order(grid_res$n_clusters < 2, -grid_res$stability_ARI), ]
write.csv(grid_res, "results/clustering_C_louvain_grid.csv", row.names = FALSE)
selected <- grid_res$param[1]
cat(sprintf("Default %s: truth ARI %.3f | selected %s: truth ARI %.3f\n",
            DEFAULT, grid_res$truth_ARI[grid_res$param == DEFAULT],
            selected, grid_res$truth_ARI[grid_res$param == selected]))

p_c <- plot_umap(cell_type, "Reference cell type", "Cell type") |
  plot_umap(factor(grid_ref[, grid_keys == DEFAULT]), paste("Default:", DEFAULT), "Cluster") |
  plot_umap(factor(grid_ref[, grid_keys == selected]), paste("Selected:", selected), "Cluster")
ggsave("figures/clustering_C_default_vs_selected.pdf", p_c, width = 15, height = 4.5)
p_c2 <- ggplot(grid_res, aes(stability_ARI, truth_ARI)) + geom_point(size = 2) +
  labs(x = "Stability (synthetic vs reference ARI)", y = "ARI vs true cell types") + theme_bw()
ggsave("figures/clustering_C_grid_stability_vs_truth.pdf", p_c2, width = 5, height = 4.5)
