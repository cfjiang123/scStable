# Validation figure (Fig. 2): synthetic samples vs the scRNA-seq reference and
# the bulk RNA-seq reference, using the lambda = 1 samples.
#   A  UMAP of reference and synthetic samples (projected), mLISI
#   B  gene-gene correlation: bulk reference vs synthetic pseudo-bulk
#   C  gene- and cell-level summary statistics
#   D  gene-wise marginal distributions: bulk vs synthetic pseudo-bulk (KS test)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

fit       <- readRDS("data/synthetic/fit.rds")
prep      <- fit$prep
ref       <- prep$sc
cell_type <- readRDS("data/cell_types.rds")[colnames(ref)]

sample_files <- list.files("data/synthetic/lambda_1", "^replicate\\d+\\.csv$", full.names = TRUE)
sample_files <- sample_files[order(as.integer(gsub("\\D", "", basename(sample_files))))]
samples <- lapply(sample_files, function(f)
  as.matrix(read.table(f, sep = "\t", header = TRUE, check.names = FALSE)))
names(samples) <- paste("Synthetic", seq_along(samples))
stopifnot(all(sapply(samples, function(m) identical(dimnames(m), dimnames(ref)))))

# ---- A: UMAP fitted on the reference, synthetic samples projected; mLISI ----
set.seed(1)
pca_fit  <- irlba::prcomp_irlba(t(log1p(ref)), n = 20, center = TRUE, scale. = TRUE)
umap_fit <- umap::umap(pca_fit$x)
project  <- function(m) predict(umap_fit, predict(pca_fit, newdata = t(log1p(m))))
umaps    <- c(list(Reference = umap_fit$layout), lapply(samples[1:2], project))

mlisi <- sapply(umaps[-1], function(u) {
  X <- rbind(umaps$Reference, u)
  meta <- data.frame(sample = rep(c("ref", "syn"), each = nrow(u)))
  median(lisi::compute_lisi(X, meta, "sample", perplexity = 30)$sample)
})
write.csv(data.frame(sample = names(mlisi), mLISI = mlisi), "results/validation_A_mlisi.csv",
          row.names = FALSE)

umap_df <- do.call(rbind, lapply(names(umaps), function(n)
  data.frame(sample = n, UMAP1 = umaps[[n]][, 1], UMAP2 = umaps[[n]][, 2], cell_type = cell_type)))
umap_df$sample <- factor(umap_df$sample, levels = names(umaps))
p_a <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~sample, nrow = 1) +
  theme_classic() + theme(aspect.ratio = 1, legend.position = "bottom")
ggsave("figures/validation_A_umap.pdf", p_a, width = 10, height = 4.5)

# ---- pseudo-bulk of each synthetic sample on the bulk model scale log(x + c_g) ----
optimal_c  <- fit$bulk_fit$optimal_c
alpha      <- sum(prep$bulk) / ncol(prep$bulk) / sum(ref)
log_pseudo <- log(sapply(samples, rowSums) * alpha + optimal_c) + fit$gen_bulk$d
log_bulk   <- fit$bulk_fit$bulk_data_counts

# ---- B: gene-gene correlation, bulk vs pseudo-bulk ----
cor_pseudo <- cor(t(log_pseudo)); cor_pseudo[is.na(cor_pseudo)] <- 0
cor_bulk   <- cor(t(log_bulk))
gene_order <- rownames(cor_pseudo)[hclust(as.dist((1 - cor_pseudo) / 2))$order]
upper      <- upper.tri(cor_bulk)
write.csv(data.frame(pearson  = cor(cor_bulk[upper], cor_pseudo[upper]),
                     spearman = cor(cor_bulk[upper], cor_pseudo[upper], method = "spearman")),
          "results/validation_B_cor_similarity.csv", row.names = FALSE)

plot_cor <- function(m, title) {
  df <- as.data.frame(as.table(m[gene_order, gene_order]))
  ggplot(df, aes(Var2, Var1, fill = Freq)) + geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-1, 1),
                         name = "Pearson\ncorrelation") +
    coord_fixed() + labs(x = NULL, y = NULL, title = title) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}
p_b <- (plot_cor(cor_bulk, "Reference bulk") | plot_cor(cor_pseudo, "Synthetic pseudo-bulk")) +
  plot_layout(guides = "collect")
ggsave("figures/validation_B_correlation.pdf", p_b, width = 10, height = 5)

# ---- C: gene- and cell-level summary statistics ----
pca_c <- irlba::prcomp_irlba(t(log1p(ref)), n = 50, center = TRUE, scale. = FALSE)
summary_stats <- function(m, pcs) {
  lc <- log1p(m)
  cc <- cor(lc)
  list(`Mean log expression` = rowMeans(lc),
       `Var log expression`  = apply(lc, 1, var),
       `Log library size`    = log(colSums(m)),
       `Cell distance`       = c(dist(pcs)),
       `Cell zero fraction`  = colMeans(m == 0),
       `Cell correlation`    = cc[upper.tri(cc)])
}
all_samples <- c(list(Reference = ref), samples)
stats_df <- do.call(rbind, lapply(names(all_samples), function(n) {
  m   <- all_samples[[n]]
  pcs <- if (n == "Reference") pca_c$x else predict(pca_c, newdata = t(log1p(m)))
  s   <- summary_stats(m, pcs)
  do.call(rbind, lapply(names(s), function(k) data.frame(sample = n, metric = k, value = s[[k]])))
}))
stats_df$sample <- factor(stats_df$sample, levels = names(all_samples))
stats_df$metric <- factor(stats_df$metric, levels = unique(stats_df$metric))
stats_df$type   <- ifelse(stats_df$sample == "Reference", "Reference", "Synthetic")
p_c <- ggplot(stats_df, aes(sample, value, color = type)) +
  geom_violin(scale = "width") +
  facet_wrap(~metric, scales = "free", nrow = 2) +
  scale_color_manual(values = c(Reference = "#00BFC4", Synthetic = "#F8766D")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())
ggsave("figures/validation_C_summary_stats.pdf", p_c, width = 12, height = 6)

# ---- D: marginal distributions of six genes, bulk vs pseudo-bulk ----
genes_sel <- names(sort(apply(log_bulk, 1, sd), decreasing = TRUE))[1:6]
ks <- data.frame(gene = genes_sel, p_value = sapply(genes_sel, function(g)
  suppressWarnings(ks.test(log_pseudo[g, ], log_bulk[g, ])$p.value)))
ks$p_adj <- p.adjust(ks$p_value, "BH")
write.csv(ks, "results/validation_D_ks_tests.csv", row.names = FALSE)

plot_gene <- function(g) {
  df <- data.frame(value = c(log_pseudo[g, ], log_bulk[g, ]),
                   type  = rep(c("Synthetic pseudo-bulk", "Reference bulk"),
                               c(ncol(log_pseudo), ncol(log_bulk))))
  ggplot(df, aes(type, value, color = type)) +
    geom_violin(scale = "width") + geom_jitter(width = 0.1, size = 0.6) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.3,
             label = sprintf("p = %.2f", ks$p_value[ks$gene == g])) +
    labs(title = g, x = NULL, y = NULL) +
    theme_bw() + theme(axis.text.x = element_blank(), panel.grid = element_blank())
}
p_d <- wrap_plots(lapply(genes_sel, plot_gene), ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")
ggsave("figures/validation_D_marginals.pdf", p_d, width = 8, height = 6)

print(mlisi)
