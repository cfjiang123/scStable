# Stable DE figure: DE between two cell types on the reference and on each
# lambda = 1 synthetic sample; q-value stability across synthetic samples.
#   A  95% quantile intervals of q for example genes; density of interval length
#   B  reference DE rank vs stability rank (interval length) for six tests
#   stable DE genes: reference q < 0.05 and upper interval end < 0.05
# (The manuscript figure also shows KEGG enrichment of stable vs other DE genes;
#  this needs real gene symbols and is omitted here.)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

TESTS <- c("wilcox", "negbinom", "bimod", "t", "LR", "poisson")
PAIR  <- c("typeA", "typeB")
ALPHA <- 0.05
CONF  <- 0.95

ref       <- readRDS("data/synthetic/fit.rds")$prep$sc
cell_type <- readRDS("data/cell_types.rds")[colnames(ref)]
sample_files <- list.files("data/synthetic/lambda_1", "^replicate\\d+\\.csv$", full.names = TRUE)
samples <- lapply(sample_files, function(f)
  as.matrix(read.table(f, sep = "\t", header = TRUE, check.names = FALSE)))

# q-values (Seurat Bonferroni-adjusted p) of each test; genes filtered by FindMarkers are NA
find_de <- function(counts) {
  seu <- CreateSeuratObject(as(counts, "CsparseMatrix"), min.cells = 0, min.features = 0)
  seu <- NormalizeData(seu, verbose = FALSE)
  Idents(seu) <- factor(cell_type[colnames(seu)])
  q <- matrix(NA_real_, nrow(seu), length(TESTS), dimnames = list(rownames(seu), TESTS))
  for (m in TESTS) {
    res <- FindMarkers(seu, ident.1 = PAIR[1], ident.2 = PAIR[2], test.use = m,
                       logfc.threshold = 0.1, min.pct = 0.1, verbose = FALSE)
    q[rownames(res), m] <- res$p_val_adj
  }
  q
}
q_ref <- find_de(ref)
q_syn <- lapply(samples, find_de)

# per test: quantile interval of q across synthetic samples
intervals <- lapply(TESTS, function(m) {
  v  <- sapply(q_syn, function(q) q[rownames(q_ref), m])
  ci <- t(apply(v, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) >= 2) quantile(x, c((1 - CONF) / 2, 1 - (1 - CONF) / 2), names = FALSE) else c(NA, NA)
  }))
  data.frame(gene = rownames(q_ref), ref_q = q_ref[, m], lower = ci[, 1], upper = ci[, 2],
             L_q = ci[, 2] - ci[, 1], row.names = rownames(q_ref))
})
names(intervals) <- TESTS

# ---- stable vs other DE genes (Wilcoxon) ----
w <- intervals$wilcox
w <- w[!is.na(w$ref_q) & w$ref_q < ALPHA, ]
w$stable <- !is.na(w$upper) & w$upper < ALPHA
write.csv(w[order(w$ref_q), ], "results/stable_de_wilcox.csv", row.names = FALSE)
cat("Reference DE genes:", nrow(w), " stable:", sum(w$stable), " other:", sum(!w$stable), "\n")

# ---- B: DE rank vs stability rank, genes significant with a proper interval under every test ----
ok <- Reduce(intersect, lapply(intervals, function(d)
  d$gene[!is.na(d$ref_q) & d$ref_q < ALPHA & !is.na(d$L_q) & d$upper < 1]))
rank_df <- do.call(rbind, lapply(TESTS, function(m) {
  d <- intervals[[m]][ok, ]
  data.frame(test = m, gene = ok, de_rank = rank(d$ref_q, ties.method = "min"),
             stability_rank = rank(d$L_q, ties.method = "min"))
}))
rank_df$test <- factor(rank_df$test, levels = TESTS)
write.csv(rank_df, "results/stable_de_ranks.csv", row.names = FALSE)
p_b <- ggplot(rank_df, aes(de_rank, stability_rank)) +
  geom_point(alpha = 0.6) + facet_wrap(~test, nrow = 2) +
  labs(x = "Reference rank", y = "Stability rank") + theme_bw()
ggsave("figures/stable_de_B_rank.pdf", p_b, width = 10, height = 6)

# ---- A: example intervals (4 shortest, 4 widest) and interval-length density ----
eps <- 1e-6
wi  <- w[!is.na(w$L_q), ]
wi  <- wi[order(wi$L_q), ]
ex  <- rbind(transform(head(wi, 4), group = "Stable (short interval)"),
             transform(tail(wi, 4), group = "Unstable (wide interval)"))
ex$gene <- factor(ex$gene, levels = ex$gene)
p_int <- ggplot(ex, aes(gene, color = group)) +
  geom_errorbar(aes(ymin = pmax(lower, eps), ymax = pmax(upper, eps)), width = 0.3, linewidth = 1) +
  scale_y_log10() + labs(x = NULL, y = "q-value (log10)", color = NULL) +
  theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
p_den <- ggplot(wi, aes(log10(L_q + eps))) + geom_density(linewidth = 1) +
  labs(x = "log10 interval length", y = "Density") + theme_classic()
ggsave("figures/stable_de_A_intervals.pdf", p_int / p_den, width = 5, height = 6)
