# scStable

<!-- badges: start -->
<!-- badges: end -->

**scStable** generates multiple synthetic scRNA-seq samples from a single-sample
scRNA-seq. All synthetic scRNA-seq samples contain the same cells as the 
reference. They differ in gene expression by between-sample variation learned 
from multi-sample bulk RNA-seq data. The synthetic samples support stability 
assessment of downstream discoveries (such as DE genes and clusters) and 
stability-driven selection of analysis methods and hyperparameters.

scStable has two modes:

* **Bulk-reference mode**: supply a tissue- and condition-matched multi-sample
  *bulk RNA-seq reference*. Its between-sample variation is estimated with
  `fit_bulk()` (Option A below).
* **Bulk-reference-free mode**: no matched bulk reference is available.
  scStable uses pre-estimated, tissue-specific parameters, for example learned
  from GTEx, via `synthreplicate_from_tissue()` (Option B below).

> This package implements the **scStable** method described in our manuscript.
> If you use scStable, please cite the paper (citation details to be added upon acceptance).

## Installation

scStable depends on Bioconductor packages and on `scDesign3`. Install the
dependencies first, then install scStable from GitHub:

```r
# 1. Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment", "SummarizedExperiment", "S4Vectors", "BiocParallel"
))

# 2. scDesign3 (Bioconductor)
BiocManager::install("scDesign3")

# 3. CRAN dependencies
install.packages(c(
  "Seurat", "foreach", "doParallel", "nortest",
  "tmvtnorm", "truncnorm", "corpcor", "MASS"
))

# 4. scStable
# install.packages("remotes")
remotes::install_github("cfjiang123/scStable")
```

## Workflow

| Step | Function | Purpose |
|------|----------|---------|
| prep | `synthreplicate_prep()` | Shared genes, highly variable genes, top PCs of the scRNA-seq reference |
| 1 | `scDesign3_fit()` | Fit a cell-label-free scRNA-seq generative model (NB marginals, Gaussian copula, PCs as covariates) |
| 2 | `fit_bulk()` | Estimate bulk-derived between-sample variation (mean, covariance, pseudo-counts) |
| 3 | `synthreplicate_gen_bulk()` + `synthreplicate_gen_sc()` | Map the variation onto the scRNA-seq model, with scale factor λ (`scaling_factor`) |
| 4 | `synthreplicate_gen_sc()` | Generate the synthetic scRNA-seq samples |

The scale factor λ sets the magnitude of injected variation: `0` means none,
`1` means realistic bulk-derived variation, and values above `1` amplify it for
stress-testing.


## Quick start

### Option A: bulk-reference mode

```r
library(scStable)

# scRNA_matrix   : genes x cells   counts (rownames = gene symbols)
# bulkRNA_matrix : genes x samples counts (rownames = gene symbols)

prep <- synthreplicate_prep(bulkRNA_matrix, scRNA_matrix, number.gene = 1500)

sc_fit   <- scDesign3_fit(prep$sc, prep$pca, save_dir = "scStable_model")   # Step 1
bulk_fit <- fit_bulk(prep$bulk, n_cores = 4)                                # Step 2

bulk_synth <- synthreplicate_gen_bulk(                                      # Step 3
  prep$bulk, bulk_fit$mu, bulk_fit$cov, bulk_fit$optimal_c, prep$sc,
  number.replicate = 100
)
synthreplicate_gen_sc(                                                      # Steps 3-4
  prep$bulk, bulk_synth$sampling_bulk, bulk_fit$mu, bulk_synth$d,
  bulk_fit$optimal_c, prep$sc, sc_fit,
  scaling_factor = 1, save.dir = "scStable_samples"
)
```

Each synthetic sample is written to `scStable_samples/replicate<r>.csv`
(tab-separated, genes x cells). By default the synthetic bulk samples are drawn
with the full covariance if the bulk reference has at least 30 samples, and
with a diagonal covariance otherwise (`use.cor`).

### Option B: bulk-reference-free mode

Estimate the tissue-specific parameters once from a large public bulk resource
and save them. Then reuse them for any scRNA-seq sample of that tissue:

```r
# once per tissue, e.g. GTEx counts for the tissue (genes x samples)
saveRDS(fit_bulk(gtex_tissue_counts, n_cores = 8), "prostate_params.rds")

res <- synthreplicate_from_tissue(
  scRNA_matrix     = scRNA_matrix,
  bulk_params      = "prostate_params.rds",
  save_dir         = "scStable_model",
  replicate_dir    = "scStable_samples",
  number.replicate = 100
)
```

Alternatively, pass `tissue_name` and `gtex_data_dir` (a folder of per-tissue
`SummarizedExperiment` `.RDS` files). The parameters are then estimated on the
fly.

## Assessing stability

The script `inst/scripts/scStable_stability.R` computes the sample stability
metrics of the manuscript from the synthetic samples. Each function documents
its input and output in the script header.

```r
source(system.file("scripts", "scStable_stability.R", package = "scStable"))
samples <- list_synthetic_samples("scStable_samples")   # replicate<r>.csv paths
ref     <- prep$sc                                       # the scRNA-seq reference

# Stable DE genes
de <- stable_de(ref, samples, group = cell_type, ident.1 = "B", ident.2 = "T",
                test.use = "wilcox")
subset(de, stable)

# Per-cell entropy stability from clusterings of the synthetic samples
labels <- sapply(samples, function(f) seurat_louvain_grid(read_synthetic_sample(f),
                 n_pc = 10, k_nn = 20, resolution = 0.8)[, 1])
S <- entropy_stability(labels)

# Hyperparameter selection by stability 
sel <- stability_select(ref, samples, run_fun = seurat_louvain_grid)
head(sel$summary)
```

| Function | Input | Output |
|----------|-------|--------|
| `quantile_interval()` | features x R matrix of any scalar summary T | `lower`, `upper`, length `L` per feature |
| `stable_de()` | reference, samples, cell labels (or two conditions via `ref2`/`samples2`) | per-gene `ref_padj`, interval, `L`, `de`, `stable` |
| `entropy_stability()` | cells x R cluster labels | per-cell stability S<sub>i</sub> in [0, 1] |
| `stability_select()` | reference, samples, `run_fun` returning one label column per method/setting | ARI, NMI, Jaccard, FMI per specification, ranked |
| `stability_path()` | named list of samples per scale factor λ | the same metrics along λ |

Scripts reproducing the manuscript figures are in `analysis/` on GitHub (not
part of the installed package);

## License

MIT © Chengfeng Jiang.
