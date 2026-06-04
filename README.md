# scStable

<!-- badges: start -->
<!-- badges: end -->

**scStable** generates synthetic single-cell RNA-seq (scRNA-seq) replicates that
are anchored to bulk RNA-seq data. It fits a multivariate model to a reference
bulk RNA-seq matrix, samples new synthetic bulk profiles, and propagates the
induced per-gene variation into newly simulated single-cell count matrices
(via [scDesign3](https://github.com/SONGDONGYUAN1994/scDesign3)), producing
biologically plausible synthetic replicates for benchmarking and method-stability
evaluation.

> This package implements the **scStable** method described in our manuscript
> submitted to *Genome Biology*. If you use scStable, please cite the paper
> (citation details to be added upon acceptance).

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
remotes::install_github("yourusername/scStable")
```

## Workflow

scStable runs in five steps. Each step has a dedicated function, and the wrapper
`synthreplicate_from_tissue()` chains them together.

| Step | Function | Purpose |
|------|----------|---------|
| 1 | `synthreplicate_prep()` | Filter/harmonise bulk + sc matrices, select HVGs, compute PCs |
| 2 | `fit_bulk()` | Fit the bulk RNA-seq distribution (mean, covariance, per-gene offset) |
| 3 | `scDesign3_fit()` | Fit the single-cell generative model (scDesign3) |
| 4 | `synthreplicate_gen_bulk()` | Sample synthetic bulk replicates |
| 5 | `synthreplicate_gen_sc()` | Generate and save synthetic single-cell replicates |

## Quick start

### Option A — step by step

```r
library(scStable)

# bulkRNA_matrix : genes x samples (rownames = gene symbols)
# scRNA_matrix   : genes x cells   (rownames = gene symbols)

prep <- synthreplicate_prep(
  bulkRNA_matrix = bulkRNA_matrix,
  scRNA_matrix   = scRNA_matrix,
  number.pc      = 10,
  number.gene    = 1500
)

bulk_fit <- fit_bulk(
  bulkRNA_matrix = prep$bulk,
  n_cores        = 20
)

sc_fit <- scDesign3_fit(
  scRNA_matrix = prep$sc,
  top_pcs      = prep$pca,
  save_dir     = "scStable_model",
  use.option   = 2
)

bulk_synth <- synthreplicate_gen_bulk(
  bulkRNA_matrix   = prep$bulk,
  mu               = bulk_fit$mu,
  cov              = bulk_fit$cov,
  optimal_c        = bulk_fit$optimal_c,
  scRNA_matrix     = prep$sc,
  use.cor          = 3,
  number.replicate = 100
)

synthreplicate_gen_sc(
  bulkRNA_matrix = prep$bulk,
  bulk_synth     = bulk_synth$sampling_bulk,
  mu             = bulk_fit$mu,
  d              = bulk_synth$d,
  optimal_c      = bulk_fit$optimal_c,
  scRNA_matrix   = prep$sc,
  scRNA_list     = sc_fit,
  match.option   = 2,
  save.dir       = "scStable_replicates"
)
```

Each synthetic replicate is written to `scStable_replicates/replicate<i>.csv`.

### Option B — one call (GTEx tissue wrapper)

```r
res <- synthreplicate_from_tissue(
  tissue_name      = "Liver",
  scRNA_matrix     = scRNA_matrix,
  gtex_data_dir    = "path/to/GTEx/tissue_data/",
  save_dir         = "scStable_model",
  replicate_dir    = "scStable_replicates",
  number.replicate = 100
)
```

## License

MIT © the scStable authors.
