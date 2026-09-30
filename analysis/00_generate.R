# Generate the inputs used by the figure scripts in analysis/.
# Single-sample scRNA-seq reference: Zheng4 (Zhengmix4eq, DuoClustering2018).
# Multi-sample bulk RNA-seq reference: GTEx whole blood (bulk-reference-free mode).
#
# Outputs (run from the repository root):
#   data/gtex_whole_blood_counts.rds    GTEx whole-blood counts, genes (symbols) x samples
#   data/cell_types.rds                 named vector of Zheng4 cell types (names = cell barcodes)
#   data/synthetic/fit.rds              list(prep, bulk_fit, gen_bulk):
#                                         prep     synthreplicate_prep() output
#                                         bulk_fit fit_bulk() on the GTEx bulk reference
#                                         gen_bulk synthreplicate_gen_bulk() output
#   data/synthetic/model/               scDesign3_fit() objects
#   data/synthetic/lambda_<l>/replicate<r>.csv   synthetic samples for each scale factor
#
# GTEx files (GTEx Portal, open access, https://gtexportal.org/home/downloads):
#   gene read counts   GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
#   sample attributes  GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
suppressPackageStartupMessages({
  library(scStable)
  library(SummarizedExperiment)
  library(data.table)
})

GTEX_GCT   <- "data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
GTEX_ATTR  <- "data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
TISSUE     <- "Whole Blood"
N_GENE     <- 1500   # HVGs; the clustering analysis in the manuscript uses 200
N_PC       <- 10
R          <- 100
LAMBDAS    <- c(0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2)
N_CORES    <- 10
OUT        <- "data/synthetic"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- GTEx whole-blood counts, Ensembl IDs -> gene symbols ----
attr    <- fread(GTEX_ATTR, select = c("SAMPID", "SMTSD"))
header  <- names(fread(GTEX_GCT, skip = 2, nrows = 0))
samples <- intersect(attr$SAMPID[attr$SMTSD == TISSUE], header)
gct     <- fread(GTEX_GCT, skip = 2, select = c("Name", samples))
bulk    <- as.matrix(gct[, samples, with = FALSE])

sym  <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
          keys = sub("\\..*$", "", gct$Name), keytype = "ENSEMBL", column = "SYMBOL"))
keep <- !is.na(sym) & !duplicated(sym)
bulk <- bulk[keep, , drop = FALSE]
rownames(bulk) <- sym[keep]
saveRDS(bulk, "data/gtex_whole_blood_counts.rds")
cat("GTEx", TISSUE, ":", nrow(bulk), "genes x", ncol(bulk), "samples\n")

# ---- Zheng4 scRNA-seq reference and cell types ----
sce <- DuoClustering2018::sce_filteredExpr10_Zhengmix4eq(metadata = FALSE)
sc  <- as.matrix(assay(sce, "counts"))
sym_sc <- rowData(sce)$symbol
bad    <- is.na(sym_sc) | sym_sc == ""
sym_sc[bad] <- rownames(sc)[bad]
rownames(sc) <- make.unique(sym_sc)
saveRDS(setNames(as.character(sce$phenoid), colnames(sc)), "data/cell_types.rds")
cat("Zheng4:", nrow(sc), "genes x", ncol(sc), "cells\n")

# ---- Steps 1-3: fit the scRNA-seq and bulk models, draw synthetic bulk ----
set.seed(1)
prep     <- synthreplicate_prep(bulk, sc, number.pc = N_PC, number.gene = N_GENE)
bulk_fit <- fit_bulk(prep$bulk, n_cores = N_CORES)
sc_fit   <- scDesign3_fit(prep$sc, prep$pca, save_dir = file.path(OUT, "model"),
                          n_cores_marginal = N_CORES)
gen_bulk <- synthreplicate_gen_bulk(prep$bulk, bulk_fit$mu, bulk_fit$cov, bulk_fit$optimal_c,
                                    prep$sc, number.replicate = R)
saveRDS(list(prep = prep, bulk_fit = bulk_fit, gen_bulk = gen_bulk), file.path(OUT, "fit.rds"))

# ---- Steps 3-4: synthetic scRNA-seq samples for each scale factor ----
for (l in LAMBDAS) {
  dir_l <- file.path(OUT, paste0("lambda_", l))
  cat("lambda =", l, "->", dir_l, "\n")
  synthreplicate_gen_sc(prep$bulk, gen_bulk$sampling_bulk, bulk_fit$mu, gen_bulk$d,
                        bulk_fit$optimal_c, prep$sc, sc_fit,
                        scaling_factor = l, save.dir = dir_l)
}
