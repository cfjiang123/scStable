#!/usr/bin/env Rscript
#
# Example script for running scStable
#
# This script demonstrates how to use the scStable package
# to generate synthetic single-cell RNA-seq replicates.
#

# =============================================================================
# 1. SETUP
# =============================================================================

# Load the package (after installation)
library(scStable)

# Set paths - MODIFY THESE TO MATCH YOUR SETUP
data_dir <- "data/"                          # Directory containing your data
output_dir <- "output/"                       # Directory for outputs
model_dir <- file.path(output_dir, "model")  # scDesign3 model outputs
replicate_dir <- file.path(output_dir, "replicates")  # Synthetic replicates

# Create output directories if they don't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(replicate_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

# Load your bulk RNA-seq matrix
# Format: genes (rows) x samples (columns), with gene names as rownames
bulk_matrix <- readRDS(file.path(data_dir, "bulk_matrix.rds"))
cat("Loaded bulk matrix:", nrow(bulk_matrix), "genes x", ncol(bulk_matrix), "samples\n")

# Load your single-cell RNA-seq count matrix
# Format: genes (rows) x cells (columns), with gene names as rownames
sc_matrix <- readRDS(file.path(data_dir, "sc_matrix.rds"))
cat("Loaded sc matrix:", nrow(sc_matrix), "genes x", ncol(sc_matrix), "cells\n")

# =============================================================================
# 3. RUN PIPELINE (Step by Step)
# =============================================================================

cat("\n========================================\n")
cat("Running scStable pipeline\n")
cat("========================================\n\n")

# Step 1: Preprocess and filter data
cat("Step 1: Preprocessing data...\n")
prep_result <- synthreplicate_prep(
    bulkRNA_matrix = bulk_matrix,
    scRNA_matrix = sc_matrix,
    use.pc = TRUE,
    number.pc = 10,           # Number of PCs to compute
    number.gene = 1500,       # Number of highly variable genes
    bulk_freq = 0.2,          # Max zero frequency for bulk genes
    sc_freq = 0.95            # Max zero frequency for sc genes
)
cat("  Selected", length(prep_result$hvg), "highly variable genes\n")

# Step 2: Fit bulk distribution
cat("\nStep 2: Fitting bulk distribution...\n")
bulk_fit <- fit_bulk(
    bulkRNA_matrix = prep_result$bulk,
    n_cores = 10,             # Adjust based on your system
    p_val_threshold = 0.5,
    normal_test = "auto"
)
cat("  Bulk distribution fitting completed\n")

# Step 3: Fit scDesign3 model
cat("\nStep 3: Fitting scDesign3 model...\n")
sc_model <- scDesign3_fit(
    scRNA_matrix = prep_result$sc,
    top_pcs = prep_result$pca,
    Cell_label = NULL,        # Optional: provide cell type labels
    save_dir = model_dir,
    use.option = 2,           # Modeling option (1, 2, or 3)
    n_cores_marginal = 10,    # Adjust based on your system
    family_copula = "nb"
)
cat("  scDesign3 model saved to:", model_dir, "\n")

# Step 4: Generate synthetic bulk samples
cat("\nStep 4: Generating synthetic bulk samples...\n")
synth_bulk <- synthreplicate_gen_bulk(
    bulkRNA_matrix = prep_result$bulk,
    mu = bulk_fit$mu,
    cov = bulk_fit$cov,
    optimal_c = bulk_fit$optimal_c,
    scRNA_matrix = prep_result$sc,
    use.cor = 3,              # Correlation method (1, 2, or 3)
    number.replicate = 50     # Number of replicates to generate
)
cat("  Generated", ncol(synth_bulk$sampling_bulk), "synthetic bulk samples\n")

# Step 5: Generate synthetic single-cell replicates
cat("\nStep 5: Generating synthetic single-cell replicates...\n")
synthreplicate_gen_sc(
    bulkRNA_matrix = prep_result$bulk,
    bulk_synth = synth_bulk$sampling_bulk,
    mu = bulk_fit$mu,
    d = synth_bulk$d,
    optimal_c = bulk_fit$optimal_c,
    scRNA_matrix = prep_result$sc,
    scRNA_list = sc_model,
    match.option = 2,         # Matching method (1 or 2)
    scaling_factor = 1,
    n_cores = 1,
    save.dir = replicate_dir
)
cat("  Synthetic replicates saved to:", replicate_dir, "\n")

# =============================================================================
# 4. VERIFY OUTPUT
# =============================================================================

cat("\n========================================\n")
cat("Pipeline completed!\n")
cat("========================================\n\n")

# List generated files
files <- list.files(replicate_dir, pattern = "replicate.*\\.csv$")
cat("Generated", length(files), "replicate files:\n")
for (f in head(files, 5)) {
    cat("  -", f, "\n")
}
if (length(files) > 5) {
    cat("  ... and", length(files) - 5, "more\n")
}

# =============================================================================
# ALTERNATIVE: Use the all-in-one wrapper (NO paired bulk required)
# =============================================================================

# If you do NOT have a paired bulk sample for your single-cell data, scStable
# can borrow a public bulk reference from the SAME TISSUE (e.g. GTEx) as a
# surrogate anchor and run the whole pipeline in one call:
#
# result <- synthreplicate_from_tissue(
#     tissue_name = "Liver",                    # tissue matching your sc data
#     scRNA_matrix = sc_matrix,
#     gtex_data_dir = "path/to/gtex/data/",    # directory with <tissue>.RDS files
#     save_dir = model_dir,
#     replicate_dir = replicate_dir,
#     number.gene = 1500,
#     number.pc = 10,
#     number.replicate = 100
# )

cat("\nDone!\n")
