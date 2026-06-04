###############################################################################
# scStable - step 1: filter the bulk RNA-seq and scRNA-seq
# input : bulk RNA-seq matrix (with names, normalized),
#         scRNA-seq matrix (with names),
#         scRNA-seq cell type (optional, not used)
# output: filtered bulk RNA-seq matrix (with names),
#         filtered scRNA-seq matrix (with names),
#         top PCs (with cell names)
###############################################################################

#' Preprocess bulk and single-cell RNA-seq data for scStable
#'
#' \code{synthreplicate_prep} is the first step of the \pkg{scStable} workflow.
#' It harmonises a bulk RNA-seq matrix and a reference scRNA-seq matrix to a
#' common, informative gene set: genes are filtered by zero-count frequency in
#' both modalities, highly variable genes (HVGs) are selected from the
#' single-cell data with \pkg{Seurat}, and (optionally) the top principal
#' components of the single-cell HVG space are computed for downstream
#' modelling.
#'
#' @param bulkRNA_matrix Numeric matrix of (normalized) bulk RNA-seq counts with
#'   gene names as \code{rownames} and samples as columns.
#' @param scRNA_matrix Numeric or sparse matrix of single-cell RNA-seq counts
#'   with gene names as \code{rownames} and cells as columns.
#' @param use.pc Logical; if \code{TRUE} (default) compute and return the top
#'   principal components of the single-cell HVG space.
#' @param number.pc Integer; number of principal components to compute / return.
#' @param number.gene Integer; number of highly variable genes to select.
#' @param bulk_freq Numeric in \[0, 1\]; keep bulk genes whose proportion of
#'   minimum-value (e.g. zero) entries is below this threshold.
#' @param sc_freq Numeric in \[0, 1\]; keep single-cell genes whose proportion
#'   of zero entries is below this threshold.
#' @param min.cells Integer passed to \code{Seurat::CreateSeuratObject}.
#' @param select_genes Optional character vector of gene names to use instead of
#'   automatically selected HVGs.
#'
#' @return A list with elements \code{bulk} (filtered bulk matrix on HVGs),
#'   \code{sc} (filtered single-cell matrix on HVGs), \code{hvg} (selected gene
#'   names), \code{pca} (PC embeddings or \code{NULL}) and \code{seurat} (the
#'   \code{Seurat} object).
#'
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures
#'   VariableFeatures ScaleData RunPCA Embeddings
#' @export
synthreplicate_prep <- function(bulkRNA_matrix,
                                scRNA_matrix,
                                use.pc       = TRUE,
                                number.pc    = 20,
                                number.gene  = 2000,
                                bulk_freq    = 0.2,
                                sc_freq = 0.95,
                                min.cells    = 2,
                                select_genes = NULL) {
  # 0. Check inputs: both matrices must have rownames = gene names
  if (is.null(rownames(bulkRNA_matrix)) ||
      is.null(rownames(scRNA_matrix))) {
    stop("bulkRNA_matrix or scRNA_matrix must contain gene names as rownames")
  }

  # 1. Filter bulk RNA-seq genes by zero-count frequency
  #    Keep genes with proportion of zeros < bulk_freq
  min_freq <- apply(bulkRNA_matrix, 1, function(x) sum(x == min(x))/length(x))
  indices <- which(min_freq < bulk_freq) # bulk_freq = 0.2, default
  keep_genes = rownames(bulkRNA_matrix)[indices]

  min_freq <- rowMeans(scRNA_matrix == 0, na.rm = TRUE)
  indices <- which(min_freq < sc_freq)
  keep_genes2 = rownames(scRNA_matrix)[indices]

  common_genes <- intersect(keep_genes, keep_genes2)
  if (length(common_genes) == 0) {
    stop("No common genes remain after filtering by bulk_freq")
  }
  bulkRNA_matrix_c <- bulkRNA_matrix[common_genes, , drop = FALSE]
  scRNA_matrix_c   <- scRNA_matrix[common_genes, , drop = FALSE]

  if(!is.null(select_genes)){
    bulkRNA_matrix_c <- bulkRNA_matrix[select_genes, , drop = FALSE]
    scRNA_matrix_c   <- scRNA_matrix[select_genes, , drop = FALSE]
  }


  # 2. Create Seurat object and normalize single-cell data
  scS <- CreateSeuratObject(
    counts       = scRNA_matrix_c,
    project      = "scS",
    min.cells    = min.cells,
    min.features = 0
  )
  scS <- NormalizeData(scS)

  # 3. Select highly variable genes (HVGs)

    max_genes <- nrow(scS)
    if (number.gene > max_genes) {
      warning(sprintf(
        "number.gene (%d) exceeds available genes (%d); resetting to %d.",
        number.gene, max_genes, max_genes
      ))
      number.gene <- max_genes
    }
    scS <- FindVariableFeatures(
      object          = scS,
      selection.method = "vst",
      nfeatures        = number.gene
    )
  if(is.null(select_genes)){
    hvg <- VariableFeatures(scS)
  }
  else{
    hvg <- select_genes
  }


  # 4. (Optional) Scale HVGs and run PCA
  if (use.pc) {
    scS <- ScaleData(scS, features = hvg)
    scS <- RunPCA(scS, features = hvg, npcs = number.pc)
    pca_embed <- Embeddings(scS, "pca")
    n_pc_avail <- ncol(pca_embed)
    pca_mat    <- pca_embed[, seq_len(min(number.pc, n_pc_avail)), drop = FALSE]
  } else {
    pca_mat <- NULL
  }

  # 5. Subset bulk & single-cell matrices to HVGs
  bulk_out <- bulkRNA_matrix_c[hvg, , drop = FALSE]
  sc_out   <- as.matrix(scRNA_matrix_c[hvg, , drop = FALSE])

  # 6. Return a list: filtered bulk, filtered sc, HVG names, PCA embeddings, and Seurat object
  list(
    bulk   = bulk_out,
    sc     = sc_out,
    hvg    = hvg,
    pca    = pca_mat,
    seurat = scS
  )
}
