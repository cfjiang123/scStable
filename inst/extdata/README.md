# Example Data Directory

This directory (`inst/extdata/`) is where you should place your data files.

## Required Data Format

### Bulk RNA-seq Matrix (`bulk_example.rds`)

Your bulk RNA-seq matrix should be an R matrix object with:
- **Rows**: Gene names (as rownames)
- **Columns**: Sample IDs
- **Values**: Normalized expression values (e.g., TPM, FPKM)

### Single-Cell RNA-seq Matrix (`sc_example.rds`)

Your scRNA-seq matrix should be a sparse or dense matrix with:
- **Rows**: Gene names (as rownames)
- **Columns**: Cell barcodes
- **Values**: Raw counts (integers)

## How to Add Your Data

1. Prepare your data in the format described above
2. Save as RDS files:
```r
saveRDS(bulk_matrix, "inst/extdata/bulk_example.rds")
saveRDS(sc_matrix, "inst/extdata/sc_example.rds")
```
