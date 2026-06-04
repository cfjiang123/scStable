# scStable 0.1.0

* Initial release of scStable.
* Implements the five-step scStable workflow for generating synthetic
  single-cell RNA-seq replicates anchored to bulk RNA-seq:
  `synthreplicate_prep()`, `fit_bulk()`, `scDesign3_fit()`,
  `synthreplicate_gen_bulk()` and `synthreplicate_gen_sc()`.
* Adds the end-to-end wrapper `synthreplicate_from_tissue()`, which supports
  running the full pipeline **without a paired bulk RNA-seq sample**: it uses a
  publicly available bulk reference from the same tissue (e.g. GTEx) as a
  surrogate anchor for the single-cell data.
