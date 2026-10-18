# Known Limitations and Validation Requirements

The repository cleanup improved readability, clarified inputs and outputs, and corrected
identifiable defects. End-to-end execution and result-consistency checks have not been
performed because the input data and an R/bioinformatics environment were unavailable.

- **SCENIC+ legacy example:** The example uses older modules and object-construction APIs.
  RNA/ATAC AnnData objects alone do not establish that the selected version's input
  requirements are met. Motif database paths are defined but are not passed into the
  workflow. Prepare cisTopic, motif-enrichment, and other required inputs using the
  official tutorial and revise the entry point accordingly. The cleanup only removed
  execution on import and global warning suppression, and made failures exit unsuccessfully.
- **Motif selection and sequence naming:** `Findmotif.r` now selects footprint targets
  using Signac's `motif.names`, replacing the previous `HasMotif`/scalar `||` filter.
  Chromosome names are no longer forcibly rewritten; sequence names incompatible with
  bosTau9 cause an error. Genome assembly, fragment coordinates, and the case of motif
  names still require validation on real data.
- **Peak processing:** `Multi-Sample_scATAC-seq_Peak_Call.r` retains the original sample
  examples. Removing the first row of `singlecell.csv`, restricting chromosomes to
  1–29/X, and fragment indexing methods must be checked against the actual files.
  BED-to-GRanges conversion now explicitly converts zero-based starts to one-based
  starts; this correction may change previous results.
- **DoubletFinder metadata:** Merging now preserves cell indices and rejects conflicting
  values. When multiple rounds produce different pANN scores or classifications, select
  the intended columns explicitly rather than silently concatenating them. The
  contamination-filtering step is now an explicit function and retains the merged columns.
- **H5 sample identifiers:** Merging uses complete filename stems as sample IDs rather
  than truncating names at underscores. If filenames include sequencing lanes, define
  the intended merging relationships using an experimental sample sheet.
- **Visium analysis parameters:** The mitochondrial gene list, liver marker panel,
  minimum counts, and Leiden resolution retain the original study settings. Review
  them for each tissue, gene identifier scheme, and sample. A fixed choice of 20 PCs
  is unsuitable for very small datasets. Preserving counts/raw and using log expression
  for differential testing changes the previous differential-expression workflow;
  review the resulting differences.
- **Multiple Visium libraries:** Spatial-neighbor construction still requires validation
  of `library_key` and the strategy for connections across samples.
- **Monocle3:** The script accesses internal slots. Missing or infinite pseudotime,
  zero variance, spline predictor values, and related edge cases require validation.
- **R plotting scripts:** Some depend on objects in the active session or predefined
  metadata. The README distinguishes function-based scripts from session-dependent snippets.
- **Reproducibility records:** A dependency lockfile from a successful analysis environment,
  raw-data provenance/accession numbers, reference-file checksums, and manuscript citation
  details are still missing.

## Official References

- [Tangram mapping API](https://github.com/broadinstitute/Tangram/blob/master/tangram/mapping_utils.py)
- [Tangram annotations API](https://github.com/broadinstitute/Tangram/blob/master/tangram/utils.py)
- [SCENIC+](https://github.com/aertslab/scenicplus)
- [Signac motif implementation](https://github.com/stuart-lab/signac/blob/master/R/motifs.R)
- [Signac footprint implementation](https://github.com/stuart-lab/signac/blob/master/R/footprinting.R)
