# Bovine Multi-Omics Analysis Scripts

Research scripts for single-nucleus RNA sequencing, single-nucleus ATAC sequencing,
and Visium spatial transcriptomics in bovine tissues.

This repository is a collection of analysis scripts. Some files define reusable functions;
others require prepared Seurat objects in an active R session. Raw data, complete reference
resources, and a validated environment are not included, so the repository does not yet
establish end-to-end reproducibility of the study results.

## Repository Guide

| File | Purpose | Usage |
| --- | --- | --- |
| `snRNA/Merge_Single-Cell_h5_Datasets.r` | Merge 10x H5 datasets | Run with `Rscript` |
| `snRNA/QC_and_Basic_Analysis.r` | DoubletFinder, decontamination, dimensionality reduction, and enrichment | Source the script and call its functions |
| `snRNA/Celltype_Identification_and_Visualization.r` | Cell-type annotation and visualization | Source the script and call its functions |
| `snRNA/CellChat_Inter.r` | Compare cell-cell communication | Source the script and call its functions |
| `snRNA/Monocle3_pipeline.r` | Trajectory analysis and pseudotime plots | Source the script and call its functions |
| `snATAC/Multi-Sample_scATAC-seq_Peak_Call.r` | Construct consensus peaks and quantify fragments | Configure inputs before running |
| `snATAC/QC_basic_analysis.r` | ATAC QC, integration, and gene activity | Source the script and call its functions |
| `snATAC/integration.r` | Transfer labels from RNA to ATAC data | Configure environment variables before running |
| `snATAC/Findmotif.r` | Motif analysis, chromVAR, and footprinting | Requires an active `sc` object and sequence names matching bosTau9 |
| `snATAC/heatmap_RNA&ATAC.r` | TF activity and expression heatmaps | Requires an active `sc` object |
| `ST/basic_analysis.py` | Visium QC, clustering, and spatial statistics | Run from the Python command line |
| `ST/tangram_single.py` | Map single cells to spatial spots | Run from the Python command line |
| `snATAC/scenicplus.py` | Legacy SCENIC+ API example | Unvalidated; see the known limitations |

## Environment

Spatial analysis requires `scanpy`, `squidpy`, `matplotlib`, `igraph`, and `leidenalg`.
Mapping additionally requires `tangram-sc` and `torch`. Use a separate environment for
SCENIC+ and prepare its inputs according to the official tutorial for the selected version.
See [DEPENDENCIES.md](docs/DEPENDENCIES.md) for the R package inventory.

No dependency lockfile from a successfully validated analysis environment is available.
After a successful run, record `pip freeze` and R's `sessionInfo()`. An untested list of
package versions should not be presented as a reproducible environment.

## Running the Scripts

Run the following commands from the repository root. Input paths are resolved relative
to the directory from which the command is invoked; `--work_dir` specifies the output
directory. Outputs use fixed filenames, so use a new output directory for each run to
preserve previous results.

```bash
Rscript snRNA/Merge_Single-Cell_h5_Datasets.r data/h5 results/merged.rds

python3 ST/basic_analysis.py \
  --adata_path data/visium/sample1 --work_dir results/spatial/sample1

python3 ST/tangram_single.py \
  --sc_path data/reference.h5ad \
  --st_path results/spatial/sample1/processed_data.h5ad \
  --annotation_key celltype --work_dir results/tangram/sample1 --seed 0
```

The Visium library ID is inferred when exactly one library is present. If multiple
libraries are present, select one explicitly with `--library_id`.

Spatial preprocessing preserves counts in `layers["counts"]` and log-normalized expression
in `raw`; differential expression uses `raw`. Tangram preferentially reads the counts
layer and saves the mapping matrix, a CSV weight table, and `obsm["tangram_ct_pred"]`.
These predictions are mapping weights and should not be interpreted directly as absolute
cell counts.

Configure ATAC inputs through environment variables:

```bash
export ATLAS_ATAC_DATA_DIR="$PWD/data/atac"
export ATLAS_RNA_REFERENCE="$PWD/data/reference.rds"
export ATLAS_ATAC_QUERY="$PWD/data/query.rds"
```

The peak-processing script's `DATA_PATHS` and `SAMPLE_IDS` must still be configured for
your samples. Genome assembly, chromosome naming, annotation, and fragment coordinates
must be consistent.

## Validation and Contributions

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q ST snATAC
```

The unit tests check filesystem handling and library selection; they do not validate
the scientific analysis. Changes to QC, species annotation, normalization, differential
analysis, or model parameters require comparisons on representative datasets.

See [KNOWN_LIMITATIONS.md](docs/KNOWN_LIMITATIONS.md) for unresolved issues and
[DATA_AND_LICENSE.md](docs/DATA_AND_LICENSE.md) for data provenance and licensing status.
