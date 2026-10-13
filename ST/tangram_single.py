"""Map annotated single cells to Visium spots with Tangram.

Prediction weights are saved in obsm['tangram_ct_pred']; they are not calibrated
cell counts. Inputs must contain nonnegative expression, not scaled values.
"""

import argparse
import logging

from pipeline_utils import input_path, output_directory, select_library

LOGGER = logging.getLogger(__name__)


def main(work_dir, sc_path, st_path, annotation_key, library_id=None, seed=0):
    reference_path = input_path(sc_path)
    spatial_path = input_path(st_path)
    output = output_directory(work_dir)

    import scanpy as sc
    import tangram as tg
    import torch
    from matplotlib import pyplot as plt

    reference = sc.read_h5ad(reference_path)
    spatial = sc.read_h5ad(spatial_path)
    if annotation_key not in reference.obs:
        raise ValueError(f"Reference obs lacks annotation column: {annotation_key}")
    if reference.obs[annotation_key].isna().any():
        raise ValueError("Cell annotations contain missing values.")
    library = select_library(spatial, library_id)
    for adata in (reference, spatial):
        adata.var_names_make_unique()
        if "counts" in adata.layers:
            adata.X = adata.layers["counts"].copy()
        if adata.n_obs == 0 or adata.n_vars == 0:
            raise ValueError("Input expression matrix is empty.")
        if adata.X.min() < 0:
            raise ValueError("Tangram needs nonnegative expression; provide a counts layer.")

    # Preserve the original intent: map individual cells, then aggregate labels.
    # All shared expressed genes are used; no unused HVG selection is performed.
    tg.pp_adatas(reference, spatial, genes=None)
    if not reference.uns["training_genes"]:
        raise ValueError("Reference and spatial data share no expressed training genes.")
    device = "cuda:0" if torch.cuda.is_available() else "cpu"
    LOGGER.info("Mapping cells on %s", device)
    mapping = tg.map_cells_to_space(
        reference, spatial, mode="cells", device=device, random_state=seed
    )
    tg.project_cell_annotations(mapping, spatial, annotation=annotation_key)
    predictions = spatial.obsm["tangram_ct_pred"]
    predictions.to_csv(output / "tangram_celltype_weights.csv")

    # Plot on a copy so reference labels cannot overwrite existing spot metadata.
    plot_data = spatial.copy()
    plot_keys = []
    for index, cell_type in enumerate(predictions.columns[:6]):
        key = f"tangram_prediction_{index}"
        plot_data.obs[key] = predictions[cell_type].reindex(plot_data.obs_names)
        plot_keys.append(key)
    sc.pl.spatial(
        plot_data, color=plot_keys, title=[str(x) for x in predictions.columns[:6]],
        cmap="viridis", library_id=library, show=False,
    )
    plt.savefig(output / "tangram_celltype_spatial_density.pdf", bbox_inches="tight")
    plt.close("all")
    spatial.uns["tangram_run"] = {"seed": seed, "annotation_key": annotation_key, "mode": "cells"}
    mapping.write(output / "tangram_mapping.h5ad")
    spatial.write(output / "spatial_with_tangram.h5ad")
    LOGGER.info("Saved mapping and predictions to %s", output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work_dir", required=True, help="Output directory")
    parser.add_argument("--sc_path", required=True, help="Reference h5ad file")
    parser.add_argument("--st_path", required=True, help="Spatial h5ad file")
    parser.add_argument("--annotation_key", default="celltype")
    parser.add_argument("--library_id", help="Visium library ID; inferred for one library")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    main(**vars(args))
