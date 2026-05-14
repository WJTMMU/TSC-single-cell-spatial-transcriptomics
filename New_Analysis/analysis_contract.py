import json
import os

import numpy as np
import pandas as pd


DATA_CONTRACT_VERSION = "2026-05-14"
DEFAULT_EXPRESSION_LAYER = "log1p"


def stamp_data_contract(adata, producer_script):
    """Persist a machine-readable contract so downstream scripts reuse the same matrices."""
    adata.uns["data_contract"] = {
        "version": DATA_CONTRACT_VERSION,
        "producer_script": producer_script,
        "X": "Scaled HVG matrix for embedding, clustering, and visualization only",
        "layers[counts]": "Raw counts on the saved feature space; never normalized or scaled",
        "layers[log1p]": "Log1p-normalized expression on the saved feature space; default for NMF/GRN/DE summaries",
        "raw": "Full-gene raw counts captured before HVG subsetting; use when complete gene space is required",
        "default_expression_layer": DEFAULT_EXPRESSION_LAYER,
        "nmf_input_layer": DEFAULT_EXPRESSION_LAYER,
        "grn_input_layer": DEFAULT_EXPRESSION_LAYER,
        "spatial_neighbor_policy": "Build neighbors within each sample only",
    }
    return adata


def require_data_contract(adata, caller, require_raw=False):
    contract = adata.uns.get("data_contract")
    if not isinstance(contract, dict):
        raise RuntimeError(
            f"{caller}: missing `adata.uns['data_contract']`. "
            "Please regenerate the object with `01_Data_Preprocessing.py`."
        )

    if "counts" not in adata.layers:
        raise RuntimeError(f"{caller}: missing `adata.layers['counts']` required by the data contract.")
    if DEFAULT_EXPRESSION_LAYER not in adata.layers:
        raise RuntimeError(
            f"{caller}: missing `adata.layers['{DEFAULT_EXPRESSION_LAYER}']` required by the data contract."
        )
    if require_raw and adata.raw is None:
        raise RuntimeError(f"{caller}: full-gene `adata.raw` is required but unavailable.")

    return contract


def _to_numpy(matrix):
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    else:
        matrix = np.asarray(matrix)
    return np.asarray(matrix)


def get_expression_matrix(adata, layer=DEFAULT_EXPRESSION_LAYER, nonnegative=False):
    """Return a contract-compliant expression matrix; never falls back to `.X` silently."""
    require_data_contract(adata, "get_expression_matrix")

    if layer == "X":
        raise RuntimeError(
            "Direct use of `adata.X` is disabled by the data contract. "
            "Use `layers['log1p']` for expression analyses."
        )

    if layer == "raw":
        if adata.raw is None:
            raise RuntimeError("Requested `raw` matrix but `adata.raw` is unavailable.")
        matrix = adata.raw[:, adata.var_names].X
    else:
        if layer not in adata.layers:
            raise RuntimeError(f"Requested layer '{layer}' is missing.")
        matrix = adata.layers[layer]

    matrix = _to_numpy(matrix)
    matrix = np.nan_to_num(matrix, nan=0.0)
    if nonnegative:
        matrix[matrix < 0] = 0
    return matrix


def compute_samplewise_spatial_neighbors(adata, spatial_key="spatial", sample_key="sample", **kwargs):
    """Always build Squidpy spatial neighbors within each sample only."""
    if spatial_key not in adata.obsm:
        raise RuntimeError(f"Missing `adata.obsm['{spatial_key}']`.")

    import squidpy as sq

    base_kwargs = {"coord_type": "generic", "spatial_key": spatial_key}
    base_kwargs.update(kwargs)

    if sample_key in adata.obs and adata.obs[sample_key].nunique() > 1:
        base_kwargs["library_key"] = sample_key

    sq.gr.spatial_neighbors(adata, **base_kwargs)


def query_radius_neighbors_by_sample(adata, target_mask, radius, spatial_key="spatial", sample_key="sample"):
    """Radius-neighbor search that never connects cells across different samples."""
    if spatial_key not in adata.obsm:
        raise RuntimeError(f"Missing `adata.obsm['{spatial_key}']`.")

    from scipy.spatial import cKDTree

    target_mask = np.asarray(target_mask, dtype=bool)
    sample_values = (
        adata.obs[sample_key].astype(str).to_numpy()
        if sample_key in adata.obs
        else np.repeat("all_samples", adata.n_obs)
    )
    coords = np.asarray(adata.obsm[spatial_key], dtype=float)

    all_indices = []
    for sample_id in pd.unique(sample_values):
        sample_indices = np.where(sample_values == sample_id)[0]
        sample_target_indices = sample_indices[target_mask[sample_indices]]
        if sample_target_indices.size == 0:
            continue

        tree = cKDTree(coords[sample_indices])
        local_target_idx = np.searchsorted(sample_indices, sample_target_indices)
        local_neighbors = tree.query_ball_point(coords[sample_target_indices], r=radius)
        if not local_neighbors:
            continue

        local_unique = np.unique(np.concatenate(local_neighbors)).astype(int)
        all_indices.extend(sample_indices[local_unique].tolist())

    return np.array(sorted(set(all_indices)), dtype=int)


def nmf_artifact_paths(output_dir, sample_id):
    prefix = os.path.join(output_dir, sample_id)
    return {
        "cell_weights": f"{prefix}_NMF_CellWeights.csv",
        "gene_weights": f"{prefix}_NMF_GeneWeights.csv",
        "metadata": f"{prefix}_NMF_Metadata.json",
        "top_genes": f"{prefix}_NMF_TopGenes.csv",
    }


def save_nmf_artifacts(output_dir, sample_id, obs_names, gene_names, W, H, source_layer, top_n=10):
    os.makedirs(output_dir, exist_ok=True)
    paths = nmf_artifact_paths(output_dir, sample_id)

    factor_cols = [f"Factor_{i}" for i in range(W.shape[1])]

    cell_weights = pd.DataFrame(W, index=pd.Index(obs_names, name="cell_id"), columns=factor_cols)
    cell_weights.to_csv(paths["cell_weights"])

    gene_weights = pd.DataFrame(H.T, index=pd.Index(gene_names, name="gene"), columns=factor_cols)
    gene_weights.to_csv(paths["gene_weights"])

    top_genes = {}
    for factor in factor_cols:
        top_genes[factor] = gene_weights[factor].sort_values(ascending=False).head(top_n).index.tolist()
    pd.DataFrame({k: pd.Series(v) for k, v in top_genes.items()}).to_csv(paths["top_genes"], index=False)

    with open(paths["metadata"], "w", encoding="utf-8") as handle:
        json.dump(
            {
                "data_contract_version": DATA_CONTRACT_VERSION,
                "source_layer": source_layer,
                "n_components": int(W.shape[1]),
                "n_cells": int(W.shape[0]),
                "n_genes": int(H.shape[1]),
            },
            handle,
            ensure_ascii=True,
            indent=2,
        )

    return paths


def load_nmf_artifacts(output_dir, sample_id):
    paths = nmf_artifact_paths(output_dir, sample_id)
    missing = [path for path in paths.values() if not os.path.exists(path)]
    if missing:
        missing_str = "\n".join(missing)
        raise FileNotFoundError(
            f"Missing NMF artifacts for '{sample_id}'. Please run `07_Spatial_Niche_NMF.py` first:\n{missing_str}"
        )

    cell_weights = pd.read_csv(paths["cell_weights"], index_col=0)
    gene_weights = pd.read_csv(paths["gene_weights"], index_col=0)
    with open(paths["metadata"], "r", encoding="utf-8") as handle:
        metadata = json.load(handle)

    return cell_weights, gene_weights, metadata, paths
