import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad


def propogate_subset_labels(parent_adata, subset_adata, key_added, key_to_add, merge_on):
    """
    Propogate labels from subset_adata.obs[key_to_add] to parent_adata.obs[key_added].

    Args:
        parent_adata: AnnData object to add the labels to
        subset_adata: AnnData object to get the labels from
        key_added: Name of the new column to add to parent_adata.obs
        key_to_add: Name of the column in subset_adata.obs to add to parent_adata.obs
        merge_on: Name of the column in subset_adata.obs and parent_adata.obs to merge on
    """

    adata_new = parent_adata.copy()

    if key_added not in adata_new.obs.columns:
        adata_new.obs[key_added] = np.nan

    # Build a lookup from merge key -> label, avoiding merge suffix collisions
    # when key_added and key_to_add are the same column name.
    subset_lookup = (
        subset_adata.obs[[merge_on, key_to_add]]
        .dropna(subset=[key_to_add])
        .drop_duplicates(subset=[merge_on], keep="last")
        .set_index(merge_on)[key_to_add]
    )

    propagated = adata_new.obs[merge_on].map(subset_lookup)
    adata_new.obs[key_added] = adata_new.obs[key_added].astype("object")
    mask = propagated.notna()
    adata_new.obs.loc[mask, key_added] = propagated.loc[mask].values
    adata_new.obs[key_added] = adata_new.obs[key_added].astype("category")

    return adata_new