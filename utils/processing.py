import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad


def propogate_subset_labels(parent_adata, subset_adata, key_added, key_to_add, merge_on='cell_id_unique'):
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



def get_ranked_genes_by_group(adata, key='rank_genes_groups'):
    result = adata.uns[key]
    groups = result['names'].dtype.names
    genes_ranked = {group:[x[i] for x in result['names']] for i,group in enumerate(groups)}
    return genes_ranked


from scipy.stats import median_abs_deviation
def flag_outliers_by_mad(adata, column, upper_mad=5, lower_mad=5):
    values = adata.obs[column].values
    med = np.median(values)
    mad_val = median_abs_deviation(values)
    
    if upper_mad != 0:
        high_cut = med + upper_mad * mad_val
    else: high_cut = max(values)

    if lower_mad != 0:
        low_cut  = med - lower_mad * mad_val
    else: low_cut = min(values)
    
    column_name = column + "_outlier"
    adata.obs[column_name] = (values > high_cut) | (values < low_cut)
    return adata