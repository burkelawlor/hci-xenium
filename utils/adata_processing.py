import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

from utils.qc import flag_outliers_by_mad

def propogate_subset_labels(dest_adata, source_adata, dest_key, source_key, merge_on='cell_id_unique'):
    """
    Propogate labels from source_adata.obs[source_key] to dest_adata.obs[dest_key].

    Args:
        dest_adata: AnnData object to add the labels to
        source_adata: AnnData object to get the labels from
        dest_key: Name of the new column to add to dest_adata.obs
        source_key: Name of the column in source_adata.obs to add to dest_adata.obs
        merge_on: Name of the column in source_adata.obs and dest_adata.obs to merge on
    """

    adata_new = dest_adata.copy()

    if dest_key not in adata_new.obs.columns:
        adata_new.obs[dest_key] = np.nan

    # Build a lookup from merge key -> label, avoiding merge suffix collisions
    # when dest_key and source_key are the same column name.
    subset_lookup = (
        source_adata.obs[[merge_on, source_key]]
        .dropna(subset=[source_key])
        .drop_duplicates(subset=[merge_on], keep="last")
        .set_index(merge_on)[source_key]
    )

    propagated = adata_new.obs[merge_on].map(subset_lookup)
    adata_new.obs[dest_key] = adata_new.obs[dest_key].astype("object")
    mask = propagated.notna()
    adata_new.obs.loc[mask, dest_key] = propagated.loc[mask].values
    adata_new.obs[dest_key] = adata_new.obs[dest_key].astype("category")

    return adata_new



def get_ranked_genes_by_group(adata, key='rank_genes_groups'):
    result = adata.uns[key]
    groups = result['names'].dtype.names
    genes_ranked = {group:[x[i] for x in result['names']] for i,group in enumerate(groups)}
    return genes_ranked


def feature_plots_from_marker_genes(adata, marker_genes_dict, save=False, prefix=None):
    for ct in marker_genes_dict:
        print(f'{ct.upper()}:')
        
        ax = sc.pl.umap(adata, color=marker_genes_dict[ct], wspace=0.1, show=False)
        plt.suptitle(ct)
        
        if save:
            if isinstance(ax, list):
                ax[0].figure.savefig(f"./figures/feature_plots/{prefix}_umap_{ct}.png", bbox_inches="tight")
            else:
                ax.figure.savefig(f"./figures/feature_plots/{prefix}_umap_{ct}.png", bbox_inches="tight")

        plt.show()