import math
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import matplotlib.pyplot as plt
import numpy as np
import squidpy as sq
import scanpy as sc
import anndata as ad

def spatial_plot_cell_types_layered(
    adata,
    ct_col,
    subset=False,
    save=False,
    figsize=None,
    n_cols=None,
    size=1,
    dpi='figure'
):
    sample_names = sorted(adata.obs["sample_name"].unique())
    n = len(sample_names)
    if n == 0:
        return
    if n_cols is None:
        ncols = min(3, n)
    else:
        ncols = max(1, min(int(n_cols), n))
    nrows = math.ceil(n / ncols)
    if figsize is None:
        figsize = (4 * ncols, 4 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.ravel(np.atleast_1d(axes))

    legend_cats = None
    legend_colors = None

    for i, sample_name in enumerate(sample_names):
        adata_sample = adata[adata.obs["sample_name"] == sample_name].copy()

        if subset:
            cats = adata_sample.obs[ct_col].cat.categories
            colors = adata_sample.uns.pop(f'{ct_col}_colors')
            colors_mask = cats.isin(subset)

            adata_sample.obs[ct_col] = np.select([adata_sample.obs[ct_col].isin(subset)], [adata_sample.obs[ct_col]], np.nan)
            adata_sample.uns[f'{ct_col}_colors'] = colors[colors_mask]

        if legend_cats is None:
            legend_cats = list(adata_sample.obs[ct_col].cat.categories)
            legend_colors = np.asarray(adata_sample.uns[f'{ct_col}_colors'])

        # Sort so that NA cells are on top and plotted first (to be on the bottom)
        adata_sorted = ad.concat([adata_sample[adata_sample.obs[ct_col].isna()], adata_sample[adata_sample.obs[ct_col].notna()]], uns_merge='first')
        adata_sorted.uns[f'{ct_col}_colors'] = adata_sample.uns[f'{ct_col}_colors']
        
        ax = sq.pl.spatial_scatter(
            adata_sorted,
            library_id="spatial",
            shape=None,
            color=[ct_col],
            size=size,
            na_color='lightgray',
            return_ax=True,
            ax=axes[i],
        )

        leg = ax.get_legend()
        if leg is not None:
            leg.remove()
        ax.set_axis_off()
        ax.set_title(sample_name)

    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    if legend_cats:
        handles = [
            Patch(facecolor=c, edgecolor="none", label=str(cat))
            for cat, c in zip(legend_cats, legend_colors)
        ]
        fig.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            title=ct_col,
            frameon=True,
        )
        fig.tight_layout(rect=[0, 0, 0.88, 1])
    else:
        fig.tight_layout()

    plt.show()

    if save:
        fig.savefig(f'figures/spatial_plots/{ct_col}_layered.png', bbox_inches="tight", dpi=dpi)




def spatial_plot_cell_types_individual(
    adata,
    sample_name,
    ct_col,
    ncols=4,
    figsize=None,
    palette=None,
    size=10,
    save=False,
    include_na=False,
):
    adata_sample = adata[adata.obs["sample_name"] == sample_name].copy()

    # Get cell types 
    ct_series = adata_sample.obs[ct_col]
    cts = list(np.sort(ct_series.dropna().unique().astype(str)))
    
    # Include NA in cell types if requested and present
    has_na = ct_series.isna().any()
    if include_na and has_na:
        cts.append('NA')

    # Plot params
    n_ct = len(cts)
    nrows = math.ceil(n_ct / ncols)
    if figsize is None:
        figsize = (4 * ncols, 4 * nrows)
    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, tight_layout=True)
    ax = np.array(ax).reshape(-1)  # flatten safely for 1 row/col cases

    try:
        palette = dict(zip(adata_sample.obs[ct_col].cat.categories, adata_sample.uns[f'{ct_col}_colors'])) 
        palette['NA'] = '#000000'
    except:
        palette = None

    # Plot each cell type
    for i, ct in enumerate(cts):
        # Build a per-panel mask (special handling for NA)
        if ct == 'NA':
            mask = ct_series.isna()
            panel_label = 'NA'
        else:
            mask = ct_series.astype(str) == ct
            panel_label = ct

        # Make an object/string series to avoid dtype promotion issues
        is_ct = np.where(mask.to_numpy(), panel_label, "Other").astype(object)
        adata_sample.obs["is_ct"] = is_ct

        # Background: all cells (no color)
        sq.pl.spatial_scatter(
            adata_sample,
            library_id="spatial",
            shape=None,
            ax=ax[i],
            na_color="lightgray",
            size=size,
        )

        # Foreground: highlight the current ct (colored by is_ct)
        adata_fg = adata_sample[mask].copy()
        if adata_fg.n_obs > 0:
            if palette is not None and panel_label in palette:
                sq.pl.spatial_scatter(
                    adata_fg,
                    library_id="spatial",
                    shape=None,
                    color="is_ct",
                    ax=ax[i],
                    size=size,
                    palette=ListedColormap([palette[panel_label]]),
                )
            else:
                sq.pl.spatial_scatter(adata_fg,
                    library_id="spatial",
                    shape=None,
                    color="is_ct",
                    ax=ax[i],
                    size=size,
                )

        ax[i].set_title(f"{panel_label}")
        ax[i].invert_yaxis()

        leg = ax[i].get_legend()
        if leg is not None:
            leg.remove()
        ax[i].set_axis_off()

    # Delete unused axes
    for j in range(n_ct, len(ax)):
        fig.delaxes(ax[j])

    plt.suptitle(f"{sample_name}")
    plt.show()

    if save:
        fig.savefig(f"figures/spatial_plots/{sample_name}_{ct_col}_individual.png", bbox_inches="tight")





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


def plot_proportions_bar(
    adata, 
    groupby, 
    ct_col, 
    palette=None, 
    save=False, 
    figsize=None, 
    title='default', 
    order=None,
):
    proportions = adata.obs[[groupby, ct_col]].groupby(groupby, observed=True).value_counts(normalize=True).unstack()

    # Reorder the index according to 'order' if provided
    if order is not None:
        # Only keep values present in proportions' index
        order = [x for x in order if x in proportions.index]
        proportions = proportions.loc[order]

    if palette is None:
        try:
            plot_ax = proportions.plot.barh(
                stacked=True, figsize=figsize, color=adata.uns[f'{ct_col}_colors'],
            )
        except Exception:
            plot_ax = proportions.plot.barh(stacked=True, figsize=figsize)
    else:
        plot_ax = proportions.plot.barh(stacked=True, figsize=figsize, color=palette)
    
    # In matplotlib, the first index appears at the bottom of the y-axis.
    # To put the first value of 'order' at the TOP, invert the y-axis.
    plot_ax.invert_yaxis()
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    if title == 'default':
        plt.title(f"{ct_col} proportions by {groupby}")
    else:
        plt.title(title)
    
    if save is True:
        plt.savefig(f"./figures/proportions/{ct_col}_by_{groupby}_bar.png", bbox_inches="tight")
    elif isinstance(save, str):
        plt.savefig(f"./figures/proportions/{save}.png", bbox_inches="tight")
    plt.show()


def plot_proportions_line(adata, groupby, ct_col, order=None, palette=None, save=False, figsize=None, title='default'):
    proportions = adata.obs[[groupby,ct_col]].groupby(groupby, observed=True).value_counts(normalize=True).unstack()

    # Ensure the x-tick order
    if order is None:
        x_ticks = list(proportions.index)
    else:   
        x_ticks = order
    x_pos = range(len(x_ticks))

    fig, ax = plt.subplots(figsize=figsize)

    for cell_type in proportions.columns:
        y = [proportions.loc[idx, cell_type] if idx in proportions.index else np.nan for idx in x_ticks]
        
        if palette is None:
            try:
                groups = adata.obs[ct_col].cat.categories
                colors = adata.uns[f'{ct_col}_colors']
                palette = dict(zip(groups, colors))
                ax.plot(x_pos, y, marker='o', label=cell_type, color=palette[cell_type])
            except:
                ax.plot(x_pos, y, marker='o', label=cell_type)
        else:
            ax.plot(x_pos, y, marker='o', label=cell_type, color=palette[cell_type])

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_ticks)
    ax.set_ylabel('Proportion')

    if title == 'default':
        plt.title(f"{ct_col} proportions by {groupby}")
    else:
        plt.title(title)

    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    if save == True:
        plt.savefig(f"./figures/proportions/{ct_col}_by_{groupby}_line.png", bbox_inches="tight")
    elif isinstance(save, str):
        plt.savefig(f"./figures/proportions/{save}.png", bbox_inches="tight")
    plt.show()