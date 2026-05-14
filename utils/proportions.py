import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory


def plot_proportions_bar(
    adata, 
    groupby, 
    ct_col, 
    palette=None, 
    save=False, 
    figsize=None, 
    title='default', 
    order=None,
    normalize=True,
    display_totals=False,
    display_totals_title='Total count',
    legend_bbox=None,
    return_proportions=False,
    ylabel=None,
):
    counts = adata.obs[[groupby, ct_col]].groupby(groupby, observed=True).value_counts(normalize=False).unstack()
    proportions = adata.obs[[groupby, ct_col]].groupby(groupby, observed=True).value_counts(normalize=normalize).unstack()
    

    # Reorder the index according to 'order' if provided
    if order is not None:
        # Only keep values present in proportions' index
        order = [x for x in order if x in proportions.index]
        proportions = proportions.loc[order]
    else: 
        order = proportions.index

    if palette is None:
        try:
            plot_ax = proportions.plot.barh(stacked=True, figsize=figsize, color=adata.uns[f'{ct_col}_colors'])
        except Exception:
            plot_ax = proportions.plot.barh(stacked=True, figsize=figsize)
    else:
        plot_ax = proportions.plot.barh(stacked=True, figsize=figsize, color=palette)
    
    # In matplotlib, the first index appears at the bottom of the y-axis.
    # To put the first value of 'order' at the TOP, invert the y-axis.
    plot_ax.invert_yaxis()
    
    
    if not legend_bbox:
        plt.legend(bbox_to_anchor=(1.05,0.5), loc='center left')
    else:
        plt.legend(bbox_to_anchor=legend_bbox, loc='center left')
    
    if title == 'default':
        plt.title(f"{ct_col} proportions by {groupby}")
    else:
        plt.title(title)

    if normalize:
        plt.xlabel('Proportion')
    else:
        plt.xlabel('Count')

    if isinstance(ylabel, str):
        plt.ylabel(ylabel)
    elif ylabel is None:
        plt.ylabel(None)
    
    if display_totals is not False:
        if isinstance(display_totals, dict):
            total_values = [str(display_totals[group]) for group in order]
        elif isinstance(display_totals, list):
            total_values = [str(v) for v in display_totals]
        else:
            computed = counts.sum(axis=1)[order]
            fmt = '{:,.0f}' if computed.dtype == 'int64' else '{:.3f}'
            total_values = [fmt.format(v) for v in computed]

        trans = blended_transform_factory(plot_ax.transAxes, plot_ax.transData)
        val_texts = []
        for i, val in enumerate(total_values):
            val_texts.append(plot_ax.text(
                1.02, i, val,
                transform=trans, va='center', ha='left',
                fontsize=plt.rcParams.get('ytick.labelsize', 'medium'),
            ))

        # Position title just past the rightmost edge of the value labels
        plot_ax.figure.canvas.draw()
        renderer = plot_ax.figure.canvas.get_renderer()
        max_x_disp = max(t.get_window_extent(renderer).x1 for t in val_texts)
        max_x_ax = plot_ax.transAxes.inverted().transform((max_x_disp, 0))[0]
        plot_ax.text(max_x_ax + 0.04, 0.5, display_totals_title, transform=plot_ax.transAxes, va='center', ha='center', rotation=-90)
        
        if not legend_bbox:
            plt.legend(bbox_to_anchor=(max_x_ax + 0.1, 0.5), loc='center left')
        else:
            plt.legend(bbox_to_anchor=legend_bbox, loc='center left')

    if save is True:
        plt.savefig(f"./figures/proportions/{ct_col}_by_{groupby}_bar.png", bbox_inches="tight")
    elif isinstance(save, str):
        plt.savefig(f"./figures/proportions/{save}.png", bbox_inches="tight")
    plt.show()

    if return_proportions:
        return proportions
        


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