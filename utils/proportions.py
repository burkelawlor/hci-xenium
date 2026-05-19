import itertools

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from scipy import stats
from statsmodels.stats.multitest import multipletests


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


def _per_sample_proportions(adata, ct_col, cell_type, groupby, sample_name_col):
    obs = adata.obs[[ct_col, groupby, sample_name_col]].copy()
    per_sample = (
        obs.groupby(sample_name_col, observed=True)
        .apply(lambda df: (df[ct_col] == cell_type).sum() / len(df), include_groups=False)
        .rename('proportion')
    )
    sample_to_group = obs.drop_duplicates(sample_name_col).set_index(sample_name_col)[groupby]
    return pd.DataFrame({'proportion': per_sample, 'group': sample_to_group})


def test_proportion_differences(
    adata,
    ct_col,
    cell_type,
    groupby,
    sample_name_col='sample_name',
    order=None,
    correction='fdr_bh',
):
    per_sample_df = _per_sample_proportions(adata, ct_col, cell_type, groupby, sample_name_col)

    groups = order if order is not None else list(per_sample_df['group'].unique())
    groups = [g for g in groups if g in per_sample_df['group'].values]
    group_vals = [per_sample_df[per_sample_df['group'] == g]['proportion'].values for g in groups]

    kruskal_stat, kruskal_p = stats.kruskal(*group_vals)

    pairs = list(itertools.combinations(groups, 2))
    raw_pvals = [
        stats.mannwhitneyu(
            per_sample_df[per_sample_df['group'] == a]['proportion'].values,
            per_sample_df[per_sample_df['group'] == b]['proportion'].values,
            alternative='two-sided',
        ).pvalue
        for a, b in pairs
    ]
    _, adj_pvals, _, _ = multipletests(raw_pvals, method=correction)

    pairwise = pd.DataFrame(np.nan, index=groups, columns=groups)
    for (a, b), p in zip(pairs, adj_pvals):
        pairwise.loc[a, b] = p
        pairwise.loc[b, a] = p

    return {
        'kruskal_stat': kruskal_stat,
        'kruskal_p': kruskal_p,
        'pairwise': pairwise,
        '_pairs': pairs,
        '_adj_pvals': list(adj_pvals),
    }


def _draw_significance_brackets(ax, groups, stats_result, means, errs, per_sample_df):
    def _stars(p):
        if p < 0.001: return '***'
        if p < 0.01: return '**'
        if p < 0.05: return '*'
        return None

    pairs = stats_result['_pairs']
    adj_pvals = stats_result['_adj_pvals']

    sig_pairs = [(pair, p) for pair, p in zip(pairs, adj_pvals) if _stars(p) is not None]
    if not sig_pairs:
        return

    # Sort by span so shorter (inner) brackets are drawn first
    sig_pairs.sort(key=lambda x: abs(groups.index(x[0][1]) - groups.index(x[0][0])))

    # Start all brackets above the tallest element anywhere in the plot
    y_max = max(
        max(means[g] + errs[g] for g in groups),
        per_sample_df['proportion'].max(),
    )
    heights = {i: y_max for i in range(len(groups))}
    step = y_max * 0.12

    n_sig = len(sig_pairs)
    ax.set_ylim(top=y_max + step * (n_sig + 2))

    for (g1, g2), p in sig_pairs:
        i1, i2 = groups.index(g1), groups.index(g2)
        y = max(heights[i1], heights[i2]) + step
        ax.plot([i1, i1, i2, i2], [heights[i1] + step * 0.3, y, y, heights[i2] + step * 0.3],
                lw=1, c='black')
        ax.text((i1 + i2) / 2, y + step * 0.05, _stars(p), ha='center', va='bottom', fontsize=10)
        for i in range(min(i1, i2), max(i1, i2) + 1):
            heights[i] = y


def plot_proportion_scatter_bar(
    adata,
    ct_col,
    cell_type,
    groupby,
    sample_name_col='sample_name',
    order=None,
    error='sem',
    color=None,
    stats=None,
    figsize=None,
    title='default',
    save=False,
    ax=None,
):
    per_sample_df = _per_sample_proportions(adata, ct_col, cell_type, groupby, sample_name_col)

    groups = order if order is not None else list(per_sample_df['group'].unique())
    groups = [g for g in groups if g in per_sample_df['group'].values]

    means = per_sample_df.groupby('group', observed=True)['proportion'].mean()
    if error == 'std':
        errs = per_sample_df.groupby('group', observed=True)['proportion'].std()
    else:
        errs = per_sample_df.groupby('group', observed=True)['proportion'].sem()

    if color is None:
        try:
            categories = list(adata.obs[ct_col].cat.categories)
            colors = adata.uns[f'{ct_col}_colors']
            color = colors[categories.index(cell_type)]
        except Exception:
            color = 'steelblue'

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    x_pos = range(len(groups))

    ax.bar(x_pos, [means[g] for g in groups], yerr=[errs[g] for g in groups],
           color=color, alpha=0.7, capsize=4, error_kw={'linewidth': 1.5})

    rng = np.random.default_rng(0)
    for i, g in enumerate(groups):
        vals = per_sample_df[per_sample_df['group'] == g]['proportion'].values
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(i + jitter, vals, color='black', s=25, zorder=3, linewidths=0)

    ax.set_xticks(list(x_pos))
    ax.set_xticklabels(groups)
    ax.set_ylabel('Proportion')
    ax.set_xlabel(groupby)

    if title == 'default':
        ax.set_title(f'Proportion of {cell_type} by {groupby}')
    else:
        ax.set_title(title)

    if stats is not None:
        _draw_significance_brackets(ax, groups, stats, means, errs, per_sample_df)

    if ax is None:
        plt.tight_layout()

        if save is True:
            plt.savefig(f'./figures/proportions/{cell_type}_proportion_by_{groupby}.png', bbox_inches='tight')
        elif isinstance(save, str):
            plt.savefig(f'./figures/proportions/{save}.png', bbox_inches='tight')
        plt.show()
    
    return ax


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