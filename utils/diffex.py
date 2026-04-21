import numpy as np
import pandas as pd
import scipy.sparse
import matplotlib.pyplot as plt
from adjustText import adjust_text
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def pseudobulk(adata, sample_col, celltype_col, condition_col, min_cells=20, counts_layer='raw_counts'):
    """
    Create pseudobulk counts by summing raw counts per (sample, cell type, condition) group.

    Groups with fewer than min_cells cells are dropped before aggregation.

    Args:
        adata: AnnData with raw integer counts in layers[counts_layer]
        sample_col: obs column for sample identity
        celltype_col: obs column for cell type labels
        condition_col: obs column for condition labels (must be string/category, not bool)
        min_cells: minimum cells required per pseudobulk group
        counts_layer: layer name containing raw counts

    Returns:
        pb_counts: DataFrame (n_groups × n_genes) of summed integer counts
        pb_meta: DataFrame (n_groups × 3) with sample, celltype, condition columns
    """
    counts = adata.layers[counts_layer]
    if scipy.sparse.issparse(counts):
        counts = counts.toarray()

    genes = adata.var_names.tolist()
    obs = adata.obs[[sample_col, celltype_col, condition_col]].copy()
    obs['__group_key__'] = (
        obs[sample_col].astype(str) + '||' +
        obs[celltype_col].astype(str) + '||' +
        obs[condition_col].astype(str)
    )

    group_sizes = obs['__group_key__'].value_counts()
    keep = group_sizes[group_sizes >= min_cells].index
    mask = obs['__group_key__'].isin(keep).to_numpy()

    obs_f = obs.loc[mask].copy()
    counts_f = counts[mask, :]

    group_keys = obs_f['__group_key__'].unique()
    key_to_row = {k: i for i, k in enumerate(group_keys)}
    row_ids = obs_f['__group_key__'].map(key_to_row).to_numpy()

    pb = np.zeros((len(group_keys), counts_f.shape[1]), dtype=np.float64)
    for i, r in enumerate(row_ids):
        pb[r] += counts_f[i]

    pb_meta = pd.DataFrame(
        [k.split('||') for k in group_keys],
        columns=[sample_col, celltype_col, condition_col],
        index=group_keys,
    ).astype('category')
    pb_counts = pd.DataFrame(pb, index=group_keys, columns=genes).astype(int)

    return pb_counts, pb_meta


def run_pydeseq2_per_celltype(pb_counts, pb_meta, celltype_col, condition_col, contrast,
                               design_formula=None, min_pseudobulk_samples=4):
    """
    Run PyDESeq2 per cell type, comparing two conditions.

    Args:
        pb_counts: DataFrame from pseudobulk()
        pb_meta: DataFrame from pseudobulk()
        celltype_col: column in pb_meta with cell type labels
        condition_col: column in pb_meta with condition labels
        contrast: tuple (condition_col, numerator, denominator),
                  e.g. ('condition', 'normal', 'RMC_EPD')
        design_formula: DESeq2 formula string, e.g. '~ condition' or '~ sample + condition'.
                        Defaults to f'~ {condition_col}'.
        min_pseudobulk_samples: minimum pseudobulk rows per cell type to attempt DE

    Returns:
        dict mapping cell type name -> results DataFrame sorted by padj
    """
    if design_formula is None:
        design_formula = f'~ {condition_col}'

    results = {}
    for ct, meta_ct in pb_meta.groupby(celltype_col):
        idx = meta_ct.index
        counts_ct = pb_counts.loc[idx]
        meta_ct = meta_ct.copy()

        if len(meta_ct) < min_pseudobulk_samples:
            print(f'[SKIP] {ct}: only {len(meta_ct)} pseudobulk samples (need {min_pseudobulk_samples})')
            continue
        if meta_ct[condition_col].nunique() < 2:
            print(f'[SKIP] {ct}: only one condition represented')
            continue
        if (meta_ct.groupby(condition_col, observed=True).size() < 2).any():
            print(f'[SKIP] {ct}: a condition has fewer than 2 samples')
            continue

        try:
            dds = DeseqDataSet(counts=counts_ct, metadata=meta_ct, design=design_formula)
            dds.deseq2()
            stats_res = DeseqStats(dds, contrast=contrast)
            stats_res.summary()
            result = stats_res.results_df.copy()
            result['gene'] = result.index
            result = result.sort_values('padj')
            results[ct] = result
            n_sig = (result['padj'] < 0.05).sum()
            print(f'[OK] {ct}: {n_sig} significant genes (padj < 0.05)')
        except Exception as e:
            print(f'[WARN] {ct}: skipped due to error — {e}')

        print('#' * 80)

    return results


def run_diffex(adata, sample_col, celltype_col, condition_col, contrast,
               design_formula=None, min_cells=20, min_pseudobulk_samples=4,
               counts_layer='raw_counts'):
    """
    Full pseudobulk + PyDESeq2 pipeline per cell type.

    Args:
        adata: AnnData with raw counts in layers[counts_layer]
        sample_col: obs column for sample identity
        celltype_col: obs column for cell type labels
        condition_col: obs column for condition labels (string/category)
        contrast: tuple (condition_col, numerator, denominator)
        design_formula: DESeq2 formula; defaults to f'~ {condition_col}'
        min_cells: min cells per pseudobulk group
        min_pseudobulk_samples: min pseudobulk rows per cell type for DE
        counts_layer: layer with raw integer counts

    Returns:
        dict mapping cell type -> results DataFrame (sorted by padj)
    """
    pb_counts, pb_meta = pseudobulk(
        adata, sample_col, celltype_col, condition_col,
        min_cells=min_cells, counts_layer=counts_layer,
    )
    return run_pydeseq2_per_celltype(
        pb_counts, pb_meta, celltype_col, condition_col, contrast,
        design_formula=design_formula, min_pseudobulk_samples=min_pseudobulk_samples,
    )


def plot_volcano(results_df, title='', padj_thresh=0.05, lfc_thresh=1.0, top_n_labels=10,
                 figsize=(7, 6)):
    """
    Volcano plot for a single cell type's DESeq2 results.

    Args:
        results_df: results DataFrame from run_diffex() (one cell type)
        title: plot title (e.g. cell type name)
        padj_thresh: adjusted p-value threshold for significance
        lfc_thresh: absolute log2 fold-change threshold for up/down coloring
        top_n_labels: number of top significant genes (by padj) to label
        figsize: figure size tuple

    Returns:
        matplotlib Figure
    """
    df = results_df.dropna(subset=['padj', 'log2FoldChange', 'pvalue']).copy()
    df['neg_log10_p'] = -np.log10(df['pvalue'].clip(lower=1e-300))

    sig_up = (df['padj'] < padj_thresh) & (df['log2FoldChange'] > lfc_thresh)
    sig_dn = (df['padj'] < padj_thresh) & (df['log2FoldChange'] < -lfc_thresh)
    ns = ~(sig_up | sig_dn)

    colors = {
        'up': '#d62728',
        'down': '#1f77b4',
        'ns': '#aaaaaa',
    }

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(df.loc[ns, 'log2FoldChange'], df.loc[ns, 'neg_log10_p'],
               color=colors['ns'], s=8, alpha=0.5, linewidths=0, label='NS')
    ax.scatter(df.loc[sig_dn, 'log2FoldChange'], df.loc[sig_dn, 'neg_log10_p'],
               color=colors['down'], s=12, alpha=0.8, linewidths=0,
               label=f'Down ({sig_dn.sum()})')
    ax.scatter(df.loc[sig_up, 'log2FoldChange'], df.loc[sig_up, 'neg_log10_p'],
               color=colors['up'], s=12, alpha=0.8, linewidths=0,
               label=f'Up ({sig_up.sum()})')

    # Threshold lines
    ax.axhline(-np.log10(padj_thresh), color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(lfc_thresh, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.axvline(-lfc_thresh, color='black', linestyle='--', linewidth=0.8, alpha=0.6)

    # Label top significant genes
    sig_genes = df[sig_up | sig_dn].nsmallest(top_n_labels, 'padj')
    texts = [
        ax.text(row['log2FoldChange'], row['neg_log10_p'], row['gene'],
                fontsize=7, ha='center', va='bottom')
        for _, row in sig_genes.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    ax.set_xlabel('log$_2$ Fold Change', fontsize=12)
    ax.set_ylabel('-log$_{10}$(p-value)', fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.legend(fontsize=9, framealpha=0.7)
    fig.tight_layout()

    return fig
