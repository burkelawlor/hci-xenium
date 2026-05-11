import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation


def run_qc_single_sample(adata, sample_name, save_fig=False):

    sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)
    adata.obs['nucleus_ratio'] = adata.obs["nucleus_area"] / adata.obs["cell_area"]

    cprobes = (
        adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
    )
    cwords = (
        adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
    )
    print("Metrics for sample: ", sample_name)
    print(f"Negative DNA probe count % : {cprobes}")
    print(f"Negative decoding count % : {cwords}")

    
    fig, axs = plt.subplots(1, 4, figsize=(15, 4))

    axs[0].set_title("Total transcripts per cell")
    sns.histplot(
        adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )

    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )


    axs[2].set_title("Area of segmented cells")
    sns.histplot(
        adata.obs["cell_area"],
        kde=False,
        ax=axs[2],
    )

    axs[3].set_title("Nucleus ratio")
    sns.histplot(
        adata.obs["nucleus_ratio"],
        kde=False,
        ax=axs[3],
    )

    if save_fig:
        plot_title = "./figures/qc_plots/qc_metrics_" + sample_name + ".png"
        plt.savefig(plot_title)
    
    plt.show()



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