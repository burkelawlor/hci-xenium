# HCI Xenium

Spatial transcriptomics analysis of 10x Genomics Xenium data performed by Burke Lawlor (advisor: Dekker Deacon). This project is part of the Judson Torres lab at the University of Utah's Hunstman Cancer Institute (HCI).

The scope of this repository covers multiple analyses which are connected by a general topic of understanding immune-related adverse events (irAE) in melanoma treatmnet.

Samples groups include:

| Disease state                     | irAE Condition?   | n samples |
| ----------                        | -----------       | --------- | 
| Normal skin (NS)                  | No                | 4 | 
| Reactive Mycosis-like Condition with erosive pustular dermatosis (RMC EPD)    | Yes      | 7 | 
| Non-treatment-related adverse event erosive pustular dermatosis (NL EPD)      | No       | 3 | 
| Reactive Mycosis-like Condition with Papulonecrotic Trunk (PNT)               | Yes      | 3 |  
| Non irAE Subacute cutaneous lupus erythematosus (SCLE)                        | No       | 5 | 
| non irAE Lichen planus (LP)               | No                | 3 | 
| irAE Lichen planus pemphigoides (LPP)     | Yes               | 4 | 
| Oral lichen planus (OLP)                  | ?                 | 6 | 
| Mastocytosis (MAST)                       | ?                 | 6 | 
| ? TC                                      | TC                | 1 |   


## Repository Layout

```
├── notebooks/          # Main analysis notebooks 
    ├── rmc-analysis/
    ├── lp-analysis/
    ├── masto-analysis/
├── utils/              # Shared Python utilities
    ├── adata_processing.py     # anndata helpers
    ├── adata_loading.py        # sample loading helpers
    ├── diffex.py               # differential expression pipeline
    ├── marker_genes.py         # maker gene maps
    ├── proportions.py          # plotting and stats for proportion analyses
    ├── qc.py                   # single cell qc helpers
    ├── spatial.py              # spatial plotting functions
├── figures/            # Output plots, organized by type
├── output/             # Tabular analysis results
├── scripts/            # Utility scripts
    ├── omeconvert.py/          # converts a TIFF to OME-TIFF (H&E preprocessing)
├── data/               # Local data — not tracked in git
    ├── processed/      # Processed data
        ├── adata/              # Processed AnnData h5ad files
        ├── he/                 # Processed h&e images and alignment files
    ├── raw/            # Original, immutable data
        ├── lymphoid_regions/   # Lymphoid agg coords drawn by Dekker
        ├── Xenium/             # Raw Xenium data plus metadata
```

**Notebooks** are organized by analysis project and numbered in order of intended run order. See README

---

## Environment

The project was developed in Python 3.12. A full pinned snapshot of the environment is in `requirements.txt`. Key top-level packages:

| Package | Version |
|---|---|
| `scanpy` | 1.12 |
| `squidpy` | 1.8.1 |
| `anndata` | 0.12.10 |
| `spatialdata` | 0.7.2 |
| `spatialdata-io` | 0.6.0 |
| `spatialdata-plot` | 0.2.14 |
| `pydeseq2` | 0.5.4 |
| `gseapy` | 1.2.1 |
| `sopa` | 2.2.5 |
| `pandas` | 2.3.3 |
| `numpy` | 2.4.2 |
| `scipy` | 1.16.3 |
| `scikit-learn` | 1.8.0 |
| `statsmodels` | 0.14.6 |
| `matplotlib` | 3.10.8 |
| `seaborn` | 0.13.2 |
| `glasbey` | 0.3.0 |
| `adjustText` | 1.3.0 |
| `geopandas` | 1.1.2 |