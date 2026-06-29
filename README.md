# HCI Xenium

Spatial transcriptomics analysis of 10x Genomics Xenium data performed by Burke Lawlor (advisor: Dekker Deacon). This project is part of the Judson Torres lab at the University of Utah's Hunstman Cancer Institute (HCI).

The scope of this repository covers multiple analyses which are connected by a general topic of understanding immune-related adverse events (irAE) in melanoma treatmnet.

### Overview of samples
| Disease state                     | Abbreviation | irAE   | n samples |
| ----------                        | ------------ | ----   | --------- | 
| Normal skin                       | NS           | No     | 4 | 
| RMC induced EPD-like condition    | RMC EPD      | Yes    | 7 | 
| Non-treatment related EPD         | NL EPD       | No     | 3 | 
| RMC-induced PNT-like condition    | PNT          | Yes    | 3 |  
| irAE SCLE                         | (SCLE)       | Yes    | 2 | 
| Non-irAE SCLE                     | (SCLE)       | No     | 3 | 
| Non-irAE LP                       | LP           | No     | 3 | 
| non-irAE OLP                      | OLP          | No     | 6 | 
| irAE LPP                          | LPP          | Yes    | 4 | 
| Mastocytosis                      | MAST         | No     | 6 | 
| Transition cases                  | TC           | No     | 1 |   

### Relevant definitions
* **immune-related adverse event (irAE)**
* **Revolutions medicine corporation (RMC)** - makes RMC-6236 aka daraxonrasib, a RAS-targeted chemotherapy drug
* **Papulonecrotic Trunk (PNT)** - rare skin condition that appears as crops of small, dusky-red bumps
* **Subacute cutaneous lupus erythematosus (SCLE)** - a distinct, highly photosensitive autoimmune skin condition
* **Lichen planus (LP)** - a non-contagious, chronic inflammatory condition where your immune system mistakenly attacks cells of your skin, hair, nails, or mucous membranes
* **Oral Lichen planus (OLP)** - LP of the mouth
* **Lichen planus pemphigoides (LPP)** - a rare, acquired autoimmune disease that features overlapping symptoms of both lichen planus and bullous pemphigoid
* **Mastocytosis** - a rare condition caused by the excessive buildup of mast cells
* **Transition cases (TC)** 



## Repository Layout

```
├── figures/            # Output plots, organized by type
├── notebooks/          # Main analysis notebooks 
    ├── rmc-analysis/
    ├── lp-analysis/
    └── masto-analysis/
├── scripts/            # Utility scripts
    └── omeconvert.py/          # converts a TIFF to OME-TIFF (H&E preprocessing)
├── utils/              # Shared Python utilities
    ├── adata_processing.py     # anndata helpers
    ├── adata_loading.py        # sample loading helpers
    ├── diffex.py               # differential expression pipeline
    ├── marker_genes.py         # maker gene maps
    ├── proportions.py          # plotting and stats for proportion analyses
    ├── qc.py                   # single cell qc helpers
    └── spatial.py              # spatial plotting functions
├── output/             # Tabular analysis results - not tracked in git
└── data/               # See below for more - not tracked in git
```

**Notebooks** are organized by analysis project and numbered in order of intended run order.


## Data

The data that powers these analyses is not tracked in git. The contents of the data/ directory can accessed at https://uofu.box.com/s/6dgjml3v3gnl9igxitjjm7fqo5icvcm9. To run the notebooks in this repository, download this folder from Box, rename it to `data`, and store it in your local hci-xenium repo. Note that ~275 gb of storage is required for the entire dataset. 

The data directory is laid out as follows:
```
├── data/
    ├── processed/      # Processed data
        ├── adata/              # Processed AnnData h5ad files
        ├── he/                 # Processed h&e images and alignment files
    ├── raw/            # Original, immutable data
        ├── lymphoid_regions/   # Lymphoid agg coords drawn by Dekker
        └── Xenium/             # Raw Xenium data plus metadata

```


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