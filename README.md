# HCI irAE

Spatial transcriptomics analyses of immune-related adverse events (irAE) using Xenium data.

## About

Analyses can be found in `notebooks/` which are laid out in the following order:

* **01_preprocessing_qc.ipynb** — Load Xenium data, preprocessing, QC  
* **02x_annotation_X.ipynb** — Cell-type annotations for RMC samples (PNT1-2, EPD1-7) according to the notebook name (general, immune, endothelial, keratinocyte).  
* **03_lymphoid_aggregates.ipynb** — Analysis of the lymphoid aggreate regions, which have been manually drawn by Dekker

Shared utilities live in `utils/` (plotting, processing, cell types). Plots are written to `figures/`. Input data is expected under `data/` (not in repo).
