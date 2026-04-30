# HCI irAE

Spatial transcriptomics analysis of immune-related adverse events (irAE) using 10x Genomics Xenium data. The project is notebook-driven, with shared utilities in `utils/`.

The project is divided into two disease analyses:

- **RMC analysis** (Reactive Mycosis-like Conditions): compares EPD (erythrodermic psoriasiform dermatitis), PNT (post-natal telogen), and NS (normal skin) samples. Actively in progress.
- **SCLE analysis** (Subacute Cutaneous Lupus Erythematosus): compares SCLE and NS samples. Not yet started.

Analyses primarily use [squidpy](https://squidpy.readthedocs.io), [spatialdata](https://spatialdata.scverse.org), [anndata](https://anndata.readthedocs.io), and [PyDESeq2](https://pydeseq2.readthedocs.io).

---

## Repository Layout

```
notebooks/          # Main analysis notebooks (run sequentially)
utils/              # Shared Python utilities
figures/            # Output plots, organized by type
output/             # Tabular analysis results (e.g. differential expression CSVs)
scripts/            # Utility scripts
tutorials/          # Reference notebooks for learning squidpy/spatialdata/scanpy
data/               # Local data — not tracked in git (see Data section below)
```

### `notebooks/`

Notebooks are numbered and intended to run in order. Notebook `01` preprocesses all samples together; `02*` notebooks are specific to the RMC analysis.

| Notebook | Purpose |
|---|---|
| `01_preprocessing_qc.ipynb` | Load Xenium data for all samples, quality control, normalization, and initial clustering |
| `02a_RMC_annotation_general.ipynb` | Broad cell-type annotation of RMC samples |
| `02b_RMC_annotation_immune.ipynb` | Immune cell subtype annotation |
| `02c_RMC_annotation_endothelial.ipynb` | Endothelial cell annotation |
| `02d_RMC_annotation_epithelial.ipynb` | Epithelial (keratinocyte) annotation |
| `02e_RMC_annotation_detailed.ipynb` | Review and clean up detailed annotations across all cell-type subsets |
| `02f_diffex.ipynb` | Differential expression analysis using PyDESeq2 |
| `02g_lymphoid_aggregates.ipynb` | Spatial analysis of manually annotated lymphoid aggregate regions |
| `02h_h&e.ipynb` | Explore H&E image cuts with spatial overlays |

`notebooks/archived/`, `notebooks/scratch/`, and `notebooks/examples/` hold exploratory or deprecated notebooks and are not part of the main workflow.

### `utils/`

| File | Contents |
|---|---|
| `data_loading.py` | `xenium_paths` — dict mapping sample IDs to raw Xenium output directories |
| `cell_types.py` | `CELL_TYPES_DETAILED` — hierarchical dict of cell types used for annotation |
| `plotting.py` | Spatial plots, UMAP feature plots, proportion bar/line plots |
| `processing.py` | Subset label propagation, ranked gene extraction, MAD-based outlier flagging |
| `diffex.py` | Differential expression helpers (PyDESeq2 wrappers) |

### `figures/`

Output plots organized into subdirectories by plot type: `spatial_plots/`, `umaps/`, `feature_plots/`, `dotplots/`, `proportions/`.

### `output/`

Tabular results from differential expression and other analyses, saved as CSV files.

### `scripts/`

- `omeconvert.py` — converts a TIFF to OME-TIFF with pyramid levels for viewing in Xenium Explorer

### `tutorials/`

Reference notebooks covering scanpy basics, spatialdata, and squidpy integration. Not part of the main analysis workflow.

---

## Data

Raw and processed data live under `data/` and are not tracked in git. The data directory contains raw Xenium outputs (two batches), processed AnnData `.h5ad` files, H&E images, manually annotated lymphoid region files, and a sample metadata spreadsheet.

---

## Environment

The project uses a local Python 3.12 virtual environment at `.env/`:

```bash
source .env/bin/activate
jupyter lab
```
