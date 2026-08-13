from pathlib import Path
import pandas as pd
import numpy as np
import re

import spatialdata as sd
from spatialdata_io import xenium
from spatialdata_io.readers.xenium import xenium_aligned_image
from spatialdata.models import TableModel, get_table_keys

DATA_ROOT = Path("/Users/burkelawlor/Library/CloudStorage/Box-Box/hci-xenium/data")

### Metadata and paths
metadata = pd.read_csv(DATA_ROOT / "metadata.csv", dtype=str).dropna(subset=["File path"])


### Xenium paths
data_root = DATA_ROOT / "raw/Xenium/"
xenium_paths = {row['Study ID']: data_root / row['File path'] for index, row in metadata.iterrows()}

### H&E paths
_he_meta = metadata.dropna(subset=["H&E Image Path"])
he_image_paths = _he_meta.set_index("Study ID")["H&E Image Path"].map(Path).to_dict()
he_alignment_paths = _he_meta.set_index("Study ID")["H&E Alignment Path"].map(Path).to_dict()


### Slide & batch information
studyid_to_slide = metadata.set_index('Study ID')['Slide'].to_dict()

batch_slides = {
    'batch_1': ['0076657', '0076678'],
    'batch_2': ['0084912', '0084920'],
    'batch_3': ['0076414', '0076660'],
    'australia': ['0029954', '0029955'],
    'ici_lp': ['0080734', '0080735'],
}
slide_to_batch = {slide: batch for batch, slides in batch_slides.items() for slide in slides}
studyid_to_batch = metadata.set_index('Study ID')['Slide'].map(slide_to_batch).to_dict()


### Disease state information
sample_sets = {
    'NL': ['NS1', 'NS2', 'NS3', 'NS4'],
    'RMC_EPD': ['EPD1', 'EPD2', 'EPD3', 'EPD4', 'EPD5', 'EPD6', 'EPD7'],
    'NL_EPD': ['EPD8', 'EPD9', 'EPD10'],
    'PNT': ['PNT1', 'PNT2'],
    'LP': ['LP1_LP2','LP3','LP4','LP5','LP6','LP7','LP8'],
    'LPP': ['LPP1', 'LPP2', 'LPP3', 'LPP4'],
    'ICI_LP': ['ICI_LP1', 'ICI_LP2', 'ICI_LP3', 'ICI_LP4'],
    'OLP': ['OLP1', 'OLP2', 'OLP3', 'OLP4', 'OLP5', 'OLP6'],
    'MAST': ['MAST1', 'MAST2', 'MAST3','MAST4', 'MAST5', 'MAST6'],
    'SCLE': ['SCLE1', 'SCLE2', 'SCLE3', 'SCLE4', 'SCLE5'],
    'TC': ['TC1','TC1_2','TC2','TC3','TC4_TC8','TC5'],
    'LM': ['LM1','LM2','LM3','LM4','LM6','LM7','LM8'],
}


def extract_XYcoords_for_cells(sdata):
    sdata["table"].obs["centroidX"] = sd.get_centroids(sdata['cell_circles'])["x"]
    sdata["table"].obs["centroidY"] = sd.get_centroids(sdata['cell_circles'])["y"]
    sdata["table"].obsm["global_coords"] = np.column_stack([sdata["table"].obs["centroidX"], sdata["table"].obs["centroidY"]])

def load_adata_from_xenium(sample_name):
    sdata = xenium(xenium_paths[sample_name], morphology_focus=False, cells_as_circles=True)
    
    extract_XYcoords_for_cells(sdata)
    
    adata = sdata.tables['table']
    adata.obs['cell_id_unique'] = sample_name + '_' + adata.obs['cell_id']
    adata.obs['sample_name'] = sample_name
    adata.obs['sample_set'] = next(k for k, v in sample_sets.items() if sample_name in v)
    adata.obs['batch'] = studyid_to_batch[sample_name]
    adata.obs['slide'] = studyid_to_slide[sample_name]
    
    return adata


def add_boundaries_annotation_table(
    sdata,
    source_key="table",
    boundaries_element="cell_boundaries",
    dest_key="table_boundaries",
):
    """Add a second AnnData that annotates only ``cell_boundaries``.

    spatialdata-plot requires unique ``instance_key`` values per region inside a single table.
    Duplicating rows in one table (Option A) triggers ValueError during plotting; use this
    separate ``table_boundaries`` and pass ``table_name='table_boundaries'`` when rendering
    boundaries.
    """
    table = sdata[source_key]
    _, region_key, instance_key = get_table_keys(table)
    if boundaries_element not in sdata.shapes:
        raise KeyError(f"Missing shapes layer {boundaries_element!r} in sdata")
    gdf = sdata.shapes[boundaries_element]
    inst = table.obs[instance_key]
    if not inst.isin(gdf.index).all():
        n_miss = int((~inst.isin(gdf.index)).sum())
        raise ValueError(
            f"{n_miss} rows have {instance_key!r} values not in {boundaries_element!r} shape index"
        )
    tb = table.copy()
    tb.obs[region_key] = pd.Categorical(
        [boundaries_element] * tb.n_obs,
        categories=[boundaries_element],
    )
    tb = TableModel.parse(
        tb,
        region=boundaries_element,
        region_key=region_key,
        instance_key=instance_key,
        overwrite_metadata=True,
    )
    sdata[dest_key] = tb
    return sdata


def load_sdata_from_xenium(sample_name, merge_adata=None, add_he_image=True):
    sdata = xenium(xenium_paths[sample_name], morphology_focus=False, cells_as_circles=True)

    extract_XYcoords_for_cells(sdata)

    sdata["table"].obs["cell_id_unique"] = sample_name + "_" + sdata["table"].obs["cell_id"]
    sdata["table"].obs["sample_name"] = sample_name
    sdata["table"].obs["sample_set"] = next(k for k, v in sample_sets.items() if sample_name in v)
    sdata["table"].obs["batch"] = studyid_to_batch[sample_name]
    sdata["table"].obs["slide"] = studyid_to_slide[sample_name]

    if merge_adata is not None:
        new_cols = [c for c in merge_adata.obs.columns if c not in sdata["table"].obs.columns]
        if new_cols:
            if "cell_id_unique" in merge_adata.obs.columns:
                right = merge_adata.obs.set_index("cell_id_unique")[new_cols]
            else:
                right = merge_adata.obs[new_cols]
                right.index.name = "cell_id_unique"
            sdata["table"].obs = sdata["table"].obs.join(right, on="cell_id_unique")
    
    add_boundaries_annotation_table(sdata)

    if add_he_image:
        if sample_name in he_image_paths:
            image = xenium_aligned_image(he_image_paths[sample_name], he_alignment_paths[sample_name])
            sdata.images["he_stain"] = image
        else:
            print(f"Warning: no H&E image available for {sample_name}")

    return sdata
