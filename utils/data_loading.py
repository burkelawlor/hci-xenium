from pathlib import Path
import pandas as pd

import spatialdata as sd
from spatialdata_io import xenium
import numpy as np
import re


### Metadata and paths
metadata = pd.read_excel("./data/raw/metadata.xlsx", dtype=str).dropna(subset=["File path"])


### Xenium paths
data_root = Path("./data/raw/Xenium")
xenium_paths = {row['Study ID']: data_root / row['File path'] for index, row in metadata.iterrows()}


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
    'LP': ['LP1', 'LP2', 'LP3'],
    'LPP': ['LPP1', 'LPP2', 'LPP3', 'LPP4'],
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
