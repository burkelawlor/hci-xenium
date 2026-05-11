from pathlib import Path
import pandas as pd

import spatialdata as sd
from spatialdata_io import xenium
import numpy as np
import re


### Metadata and paths

metadata = pd.read_excel("./data/metadata.xlsx", dtype=str).dropna(subset=["GnomEx ID"])

batch_1 = ['TC1', 'SCLE1', 'SCLE2', 'PNT1', 'PNT2', 'EPD1', 'EPD2', 'EPD3', 'EPD4', 'EPD5', 'EPD6', 'EPD7']
batch_2 = ['NS1', 'NS2', 'NS3', 'NS4', 'SCLE3', 'SCLE4', 'SCLE5', 'EPD8', 'EPD9', 'EPD10', 'LPP1']
batch_3 = ['LPP2', 'LPP3', 'LPP4', 'LP1', 'LP2', 'LP3', 'OLP1', 'OLP2', 'MAST1', 'MAST2', 'MAST3', 'OLP3', 'OLP4', 'OLP5', 'MAST6', 'OLP6']

data_root = Path("./data/Xenium")
batch_1_root = data_root / '26697R'
batch_2_root = data_root / '20260401__191124__28129'
batch_3_root = data_root / '20260506__195908__28200'

xenium_paths = {
    "TC1":   batch_1_root / "output-XETG00516__0076657__Region_6__20251017__213834",
    "SCLE1": batch_1_root / "output-XETG00516__0076678__Region_1__20251017__213834",
    "SCLE2": batch_1_root / "output-XETG00516__0076678__Region_2__20251017__213834",
    "PNT1":  batch_1_root / "output-XETG00516__0076678__Region_4__20251017__213834",
    "PNT2":  batch_1_root / "output-XETG00516__0076678__Region_3__20251017__213834",
    "EPD1":  batch_1_root / "output-XETG00516__0076657__Region_5__20251017__213834",
    "EPD2":  batch_1_root / "output-XETG00516__0076657__Region_2__20251017__213834",
    "EPD3":  batch_1_root / "output-XETG00516__0076678__Region_5__20251017__213834",
    "EPD4":  batch_1_root / "output-XETG00516__0076678__Region_6__20251017__213834",
    "EPD5":  batch_1_root / "output-XETG00516__0076657__Region_3__20251017__213834",
    "EPD6":  batch_1_root / "output-XETG00516__0076657__Region_1__20251017__213833",
    "EPD7":  batch_1_root / "output-XETG00516__0076657__Region_4__20251017__213834",

    'NS1':   batch_2_root / 'output-XETG00516__0084912__Region_2__20260401__191409',
    'NS2':   batch_2_root / 'output-XETG00516__0084920__Region_2__20260401__191409',
    'NS3':   batch_2_root / 'output-XETG00516__0084912__Region_3__20260401__191409',
    'NS4':   batch_2_root / 'output-XETG00516__0084912__Region_5__20260401__191409',
    'SCLE3': batch_2_root / 'output-XETG00516__0084920__Region_3__20260401__191409',
    'SCLE4': batch_2_root / 'output-XETG00516__0084912__Region_4__20260401__191409',
    'SCLE5': batch_2_root / 'output-XETG00516__0084920__Region_1__20260401__191409',
    'EPD8':  batch_2_root / 'output-XETG00516__0084912__Region_1__20260401__191409',
    'EPD9':  batch_2_root / 'output-XETG00516__0084920__Region_5__20260401__191409',
    'EPD10': batch_2_root / 'output-XETG00516__0084920__Region_4__20260401__191409',
    'LPP1':  batch_2_root / 'output-XETG00516__0084912__Region_6__20260401__191409',

    'LPP2':  batch_3_root / 'output-XETG00516__0076414__Region_7__20260506__200316',
    'LPP3':  batch_3_root / 'output-XETG00516__0076660__Region_8__20260506__200317',
    'LPP4':  batch_3_root / 'output-XETG00516__0076414__Region_3__20260506__200316',
    'LP1':   batch_3_root / 'output-XETG00516__0076414__Region_4__20260506__200316',
    'LP2':   batch_3_root / 'output-XETG00516__0076414__Region_4__20260506__200316',
    'LP3':   batch_3_root / 'output-XETG00516__0076414__Region_5__20260506__200316',
    'OLP1':  batch_3_root / 'output-XETG00516__0076660__Region_7__20260506__200317',
    'OLP2':  batch_3_root / 'output-XETG00516__0076414__Region_6__20260506__200316',
    'MAST1': batch_3_root / 'output-XETG00516__0076660__Region_5__20260506__200317',
    'MAST2': batch_3_root / 'output-XETG00516__0076660__Region_3__20260506__200317',
    'MAST3': batch_3_root / 'output-XETG00516__0076660__Region_4__20260506__200317',
    'OLP3':  batch_3_root / 'output-XETG00516__0076660__Region_6__20260506__200317',
    'OLP4':  batch_3_root / 'output-XETG00516__0076660__Region_1__20260506__200316',
    'OLP5':  batch_3_root / 'output-XETG00516__0076414__Region_1__20260506__200316',
    'MAST6': batch_3_root / 'output-XETG00516__0076414__Region_2__20260506__200316',
    'OLP6':  batch_3_root / 'output-XETG00516__0076660__Region_2__20260506__200317',
}



### Functions for loading Xenium data

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
    adata.obs['sample_set'] = re.sub(r'[^a-zA-Z]', '', sample_name)
    adata.obs['batch'] = 'batch_1' if sample_name in batch_1 else 'batch_2' if sample_name in batch_2 else 'batch_3' if sample_name in batch_3 else None
    
    return adata