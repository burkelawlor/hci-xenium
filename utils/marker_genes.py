
"""
Marker genes for general cell types.

A lot of these came from Melissa's example notebooks. Other sources are noted below.
"""

# ---------------

marker_genes_general = {
    "Melanocyte": ["MLANA", "MITF", "SOX10"],
    "Endothelial": ["PECAM1",'KDR'],
    "Immune": ["CD3E", "CD4", "CD8A", "CD14", "CD68", "KIT", "CD207", "CLEC10A", "CLEC9A", "LAMP3", "CSF3R"],
    "B cells": ['CD19','MS4A1'],
    "Plasma cells":['XBP1','MZB1','IRF4'],
    "Plasmacytoid DC": ['IL3RA', 'CLEC4A'],
    "Mast cells": ['KIT'],
    "Fibroblast": ["PDGFRA", "PDGFRB",'THY1','FAP'], #THY1=CD90
    'Pericytes': ['PDGFRB','KCNJ8','ABCC9','NOTCH3'],
    "Keratinocytes": ["GRN", "DSG1", "TGM3"],
    "SOX10+": ["SOX10"],
    "Neural": ["MPZ", 'PYGB'],
    "Myocytes": ['CSRP1','MYLK','SYNM','SMTN'],
}

''' General marker gene sources
Neural (via claude)
- PYGB: a glycogen phosphorylase found predominantly in the brain
- APP: provides instructions for making a protein concentrated in neural synapses

Muscular (via Claude)
- CSRP1: crucial for cell development, differentiation, and structural integrity, particularly within the cytoskeleton
- MYLK: encodes an enzyme that regulates smooth muscle contraction, vascular permeability, and cell movement by phosphorylating myosin
- SYNM: plays a crucial role in maintaining the structural integrity of muscle cells
- SMTN: encodes a structural protein that is found exclusively in contractile smooth muscle cells.

'''

# -------------------------

marker_genes_immune_general = {
    "T cell": ["CD3E", "CD3G", "CD8A", "CD4", "TRDC", "FOXP3"],
    "NK cell": ["NCR1", "KLRB1", "NCAM1"], #NCR1 = NKp46, KLRB1 = NK1.1, NCAM1 = CD56
    "B cell": ["CD19",'MS4A1'],
    "Plasma cells":['XBP1','MZB1','IRF4'],
    "Macrophage": ["CD68", "CD14", "CD163", "MRC1"],
    "Langerhans": ["CD207"],
    "DC": ["CLEC10A", "CLEC9A", "LAMP3", "ITGAX", "CIITA", "CD1C", "CXCL10", "CXCL9"],
    "Mast": ["KIT"],
    "Neutrophils": ["MPO", "CSF3R", "FUT4"], #FUT4 = CD15
    "Plasmacytoid DC": ["IL3RA", "CLEC4A"],
    "Eosinophil": ["CEACAM8", "ITGA4"],
    "Basophil": ["ENPP3"]
}

# -------------------------

marker_genes_myeloid = {
    "Macrophage": ["CD68", "CD14", "CD163", "MRC1"],
    "Langerhans": ["CD207"],
    "DC": ["CLEC10A", "CLEC9A", "LAMP3", "ITGAX", "CIITA", "CD1C", "CXCL10", "CXCL9"],
    "Neutrophils": ["MPO", "CSF3R", "FUT4", "CXCR1"], #FUT4 = CD15
    "Plasmacytoid DC": ["IL3RA", "CLEC4A"],
    "Eosinophil": ["CEACAM8", "ITGA4"],
    "Basophil": ["ENPP3"],
    "Mast": ["KIT"],
}

# -------------------------

marker_genes_t_cell = {
    "CD8": ["CD8A", "CD8B", "GZMB", "PRF1"],
    "CD4": ["CD4", "FOXP3"],
    "T cell state": ["TCF7", "SELL", "TOX", "PDCD1", "HAVCR2"],
    "Tregs": ['CTLA4'],
    "Gamma delta T cells": ["TRDC"],
    "NK cells": ["NCR1", "KLRB1", "NCAM1"],
    'General T cells': ['CD3E','CD3G']
}

# -------------------------

marker_genes_endothelial = {
    'General endothelial': ['PECAM1','KDR'],
    'Vascular EC': ['MYLK','CD34','PLVAP'],
    'Lymphatic EC': ['PROX1','FLT4','PDPN']
}

''' Endothelial marker gene sources
He Y, Tacconi C, Dieterich LC, et al. Novel Blood Vascular Endothelial Subtype-Specific Markers in Human Skin Unearthed by Single-Cell Transcriptomic Profiling. Cells. 2022;11(7):1111. doi:10.3390/cells11071111
'''

# -------------------------

marker_genes_epithelial = {
    'Basal KC': ['COL17A1','ITGA3','ITGA6','CD46','DLL1','ITGA2', 'CDH3'],
    'Spinous KC': ['DSG1','CDH1'], 
    'Granular KC': ['TGM3','TGM1', 'CDSN'],
    'Sebocyte': ['FADS2','FASN'],
    "Melanocyte": ["MLANA", "MITF", "SOX10"],
}

''' Epithelial marker gene sources:

Negri VA, Watt FM. Understanding Human Epidermal Stem Cells at Single-Cell Resolution. J Invest Dermatol. 2022;142(8):2061-2067. doi:10.1016/j.jid.2022.04.003
- Proliferation: **CDK1**, **MKI67**
- Basal I: CAV1
- Basal II: **ITGA3**, **ITGA6**, **CD46**, **DLL1**  
- Basal III: **ITGA2**, ITGB1  
- Basal: **COL17A1**
- Transition I: **DUSP6**, DUSP10
- Transition II: DUSP10, PPTC7, **HES1**, FOS, JUN  
- Spinous: K1, K10
- Granular: FLG
- Immune: CD74

Reynolds G, Vegh P, Fletcher J, et al. Developmental cell programs are co-opted in inflammatory skin disease. Science. 2021;371(6527):eaba6500. doi:10.1126/science.aba6500
- Undifferentiated: KRT5, KRT14
- Proliferating: KRT1, KRT10
- Differentiated: 
- Inflammatory differentiated: **TP63**, **ITGA6**, ICAM1, **TNF**,  **CCL20**
- Basal: KRT5, KRT14

Wang S, Drummond ML, Guerrero-Juarez CF, et al. Single cell transcriptomics of human epidermis identifies basal stem cell transition states. Nat Commun. 2020;11(1):4239. doi:10.1038/s41467-020-18075-7
- Basal: KRT14, KRT5, **CDH3**
- Spinous: KRT1, KRT10, **DSG1**, **CDH1**
- Granular: DSC1, KRT2, IVL, **TGM3**

'''

# -------------------------

marker_genes_detailed = {
    "Fibroblasts": ["PDGFRA", "PDGFRB",'THY1'],
    'Pericytes': ['PDGFRB','KCNJ8','ABCC9','NOTCH3'],
    "Sweat duct": ["SOX10"],

    # Myeloid
    'Mast cells':['KIT'],
    'Neutrophils':['CSF3R'],
    'Macrophages':['CD68','CD14','CD163','MRC1'],
    'Langerhans':['CD207'],
    'CLEC10A+ DCs':['CLEC10A',],#'CIITA'
    'CLEC9A+ DCs':['CLEC9A',],#'CIITA'
    # 'CXCL9+ DCs':['CXCL9',],#'CIITA'
    'LAMP3+ DCs':['LAMP3',],#'CIITA'
    'Plasmacytoid DC':['IL3RA',],#'CIITA'

    # B cells
    'B cells': ['CD19'],
    'Plasma cells':['XBP1','MZB1','IRF4'],

    # T cells
    # 'T cells':['CD3E'],
    'CD8+ T cells':['CD8A',],#
    'CD4+ T cells':['CD4'],#'CD3E'
    'Tregs':['CTLA4'],#'CD3E'
    'Innate lymphocytes':['TRDC','CD3E'],#'CD3E'
    # 'NK cells':['KLRB1'],#'CD3E'

    # Endothelial
    'Vascular ECs': ['PECAM1','CD34','PLVAP'],
    'Lymphatic ECs': ['PECAM1','PROX1','FLT4','PDPN'],

    # Epithelial
    'Sebocytes':['FADS2','FASN'],
    'Basal KCs':['COL17A1','ITGA3'],#GRN
    'Spinous KCs':['DSG1','CDH1'],#GRN
    'Granular KCs':['TGM1','CDSN'],#GRN
    'Melanocytes': ["MLANA", "MITF", "SOX10"],
    
    'Unknown/dead': ['SOX2-OT'],
}