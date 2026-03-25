from __future__ import annotations
from typing import Dict, List, Union

CellTree = Dict[str, Union["CellTree", List[str]]]

CELL_TYPES_DETAILED: CellTree = {
    "Endothelial": {},
    "Fibroblast": {},
    "Keratinocytes": {},
    "Melanocyte": {},
    "Immune": {
        "B cell": [],
        "T cell": [
            "CD4+ T cells",
            "CD8+ T cells",
            "Tregs",
            "Gamma delta T cells",
            "NK cells",
            "T cell state",
        ],
        "Myeloid cell": [
            "Mast cells",
            "Macrophage",
            "Neutrophil",
            "Plasmacytoid DC",
            "CLEC10A+ DCs",
            "CLEC9A+ DCs",
            "LAMP3+ DCs",
            "CXCL9+ DCs",
            "Langerhans",
        ],
    },
}