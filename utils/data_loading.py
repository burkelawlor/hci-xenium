from pathlib import Path

batch_1_root = Path("./data/Xenium/26697R")
batch_2_root = Path("./data/Xenium/28129")

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
}