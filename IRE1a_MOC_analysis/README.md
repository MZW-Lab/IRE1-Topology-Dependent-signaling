# IRE1-mNG / ER-Tracker Co-localization Analysis

Code for the nucleus segmentation, cell body segmentation, and Manders' Overlap Coefficient (MOC) analysis described in:

> **A topology-dependent IRE1α signaling axis generates a cryptic XBP1 isoform present in human cancers**
> Rodrigues Reyes JR, Vergel De Dios C, Datta D, Batjargal T, Pecchia S, Papa FR, Wilson MZ.
---

## Overview

This repository reproduces the image analysis pipeline used to quantify co-localization between FL-IRE1-mNG and ER-Tracker Red in SNB-19 cells treated with Thapsigargin at 4 h and 24 h post treatment

The pipeline has two steps:

1. **`segment.py`** -- Segments nuclei (Otsu threshold + watershed) and cell bodies (Cellpose CPSAM) from each TIFF.
2. **`moc_analysis.py`** -- Extracts cytoplasm masks, subtracts background, computes MOC per cell, and generates the violin plot figure.

---

## Requirements

- Python 3.9 or later
- See `requirements.txt` for package versions

Install dependencies:

```bash
pip install -r requirements.txt
```

Cellpose will download the CPSAM model weights on first run (~500 MB). An internet connection is required for the first run only.

---

## Input format

Each TIFF must be a 3-channel, 16-bit image with the following channel layout (0-indexed):

| Index | Channel        |
|-------|----------------|
| 0     | ER-Tracker Red |
| 1     | FL-IRE1a-mNG    |
| 2     | NucBlue        |

4-D files (T, C, Y, X) are accepted; only the first time frame is used.

Place images in:

```
data/4h/    <-- early timepoint TIFFs
data/24h/   <-- late timepoint TIFFs
```

Example images from the paper are included in these folders.

---

## Usage

Run from the repository root:

```bash
# Step 1: segmentation
python scripts/segment.py

# Step 2: MOC analysis and figure
python scripts/moc_analysis.py
```

---

## Outputs

All outputs are written to `results/` (created automatically):

| File | Description |
|------|-------------|
| `results/moc_results.csv` | Per-cell MOC values with metadata |
| `results/MOC_violin.pdf` | Violin + strip plot comparing 4 h vs 24 h (Mann-Whitney U) |
| `results/masks/` | Nucleus and cell body label arrays (.npy) |
| `results/segmentation_qc/` | Per-image QC overlay images (.png) |

---

## Key parameters

Parameters are defined as constants at the top of each script and can be edited directly:

| Parameter | Default | Location | Description |
|-----------|---------|----------|-------------|
| `DIAMETER_FACTOR` | 2.2 | `segment.py` | Cell diameter = factor x median nucleus diameter |
| `DIAMETER_MIN/MAX` | 80 / 200 px | `segment.py` | Cellpose diameter clamp range |
| `MIN_CYTO_AREA_PX` | 200 px | `moc_analysis.py` | Minimum cytoplasm area to include a cell |
| `NUC_EROSION_PX` | 2 px | `moc_analysis.py` | Nucleus boundary erosion before cytoplasm extraction |

---

## Software versions

| Package | Version |
|---------|---------|
| cellpose | 4.1.1 |
| numpy | 2.1.0 |
| scipy | 1.15.0 |
| scikit-image | 0.25.2 |
| tifffile | 2025.6.11 |
| matplotlib | 3.10.5 |
| pandas | 2.3.1 |

---

## License

[License placeholder]
