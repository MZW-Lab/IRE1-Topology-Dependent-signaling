"""
segment.py
Nucleus and cell body segmentation for the IRE1-mNG / ER-Tracker co-localization pipeline.

Inputs
------
    data/4h/   -- 3-channel TIFFs (C, Y, X), 16-bit
    data/24h/  -- 3-channel TIFFs (C, Y, X), 16-bit

    Channel layout (0-indexed):
        CH0 = ER-Tracker Red
        CH1 = FL-IRE1-mNG
        CH2 = NucBlue

Outputs (written to results/)
------
    results/masks/<image_stem>_nuc_labels.npy    -- nucleus label array (int32)
    results/masks/<image_stem>_cell_labels.npy   -- Cellpose cell body label array (int32)
    results/segmentation_qc/<image_stem>_qc.png  -- overlay QC image

Usage
-----
    python scripts/segment.py
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

# Allow imports from scripts/ regardless of working directory
sys.path.insert(0, str(Path(__file__).parent))
from utils import load_tiff, nucleus_segment

# ---------------------------------------------------------------------------
# Paths -- edit DATA_DIRS and RESULTS_DIR if your folder layout differs
# ---------------------------------------------------------------------------
REPO_ROOT   = Path(__file__).parent.parent
DATA_DIRS   = {
    '4h':  REPO_ROOT / 'data' / '4h',
    '24h': REPO_ROOT / 'data' / '24h',
}
RESULTS_DIR = REPO_ROOT / 'results'
MASKS_DIR   = RESULTS_DIR / 'masks'
QC_DIR      = RESULTS_DIR / 'segmentation_qc'

# ---------------------------------------------------------------------------
# Cellpose parameters
# ---------------------------------------------------------------------------
CELLPOSE_MODEL    = 'cpsam'          # Cellpose CPSAM model
DIAMETER_FACTOR   = 2.2              # cell diameter = DIAMETER_FACTOR * median nucleus diameter
DIAMETER_MIN      = 80               # px
DIAMETER_MAX      = 200              # px
MIN_CYTO_AREA_PX  = 200             # cells with smaller cytoplasm are skipped in moc_analysis.py

# ---------------------------------------------------------------------------
# Channel indices (0-based)
# ---------------------------------------------------------------------------
CH_ER    = 0   # ER-Tracker Red
CH_IRE1  = 1   # FL-IRE1-mNG
CH_NUC   = 2   # NucBlue


def estimate_cell_diameter(nuc_labels):
    """Estimate cell diameter as DIAMETER_FACTOR * median nucleus equivalent diameter."""
    areas = []
    for lbl in np.unique(nuc_labels):
        if lbl == 0:
            continue
        areas.append(np.sum(nuc_labels == lbl))
    if not areas:
        return (DIAMETER_MIN + DIAMETER_MAX) // 2
    median_area = np.median(areas)
    nuc_diam    = 2 * np.sqrt(median_area / np.pi)
    cell_diam   = DIAMETER_FACTOR * nuc_diam
    return float(np.clip(cell_diam, DIAMETER_MIN, DIAMETER_MAX))


def run_cellpose(ch_ire1, ch_nuc, diameter):
    """Run Cellpose CPSAM on a two-channel input (IRE1 + NucBlue)."""
    from cellpose import models
    model = models.CellposeModel(gpu=False, model_type=CELLPOSE_MODEL)
    # Stack channels: Cellpose expects (Y, X, C) or list; use two-channel input
    img = np.stack([ch_ire1, ch_nuc], axis=0)   # (2, Y, X)
    masks, _, _ = model.eval(
        img,
        diameter=diameter,
        channels=[1, 2],   # channel 1 = cytoplasm (IRE1), channel 2 = nucleus (NucBlue)
        flow_threshold=0.4,
        cellprob_threshold=0.0,
    )
    return masks.astype(np.int32)


def save_qc_overlay(img_cyx, nuc_labels, cell_labels, out_path):
    """Save a 3-panel QC overlay: CH1 (IRE1) | nucleus outlines | cell outlines."""
    from skimage.segmentation import find_boundaries

    ch_ire1 = img_cyx[CH_IRE1].astype(np.float32)
    ch_ire1 = (ch_ire1 - ch_ire1.min()) / (ch_ire1.max() - ch_ire1.min() + 1e-8)

    nuc_bounds  = find_boundaries(nuc_labels,  mode='outer')
    cell_bounds = find_boundaries(cell_labels, mode='outer')

    fig, axes = plt.subplots(1, 3, figsize=(9, 3))
    axes[0].imshow(ch_ire1, cmap='gray');         axes[0].set_title('FL-IRE1-mNG (CH1)', fontsize=7)
    axes[1].imshow(ch_ire1, cmap='gray')
    axes[1].contour(nuc_labels,  colors='cyan',   linewidths=0.5)
    axes[1].set_title('Nucleus segmentation', fontsize=7)
    axes[2].imshow(ch_ire1, cmap='gray')
    axes[2].contour(cell_labels, colors='yellow', linewidths=0.5)
    axes[2].contour(nuc_labels,  colors='cyan',   linewidths=0.5)
    axes[2].set_title('Cell body segmentation', fontsize=7)

    for ax in axes:
        ax.axis('off')

    plt.rcParams.update({'pdf.fonttype': 42, 'font.family': 'DejaVu Sans'})
    plt.tight_layout(pad=0.3)
    plt.savefig(str(out_path), dpi=150, bbox_inches='tight')
    plt.close()


def process_image(tiff_path, timepoint):
    stem = tiff_path.stem
    print(f"  Processing {tiff_path.name} ...")

    img = load_tiff(tiff_path)   # (C, Y, X)

    # 1. Nucleus segmentation from CH2 (NucBlue)
    nuc_labels = nucleus_segment(img[CH_NUC])
    n_nuclei   = nuc_labels.max()
    print(f"    Nuclei detected: {n_nuclei}")

    # 2. Estimate cell diameter from nucleus size
    diameter = estimate_cell_diameter(nuc_labels)
    print(f"    Cellpose diameter: {diameter:.1f} px")

    # 3. Cellpose cell body segmentation
    cell_labels = run_cellpose(img[CH_IRE1], img[CH_NUC], diameter)
    n_cells     = cell_labels.max()
    print(f"    Cells detected:   {n_cells}")

    # 4. Save masks
    np.save(str(MASKS_DIR / f"{stem}_nuc_labels.npy"),  nuc_labels)
    np.save(str(MASKS_DIR / f"{stem}_cell_labels.npy"), cell_labels)

    # 5. Save QC overlay
    save_qc_overlay(img, nuc_labels, cell_labels,
                    QC_DIR / f"{stem}_qc.png")

    return n_cells


def main():
    MASKS_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)

    total_cells = 0
    for timepoint, data_dir in DATA_DIRS.items():
        tiff_files = sorted(data_dir.glob('*.tif')) + sorted(data_dir.glob('*.tiff'))
        print(f"\n[{timepoint}] Found {len(tiff_files)} TIFF files in {data_dir}")
        for tiff_path in tiff_files:
            n = process_image(tiff_path, timepoint)
            total_cells += n

    print(f"\nSegmentation complete. Total cells detected across all images: {total_cells}")
    print(f"Masks saved to:  {MASKS_DIR}")
    print(f"QC images saved: {QC_DIR}")


if __name__ == '__main__':
    main()
