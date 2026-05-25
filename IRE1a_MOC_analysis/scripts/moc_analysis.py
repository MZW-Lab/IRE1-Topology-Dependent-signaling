"""
moc_analysis.py
Cytoplasm mask extraction, background subtraction, Manders' Overlap Coefficient
(MOC) computation, and figure generation.

Must be run after segment.py has produced mask files in results/masks/.

Inputs
------
    data/4h/   -- original 3-channel TIFFs (C, Y, X), 16-bit
    data/24h/  -- original 3-channel TIFFs (C, Y, X), 16-bit
    results/masks/<stem>_nuc_labels.npy   -- from segment.py
    results/masks/<stem>_cell_labels.npy  -- from segment.py

Outputs
-------
    results/moc_results.csv   -- per-cell MOC values
    results/MOC_violin.pdf    -- violin + strip plot comparing 4h vs 24h

Usage
-----
    python scripts/moc_analysis.py
"""

import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from scipy import ndimage
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).parent))
from utils import load_tiff, bg_subtract, compute_moc

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT   = Path(__file__).parent.parent
DATA_DIRS   = {
    '4h':  REPO_ROOT / 'data' / '4h',
    '24h': REPO_ROOT / 'data' / '24h',
}
RESULTS_DIR = REPO_ROOT / 'results'
MASKS_DIR   = RESULTS_DIR / 'masks'

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
MIN_CYTO_AREA_PX = 200    # cells with fewer cytoplasm pixels are skipped
NUC_EROSION_PX   = 2      # pixels eroded from nucleus boundary before subtraction

# Channel indices (0-based)
CH_ER   = 0   # ER-Tracker Red
CH_IRE1 = 1   # FL-IRE1-mNG
CH_NUC  = 2   # NucBlue

# Figure aesthetics
C_4H  = '#AAAAAA'   # light gray
C_24H = '#444444'   # dark gray


# ---------------------------------------------------------------------------
# Cytoplasm mask extraction
# ---------------------------------------------------------------------------

def get_cytoplasm_masks(nuc_labels, cell_labels):
    """
    Return a dict mapping cell_id -> boolean cytoplasm mask.

    Cytoplasm = Cellpose cell body minus nucleus (with NUC_EROSION_PX erosion
    applied to the nucleus boundary to exclude ambiguous border pixels).
    """
    struct = ndimage.generate_binary_structure(2, 1)
    cyto_masks = {}
    for cell_id in np.unique(cell_labels):
        if cell_id == 0:
            continue
        cell_mask = cell_labels == cell_id

        # Find the nucleus that overlaps most with this cell
        nuc_ids, counts = np.unique(nuc_labels[cell_mask], return_counts=True)
        nuc_ids = nuc_ids[nuc_ids != 0]
        if nuc_ids.size == 0:
            nuc_mask = np.zeros_like(cell_mask)
        else:
            best_nuc = nuc_ids[np.argmax(counts[nuc_ids != 0])]
            nuc_mask = nuc_labels == best_nuc
            # Erode nucleus boundary
            for _ in range(NUC_EROSION_PX):
                nuc_mask = ndimage.binary_erosion(nuc_mask, structure=struct)

        cyto_mask = cell_mask & ~nuc_mask
        cyto_masks[cell_id] = cyto_mask

    return cyto_masks


# ---------------------------------------------------------------------------
# Per-image MOC computation
# ---------------------------------------------------------------------------

def process_image(tiff_path, timepoint):
    """
    Compute MOC for all cells in one image.

    Returns a list of dicts, one per valid cell.
    """
    stem = tiff_path.stem

    nuc_path  = MASKS_DIR / f"{stem}_nuc_labels.npy"
    cell_path = MASKS_DIR / f"{stem}_cell_labels.npy"

    if not nuc_path.exists() or not cell_path.exists():
        print(f"  WARNING: masks not found for {stem}, skipping.")
        return []

    img         = load_tiff(tiff_path)          # (C, Y, X)
    nuc_labels  = np.load(str(nuc_path))
    cell_labels = np.load(str(cell_path))

    cyto_masks = get_cytoplasm_masks(nuc_labels, cell_labels)

    # Build a combined cytoplasm mask for background estimation
    all_cyto = np.zeros(img.shape[1:], dtype=bool)
    for mask in cyto_masks.values():
        all_cyto |= mask

    # Background-subtract both channels
    ch_er_bg   = bg_subtract(img[CH_ER],   all_cyto)
    ch_ire1_bg = bg_subtract(img[CH_IRE1], all_cyto)

    records = []
    for cell_id, cyto_mask in cyto_masks.items():
        cyto_area = int(np.sum(cyto_mask))
        if cyto_area < MIN_CYTO_AREA_PX:
            print(f"    Skipping cell {cell_id} in {stem}: cytoplasm area {cyto_area} px < {MIN_CYTO_AREA_PX}")
            continue

        # Only pixels with non-zero intensity in both channels
        valid = cyto_mask & (ch_er_bg > 0) & (ch_ire1_bg > 0)
        if valid.sum() < 10:
            continue

        ch1_px = ch_er_bg[valid]
        ch2_px = ch_ire1_bg[valid]

        moc = compute_moc(ch1_px, ch2_px)

        records.append({
            'cell_id':      f"{stem}_{cell_id}",
            'image':        stem,
            'timepoint':    timepoint,
            'cyto_area_px': cyto_area,
            'n_pixels':     int(valid.sum()),
            'MOC':          moc,
        })

    print(f"  {stem}: {len(records)} cells quantified")
    return records


# ---------------------------------------------------------------------------
# Violin figure
# ---------------------------------------------------------------------------

def make_violin(df, out_path):
    mm_to_in = 1 / 25.4
    fig_w = 40 * mm_to_in
    fig_h = 34 * mm_to_in

    plt.rcParams.update({
        'font.family':        'DejaVu Sans',
        'font.size':           6,
        'axes.labelsize':      6,
        'axes.titlesize':      6.5,
        'xtick.labelsize':     6,
        'ytick.labelsize':     6,
        'pdf.fonttype':       42,
        'axes.spines.top':    False,
        'axes.spines.right':  False,
    })

    v4  = df.loc[df['timepoint'] == '4h',  'MOC'].dropna().values
    v24 = df.loc[df['timepoint'] == '24h', 'MOC'].dropna().values
    _, p = mannwhitneyu(v4, v24, alternative='two-sided')

    rng = np.random.default_rng(42)
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h), constrained_layout=True)

    for pos, color, vals in [(1, C_4H, v4), (2, C_24H, v24)]:
        vp = ax.violinplot([vals], positions=[pos], showmedians=True,
                           showextrema=True, widths=0.50)
        for body in vp['bodies']:
            body.set_facecolor(color); body.set_edgecolor(color)
            body.set_alpha(0.35);     body.set_linewidth(0.6)
        vp['cmedians'].set_color('black'); vp['cmedians'].set_linewidth(1.4)
        for part in ['cbars', 'cmins', 'cmaxes']:
            vp[part].set_color(color); vp[part].set_linewidth(0.7)
        jitter = rng.uniform(-0.09, 0.09, len(vals))
        ax.scatter(pos + jitter, vals, color=color, s=3.5, alpha=0.85,
                   zorder=3, linewidths=0)

    # Headroom
    ylo, yhi = ax.get_ylim()
    span = yhi - ylo
    ax.set_ylim(ylo, yhi + span * 0.28)
    ylo, yhi = ax.get_ylim()
    span = yhi - ylo

    # Median labels (to the right of each violin)
    for pos, vals in [(1, v4), (2, v24)]:
        med = np.median(vals)
        ax.text(pos + 0.28, med, f'{med:.3f}',
                ha='left', va='center', fontsize=5,
                fontweight='bold', color='black')

    # p-value bracket
    y_bar = ylo + span * 0.86
    ax.plot([1, 1, 2, 2],
            [y_bar - span * 0.022, y_bar, y_bar, y_bar - span * 0.022],
            color='black', lw=0.7)
    p_str = 'p < 0.0001' if p < 0.0001 else (
            'n.s.' if p >= 0.05 else f'p = {p:.4f}')
    ax.text(1.5, y_bar + span * 0.010, p_str,
            ha='center', va='bottom', fontsize=5.5)

    ax.set_xlim(0.30, 2.80)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'4 h\n(n={len(v4)})', f'24 h\n(n={len(v24)})'], fontsize=6)
    ax.set_ylabel('MOC', fontsize=6)
    ax.tick_params(axis='y', labelsize=6, length=2, pad=1)
    ax.tick_params(axis='x', length=2, pad=1)

    plt.savefig(str(out_path), bbox_inches='tight')
    plt.close()
    print(f"  Figure saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    all_records = []
    for timepoint, data_dir in DATA_DIRS.items():
        tiff_files = sorted(data_dir.glob('*.tif')) + sorted(data_dir.glob('*.tiff'))
        print(f"\n[{timepoint}] Processing {len(tiff_files)} images ...")
        for tiff_path in tiff_files:
            records = process_image(tiff_path, timepoint)
            all_records.extend(records)

    df = pd.DataFrame(all_records)
    csv_path = RESULTS_DIR / 'moc_results.csv'
    df.to_csv(str(csv_path), index=False)
    print(f"\nResults saved: {csv_path}  ({len(df)} cells)")

    # Summary
    for tp, grp in df.groupby('timepoint'):
        v = grp['MOC'].dropna()
        print(f"  {tp}: n={len(v)}, median MOC={v.median():.3f}, "
              f"IQR=[{v.quantile(0.25):.3f}, {v.quantile(0.75):.3f}]")

    # Figure
    fig_path = RESULTS_DIR / 'MOC_violin.pdf'
    make_violin(df, fig_path)
    print("\nDone.")


if __name__ == '__main__':
    main()
