"""
Script 04: LFC Scatter Plot — XBP1-s vs XBP1-es (Fig 8)
=========================================================
Description:
    Generates a genome-wide scatter plot of log2 fold changes comparing
    XBP1-s vs WT (x-axis) and XBP1-es vs WT (y-axis) for all expressed genes.
    Highlights ISR/ATF4 target genes and canonical XBP1-s downstream targets.
    Quadrant labels indicate directional regulation relative to WT.

Input:
    data/DEG_results/XBP1s_vs_WT_full.csv
    data/DEG_results/XBP1es_vs_WT_full.csv

Output:
    figures/fig8_lfc_scatter_quadrants.pdf
    figures/fig8_lfc_scatter_quadrants.png

Requirements:
    Python >= 3.9
    pandas, numpy, matplotlib, adjustText
    Install: pip install pandas numpy matplotlib adjustText
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
from adjustText import adjust_text
import os

# ── Illustrator-compatible font settings ──────────────────────────────────────
matplotlib.rcParams.update({
    "pdf.fonttype": 42,       # TrueType fonts — editable in Illustrator
    "ps.fonttype":  42,
    "font.family":  "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size":    6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size":  2.5,
    "ytick.major.size":  2.5,
})

# ── Load DEG results ──────────────────────────────────────────────────────────
s_wt  = pd.read_csv("data/DEG_results/XBP1s_vs_WT_full.csv",
                    usecols=["gene_name", "log2FoldChange", "padj", "baseMean"])
es_wt = pd.read_csv("data/DEG_results/XBP1es_vs_WT_full.csv",
                    usecols=["gene_name", "log2FoldChange", "padj", "baseMean"])

# Deduplicate by highest baseMean (handles multi-mapping gene names)
s_dedup  = (s_wt.dropna(subset=["gene_name"])
                .sort_values("baseMean", ascending=False)
                .drop_duplicates("gene_name", keep="first"))
es_dedup = (es_wt.dropna(subset=["gene_name"])
                 .sort_values("baseMean", ascending=False)
                 .drop_duplicates("gene_name", keep="first"))

# Merge on gene_name
df = (s_dedup[["gene_name", "log2FoldChange", "padj"]]
      .merge(es_dedup[["gene_name", "log2FoldChange", "padj"]],
             on="gene_name", suffixes=("_s", "_es"))
      .dropna(subset=["log2FoldChange_s", "log2FoldChange_es"]))

print(f"Total genes in scatter: {len(df):,}")

# ── Gene sets to highlight ────────────────────────────────────────────────────
# ISR / ATF4 integrated stress response targets
isr_genes = [
    "ATF4", "ATF3", "DDIT3", "CHAC1", "ASNSP1", "SLC7A11",
    "SLC6A9", "SESN2", "DDIT4", "PMAIP1", "EIF5"
]

# Canonical XBP1-s downstream targets (ERAD / secretory pathway)
xbp1_targets = ["SELENOK", "DERL3", "SEC24D", "KDELR3"]

color_isr  = "#E05C2A"   # burnt orange
color_xbp1 = "#2A7DB5"  # steel blue

x = df["log2FoldChange_s"].values
y = df["log2FoldChange_es"].values

is_isr  = df["gene_name"].isin(isr_genes)
is_xbp1 = df["gene_name"].isin(xbp1_targets)
is_any  = is_isr | is_xbp1

# ── Axis limits ───────────────────────────────────────────────────────────────
xmin, xmax = -3.5, 6.0
ymin, ymax = -3.5, 6.0

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(3.2, 3.2))

# Background genes (clipped to axis range)
x_bg = np.clip(x[~is_any], xmin, xmax)
y_bg = np.clip(y[~is_any], ymin, ymax)
ax.scatter(x_bg, y_bg, s=1.5, color="#CCCCCC", alpha=0.45,
           linewidths=0, rasterized=True, zorder=1)

# Reference lines
ax.plot([xmin, xmax], [ymin, ymax], color="#999999", lw=0.5, ls="--", zorder=0)
ax.axhline(0, color="#AAAAAA", lw=0.4, ls=":", zorder=0)
ax.axvline(0, color="#AAAAAA", lw=0.4, ls=":", zorder=0)

# Highlighted gene sets
sub_isr  = df[is_isr].copy()
sub_xbp1 = df[is_xbp1].copy()
for sub in [sub_isr, sub_xbp1]:
    sub["log2FoldChange_s"]  = sub["log2FoldChange_s"].clip(xmin, xmax)
    sub["log2FoldChange_es"] = sub["log2FoldChange_es"].clip(ymin, ymax)

ax.scatter(sub_isr["log2FoldChange_s"],  sub_isr["log2FoldChange_es"],
           s=18, color=color_isr,  linewidths=0.4, edgecolors="white",
           zorder=4, label="ISR / ATF4 targets")
ax.scatter(sub_xbp1["log2FoldChange_s"], sub_xbp1["log2FoldChange_es"],
           s=18, color=color_xbp1, linewidths=0.4, edgecolors="white",
           zorder=4, label="XBP1-s canonical targets")

# Gene labels with auto-adjustment to avoid overlap
texts = []
label_df = pd.concat([sub_isr, sub_xbp1])
for _, row in label_df.iterrows():
    color = color_isr if row["gene_name"] in isr_genes else color_xbp1
    t = ax.text(row["log2FoldChange_s"], row["log2FoldChange_es"],
                row["gene_name"], fontsize=5.5, color=color,
                fontweight="bold", ha="left", va="bottom")
    texts.append(t)

adjust_text(
    texts,
    x=label_df["log2FoldChange_s"].values,
    y=label_df["log2FoldChange_es"].values,
    ax=ax,
    arrowprops=dict(arrowstyle="-", color="#888888", lw=0.4),
    expand_points=(1.5, 1.5),
    expand_text=(1.3, 1.3),
    force_points=(0.35, 0.35),
    force_text=(0.6, 0.6),
    lim=600,
)

# ── Quadrant labels ───────────────────────────────────────────────────────────
x_zero_frac = (0 - xmin) / (xmax - xmin)
y_zero_frac = (0 - ymin) / (ymax - ymin)

quad_kw = dict(fontsize=5.5, color="#888888", style="italic",
               ha="center", va="center", transform=ax.transAxes)
ax.text((x_zero_frac + 1) / 2,      (y_zero_frac + 1) / 2 + 0.05, "Up in both",     **quad_kw)
ax.text( x_zero_frac / 2,           (y_zero_frac + 1) / 2 + 0.05, "Up in XBP1-es",  **quad_kw)
ax.text((x_zero_frac + 1) / 2,       y_zero_frac / 2 - 0.04,      "Up in XBP1-s",   **quad_kw)
ax.text( x_zero_frac / 2,            y_zero_frac / 2 - 0.04,      "Down in both",    **quad_kw)

# ── Pearson r annotation ──────────────────────────────────────────────────────
r = np.corrcoef(x, y)[0, 1]
ax.text(0.04, 0.60, f"r = {r:.3f}  n = {len(df):,}",
        transform=ax.transAxes, fontsize=5.5, va="top", color="#555555")

# ── Legend ────────────────────────────────────────────────────────────────────
leg = ax.legend(loc="lower right", fontsize=5, frameon=True, framealpha=0.85,
                edgecolor="#CCCCCC", handlelength=0.8, handletextpad=0.4,
                borderpad=0.5, labelspacing=0.35, markerscale=0.9,
                bbox_to_anchor=(0.99, 0.01))
leg.get_frame().set_linewidth(0.4)

# ── Axes ──────────────────────────────────────────────────────────────────────
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks([t for t in range(-3, 7) if xmin <= t <= xmax])
ax.set_yticks([t for t in range(-3, 7) if ymin <= t <= ymax])
ax.set_xlabel("log\u2082 fold change  (XBP1-s vs WT)", fontsize=6)
ax.set_ylabel("log\u2082 fold change  (XBP1-es vs WT)", fontsize=6)
ax.tick_params(labelsize=5.5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout(pad=0.4)

# ── Save ──────────────────────────────────────────────────────────────────────
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/fig8_lfc_scatter_quadrants.pdf", dpi=300, bbox_inches="tight")
fig.savefig("figures/fig8_lfc_scatter_quadrants.png", dpi=300, bbox_inches="tight")
plt.close()
print("Saved: figures/fig8_lfc_scatter_quadrants.pdf")
