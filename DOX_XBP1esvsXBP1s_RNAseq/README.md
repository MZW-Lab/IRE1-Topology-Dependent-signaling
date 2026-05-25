# Supplementary Data — XBP1-s vs XBP1-es Bulk RNA-seq Analysis

This repository contains all data, code, and figures associated with the bulk RNA-seq differential expression analysis comparing XBP1-s and XBP1-es isoforms against wild-type (WT) controls.

---

## Repository Structure

```
supplement/
├── README.md                          # This file
├── data/
│   ├── raw/
│   │   ├── expression_matrix_raw_counts_cpm.tsv   # Raw count matrix (all 9 samples)
│   │   └── sample_metadata_runB.csv               # Sample metadata (Run B)
│   ├── DEG_results/
│   │   ├── XBP1s_vs_WT_full.csv                   # All genes, XBP1-s vs WT
│   │   ├── XBP1s_vs_WT_significant.csv            # Significant DEGs, XBP1-s vs WT
│   │   ├── XBP1es_vs_WT_full.csv                  # All genes, XBP1-es vs WT
│   │   ├── XBP1es_vs_WT_significant.csv           # Significant DEGs, XBP1-es vs WT
│   │   ├── XBP1s_vs_XBP1es_full.csv               # All genes, XBP1-s vs XBP1-es
│   │   └── XBP1s_vs_XBP1es_significant.csv        # Significant DEGs, XBP1-s vs XBP1-es
│   └── GSEA_results/
│       ├── XBP1s_vs_WT_gsea_results.csv
│       ├── XBP1es_vs_WT_gsea_results.csv
│       └── XBP1s_vs_XBP1es_gsea_results.csv
├── figures/
│   ├── fig1_pca.pdf                               # PCA of all 9 samples
│   ├── fig2a_volcano_XBP1s_vs_WT.pdf             # Volcano: XBP1-s vs WT
│   ├── fig2b_volcano_XBP1es_vs_WT.pdf            # Volcano: XBP1-es vs WT
│   ├── fig2c_volcano_XBP1s_vs_XBP1es.pdf         # Volcano: XBP1-s vs XBP1-es
│   ├── fig3_heatmap_XBP1s_vs_XBP1es_DEGs.pdf     # Heatmap: 22 isoform-specific DEGs
│   ├── fig4a_gsea_XBP1s_vs_WT.pdf                # GSEA dotplot: XBP1-s vs WT
│   ├── fig4b_gsea_XBP1es_vs_WT.pdf               # GSEA dotplot: XBP1-es vs WT
│   ├── fig4c_gsea_XBP1s_vs_XBP1es.pdf            # GSEA dotplot: XBP1-s vs XBP1-es
│   ├── fig5a_spliceosome_es_only.pdf              # Spliceosome heatmap: XBP1-es-specific
│   ├── fig5b_spliceosome_s_only.pdf               # Spliceosome heatmap: XBP1-s-specific
│   ├── fig5c_spliceosome_both_up.pdf              # Spliceosome heatmap: both up
│   ├── fig6_spliceosome_scatter.pdf               # Spliceosome gene scatter
│   ├── fig7_genome_scatter.pdf                    # Genome-wide LFC scatter (significant genes)
│   └── fig8_lfc_scatter_quadrants.pdf             # Genome-wide LFC scatter (all genes, quadrant labels)
└── code/
    ├── 01_DESeq2_differential_expression.R        # DESeq2 pipeline (Run B)
    ├── 02_GSEA_pathway_analysis.R                 # GSEA (Hallmark + KEGG)
    ├── 03_figures_PCA_volcano_heatmap.R           # PCA, volcano, heatmap figures
    └── 04_figure_lfc_scatter.py                   # LFC scatter plot (Python)
```

---

## Experimental Design

| Condition     | Replicates | Description                              |
|---------------|-----------|------------------------------------------|
| WT + DOX      | 3         | Wild-type cells treated with doxycycline |
| XBP1-s + DOX  | 3         | XBP1-s-expressing cells + doxycycline   |
| XBP1-es + DOX | 3         | XBP1-es-expressing cells + doxycycline  |

**Total samples**: 9 (3 conditions × 3 replicates)

---

## Run B: Pseudo-replicate Rationale

XBP1-s replicate 2 was a failed sequencing run (~19,000 total counts vs 5.4–8.7 million for all other samples). A re-sequenced library was provided but remained a transcriptomic outlier:

- PC1 score: 25.1 (vs 0.0–1.6 for group mates)
- Within-group Pearson r: 0.935–0.943 (vs 0.988–0.999 for all other within-group pairs)

**Run B** replaces this sample with a pseudo-replicate computed as:

```
XBP1s_DOX_2_pseudo = floor((XBP1s_DOX_1 + XBP1s_DOX_3) / 2)
```

This approach is conservative (does not add independent biological information) and is transparently documented here. Run A results (with the re-sequenced library retained) are available upon request for robustness comparison.

---

## Analysis Pipeline

### 1. Differential Expression (DESeq2)

- **Pre-filtering**: genes with ≥ 10 counts in ≥ 3 samples retained (15,777 genes)
- **Model**: `~ condition` (WT_DOX as reference)
- **LFC shrinkage**:
  - XBP1-s vs WT, XBP1-es vs WT: `apeglm`
  - XBP1-s vs XBP1-es: `ashr` (contrast-based)
- **Significance thresholds**:

| Comparison          | \|log2FC\| threshold | padj threshold | DEGs |
|---------------------|---------------------|----------------|------|
| XBP1-s vs WT        | ≥ 1.0               | ≤ 0.05         | 734  |
| XBP1-es vs WT       | ≥ 1.0               | ≤ 0.05         | 620  |
| XBP1-s vs XBP1-es   | ≥ 0.5               | ≤ 0.05         | 22   |

### 2. Gene Set Enrichment Analysis (GSEA)

- **Gene sets**: MSigDB Hallmark (50 sets) + KEGG (186 sets) via `msigdbr`
- **Ranking metric**: sign(log2FC) × −log10(pvalue)
- **Parameters**: minGSSize = 15, maxGSSize = 500, padj method = BH
- **Tool**: `clusterProfiler::GSEA()`

### 3. Figures

All PDF figures use TrueType fonts (pdf.fonttype = 42 / cairo_pdf) and are fully editable in Adobe Illustrator. Font size: 5.5–6 pt throughout.

---

## Key Findings

1. **Genome-wide similarity**: XBP1-s and XBP1-es drive nearly identical transcriptional programs (Pearson r = 0.943 across 5,837 significant genes), both strongly upregulating ERAD, UPR, and secretory pathway genes.

2. **ISR/ATF4 targets preferentially induced by XBP1-s**: Among the 22 isoform-specific DEGs, ATF4, ATF3, DDIT3/CHOP, CHAC1, SESN2, SLC7A11, and DDIT4 are significantly higher in XBP1-s vs XBP1-es, suggesting XBP1-s more strongly activates the integrated stress response arm.

3. **Spliceosome regulation**: XBP1-s preferentially suppresses SR proteins (SRSF1/2/3/6/7) and hnRNP/RBP genes (GSEA NES = −2.28, padj = 1.2×10⁻⁷), while XBP1-es preferentially upregulates NTC/Prp19 and U4/U6 snRNP components.

4. **JUN is XBP1-s-specific**: JUN is significantly upregulated in XBP1-s but unchanged in XBP1-es, consistent with stronger activation of stress-adaptive AP-1 programs.

---

## Software Versions

| Software       | Version  |
|----------------|----------|
| R              | ≥ 4.2    |
| DESeq2         | ≥ 1.38   |
| apeglm         | ≥ 1.20   |
| ashr           | ≥ 2.2    |
| clusterProfiler| ≥ 4.6    |
| msigdbr        | ≥ 7.5    |
| ComplexHeatmap | ≥ 2.14   |
| ggplot2        | ≥ 3.4    |
| Python         | ≥ 3.9    |
| pandas         | ≥ 1.5    |
| matplotlib     | ≥ 3.6    |
| adjustText     | ≥ 0.8    |

---

## Running the Code

Scripts should be run from the `supplement/` directory in order:

```bash
# R scripts (run in R or RStudio)
Rscript code/01_DESeq2_differential_expression.R
Rscript code/02_GSEA_pathway_analysis.R
Rscript code/03_figures_PCA_volcano_heatmap.R

# Python script
pip install pandas numpy matplotlib adjustText
python code/04_figure_lfc_scatter.py
```

> **Note**: Scripts 02 and 03 source script 01 internally. Ensure all dependencies are installed before running.

---

## Contact

For questions regarding this analysis, please contact the corresponding author.
