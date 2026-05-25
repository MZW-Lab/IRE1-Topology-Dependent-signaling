# =============================================================================
# Script 03: Publication Figures — PCA, Volcano Plots, DEG Heatmap
# =============================================================================
# Description:
#   Generates publication-ready figures:
#     Fig 1  — PCA of all 9 Run B samples
#     Fig 2a — Volcano plot: XBP1-s vs WT
#     Fig 2b — Volcano plot: XBP1-es vs WT
#     Fig 2c — Volcano plot: XBP1-s vs XBP1-es
#     Fig 3  — Heatmap of 22 XBP1-s vs XBP1-es DEGs
#
# Input:
#   data/raw/expression_matrix_raw_counts_cpm.tsv
#   data/raw/sample_metadata_runB.csv
#   data/DEG_results/*.csv
#
# Output:
#   figures/fig1_pca.pdf
#   figures/fig2a_volcano_XBP1s_vs_WT.pdf
#   figures/fig2b_volcano_XBP1es_vs_WT.pdf
#   figures/fig2c_volcano_XBP1s_vs_XBP1es.pdf
#   figures/fig3_heatmap_XBP1s_vs_XBP1es_DEGs.pdf
#
# Software: R >= 4.2, DESeq2, ggplot2, ggrepel, ComplexHeatmap, circlize
# =============================================================================

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tibble)

# Font settings for Illustrator compatibility
# Use cairo_pdf device for editable text (equivalent to pdf.fonttype=42 in Python)

# ── Shared theme ──────────────────────────────────────────────────────────────
theme_pub <- function(base_size = 6) {
  theme_classic(base_size = base_size) +
    theme(
      text          = element_text(family = "sans", size = base_size),
      axis.text     = element_text(size = base_size - 0.5),
      axis.title    = element_text(size = base_size),
      legend.text   = element_text(size = base_size - 1),
      legend.title  = element_text(size = base_size),
      plot.title    = element_text(size = base_size + 0.5, face = "bold"),
      axis.line     = element_line(linewidth = 0.4),
      axis.ticks    = element_line(linewidth = 0.4)
    )
}

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# ── Load data and rebuild DESeq2 object ───────────────────────────────────────
# (Re-run steps from script 01 to get dds and rld objects)
source("01_DESeq2_differential_expression.R")  # loads dds, count_mat_B, meta

# Regularized log transform for PCA/heatmap
rld <- rlog(dds, blind = FALSE)

# ── Fig 1: PCA ────────────────────────────────────────────────────────────────
pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

condition_colors <- c(
  "WT_DOX"     = "#555555",
  "XBP1s_DOX"  = "#2A7DB5",
  "XBP1es_DOX" = "#E05C2A"
)
condition_labels <- c(
  "WT_DOX"     = "WT + DOX",
  "XBP1s_DOX"  = "XBP1-s + DOX",
  "XBP1es_DOX" = "XBP1-es + DOX"
)

p_pca <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = condition_colors, labels = condition_labels,
                     name = NULL) +
  labs(x = paste0("PC1 (", pct_var[1], "%)"),
       y = paste0("PC2 (", pct_var[2], "%)")) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"))

ggsave("figures/fig1_pca.pdf", p_pca,
       width = 2.8, height = 2.8, device = cairo_pdf)

# ── Fig 2: Volcano plots ──────────────────────────────────────────────────────
make_volcano <- function(df, title,
                         lfc_thresh, padj_thresh,
                         label_genes = NULL,
                         xlim = c(-6, 6)) {

  df <- df %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(
      neg_log10_padj = pmin(-log10(padj), 50),
      sig = case_when(
        padj <= padj_thresh & log2FoldChange >=  lfc_thresh ~ "Up",
        padj <= padj_thresh & log2FoldChange <= -lfc_thresh ~ "Down",
        TRUE ~ "NS"
      )
    )

  sig_colors <- c("Up" = "#E05C2A", "Down" = "#2A7DB5", "NS" = "#CCCCCC")

  p <- ggplot(df, aes(log2FoldChange, neg_log10_padj, color = sig)) +
    geom_point(size = 0.6, alpha = 0.6, stroke = 0) +
    scale_color_manual(values = sig_colors, guide = "none") +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
               lty = 2, lwd = 0.3, color = "grey50") +
    geom_hline(yintercept = -log10(padj_thresh),
               lty = 2, lwd = 0.3, color = "grey50") +
    xlim(xlim) +
    labs(title = title,
         x = "log\u2082 fold change",
         y = "-log\u2081\u2080(padj)") +
    theme_pub()

  if (!is.null(label_genes)) {
    label_df <- df %>% filter(gene_name %in% label_genes, sig != "NS")
    p <- p + geom_text_repel(
      data = label_df,
      aes(label = gene_name),
      size = 1.8, max.overlaps = 20,
      segment.size = 0.2, segment.color = "grey50",
      color = "black", fontface = "bold"
    )
  }
  p
}

s_wt  <- read_csv("data/DEG_results/XBP1s_vs_WT_full.csv")
es_wt <- read_csv("data/DEG_results/XBP1es_vs_WT_full.csv")
s_es  <- read_csv("data/DEG_results/XBP1s_vs_XBP1es_full.csv")

highlight_genes <- c("ATF4", "ATF3", "DDIT3", "HSPA5", "DERL3",
                     "SEC24D", "SELENOK", "KDELR3", "JUN", "SRSF1")

ggsave("figures/fig2a_volcano_XBP1s_vs_WT.pdf",
       make_volcano(s_wt,  "XBP1-s vs WT",
                   lfc_thresh = 1.0, padj_thresh = 0.05,
                   label_genes = highlight_genes),
       width = 2.8, height = 3.0, device = cairo_pdf)

ggsave("figures/fig2b_volcano_XBP1es_vs_WT.pdf",
       make_volcano(es_wt, "XBP1-es vs WT",
                   lfc_thresh = 1.0, padj_thresh = 0.05,
                   label_genes = highlight_genes),
       width = 2.8, height = 3.0, device = cairo_pdf)

ggsave("figures/fig2c_volcano_XBP1s_vs_XBP1es.pdf",
       make_volcano(s_es,  "XBP1-s vs XBP1-es",
                   lfc_thresh = 0.5, padj_thresh = 0.05,
                   label_genes = highlight_genes,
                   xlim = c(-3, 3)),
       width = 2.8, height = 3.0, device = cairo_pdf)

# ── Fig 3: Heatmap of 22 XBP1-s vs XBP1-es DEGs ─────────────────────────────
sig_22 <- read_csv("data/DEG_results/XBP1s_vs_XBP1es_significant.csv")
deg_genes <- sig_22$gene_name

# Get rlog expression for these genes
rld_mat <- assay(rld)
rownames(rld_mat) <- rowData(rld)$gene_name  # requires gene_name in rowData
# If gene_name not in rowData, merge via gene_id
# rld_mat_df <- as.data.frame(rld_mat) %>% rownames_to_column("gene_id") %>%
#   left_join(gene_map, by = "gene_id")

mat_22 <- rld_mat[rownames(rld_mat) %in% deg_genes, ]

# Z-score per gene
mat_z <- t(scale(t(mat_22)))

# Column annotation
col_ann <- HeatmapAnnotation(
  Condition = meta$condition,
  col = list(Condition = condition_colors),
  annotation_name_gp = gpar(fontsize = 5),
  simple_anno_size = unit(0.25, "cm")
)

ht <- Heatmap(
  mat_z,
  name            = "Z-score",
  col             = colorRamp2(c(-2, 0, 2), c("#2A7DB5", "white", "#E05C2A")),
  top_annotation  = col_ann,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 5),
  row_names_gp    = gpar(fontsize = 5),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 4.5)
  ),
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  row_dend_width  = unit(0.5, "cm"),
  width           = unit(4, "cm"),
  height          = unit(6, "cm")
)

pdf("figures/fig3_heatmap_XBP1s_vs_XBP1es_DEGs.pdf",
    width = 3.0, height = 4.0)
draw(ht)
dev.off()

cat("Figures 1-3 saved to figures/\n")
