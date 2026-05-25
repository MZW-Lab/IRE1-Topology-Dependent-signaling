# =============================================================================
# Script 03: Genome-wide Log2FC Heatmap with K-means Clusters (Fig. 3B)
# =============================================================================
# Manuscript: "A topology-dependent IRE1α signaling axis generates a cryptic
#              XBP1 isoform present in human cancers"
# Journal: Nature Chemical Biology
#
# Purpose:
#   Generate a genome-wide heatmap of log2 fold changes (relative to t=0)
#   for all expressed genes across the 345-min time course, split by the
#   pre-computed k-means cluster assignments (k=6). Cyto and ER are shown
#   side by side. Cluster annotations reflect GO Biological Process terms
#   as reported in the manuscript.
#
# Input:
#   data/vst_combined.csv        — VST-normalized expression (Script 01 output)
#   data/cluster_assignments.csv — Pre-computed k-means cluster labels (k=6)
#
# Output:
#   figures/Fig3B_heatmap.pdf
#
# Cluster identity (verified by GO enrichment; see manuscript Fig. 3B):
#   Cluster 0 — Cytoskeletal organization / cell adhesion / small GTPase
#               signaling; selectively downregulated in Cyto-opto-IRE1α
#               (manuscript "C1")
#   Cluster 1 — Synaptic/secretion-related; modest change
#   Cluster 2 — Extracellular matrix / ER-luminal / secreted transcripts;
#               canonical RIDD targets decayed in both conditions
#               (manuscript "C3")
#   Cluster 3 — ER stress response / protein folding; upregulated in both
#               (manuscript "C0/C4" ER stress clusters)
#   Cluster 4 — Gland/liver development; modest change
#   Cluster 5 — Ribosome biogenesis / RNA processing; upregulated in Cyto
#               (manuscript "C2" spliceosome/RNA processing cluster)
# =============================================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
})

set.seed(42)

# ── 0. Paths ──────────────────────────────────────────────────────────────────
# Run from repo root: Rscript R/03_heatmap_3B.R
data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
message("Loading VST matrix and cluster assignments...")
vst_df <- read.csv(file.path(data_dir, "vst_combined.csv"), check.names = FALSE)
rownames(vst_df) <- vst_df$gene; vst_df$gene <- NULL
vst_mat <- as.matrix(vst_df)

clusters_df <- read.csv(file.path(data_dir, "cluster_assignments.csv"),
                        stringsAsFactors = FALSE)
rownames(clusters_df) <- clusters_df$gene

# ── 2. Compute log2FC relative to t=0 ────────────────────────────────────────
message("Computing log2FC relative to t=0...")
cyto_cols <- grep("^C", colnames(vst_mat), value = TRUE)
er_cols   <- grep("^E", colnames(vst_mat), value = TRUE)

cyto_lfc <- vst_mat[, cyto_cols] - vst_mat[, "C0"]
er_lfc   <- vst_mat[, er_cols]   - vst_mat[, "E0"]

timepoints <- seq(0, 345, by = 15)
colnames(cyto_lfc) <- paste0(timepoints, "min")
colnames(er_lfc)   <- paste0(timepoints, "min")

# ── 3. Filter to clustered genes ──────────────────────────────────────────────
common_genes  <- intersect(rownames(clusters_df), rownames(cyto_lfc))
message(sprintf("Genes with cluster assignments: %d", length(common_genes)))

cyto_lfc_filt <- cyto_lfc[common_genes, ]
er_lfc_filt   <- er_lfc[common_genes, ]
cluster_vec   <- clusters_df[common_genes, "cluster"]

# ── 4. Cluster annotation labels (GO-verified) ────────────────────────────────
cluster_labels <- c(
  "0" = "C0: Cytoskeletal/adhesion/GTPase (Cyto-selective decay)",
  "1" = "C1: Synaptic/secretion-related",
  "2" = "C2: ECM/ER-luminal/secreted (RIDD targets, both)",
  "3" = "C3: ER stress/protein folding (upregulated, both)",
  "4" = "C4: Gland/liver development",
  "5" = "C5: Ribosome biogenesis/RNA processing (Cyto-selective up)"
)

cluster_colors <- c(
  "0" = "#d6604d",   # red    — Cyto-selective decay (manuscript C1)
  "1" = "#b2abd2",   # purple — modest
  "2" = "#92c5de",   # light blue — RIDD targets (manuscript C3)
  "3" = "#4393c3",   # blue   — ER stress up (both)
  "4" = "#d9d9d9",   # grey   — modest
  "5" = "#f4a582"    # salmon — Cyto-selective up (manuscript C2)
)

# ── 5. Build combined matrix: Cyto | gap | ER ─────────────────────────────────
gap_col      <- matrix(NA, nrow = length(common_genes), ncol = 1,
                       dimnames = list(common_genes, "gap"))
combined_lfc <- cbind(cyto_lfc_filt, gap_col, er_lfc_filt)

# ── 6. Color scale (symmetric, capped at ±2) ──────────────────────────────────
lfc_cap <- 2
col_fun <- colorRamp2(
  c(-lfc_cap, -1, 0, 1, lfc_cap),
  c("#2166ac", "#92c5de", "white", "#f4a582", "#d6604d")
)

# ── 7. Row annotation ─────────────────────────────────────────────────────────
row_ann <- rowAnnotation(
  Cluster = factor(as.character(cluster_vec), levels = as.character(0:5)),
  col     = list(Cluster = cluster_colors),
  annotation_name_gp = gpar(fontsize = 6),
  width = unit(0.3, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    Cluster = list(
      title      = "Cluster",
      title_gp   = gpar(fontsize = 6),
      labels     = cluster_labels[as.character(0:5)],
      labels_gp  = gpar(fontsize = 5)
    )
  )
)

# ── 8. Column annotation (condition bar) ──────────────────────────────────────
col_condition <- c(rep("Cyto", 24), "gap", rep("ER", 24))
col_ann <- HeatmapAnnotation(
  Condition = col_condition,
  col = list(Condition = c(
    "Cyto" = "#85ab8e",
    "ER"   = "#f6992d",
    "gap"  = "white"
  )),
  annotation_name_gp = gpar(fontsize = 6),
  height = unit(0.3, "cm"),
  show_legend = TRUE,
  annotation_legend_param = list(
    Condition = list(
      title     = "Condition",
      labels    = c("Cyto" = "Cyto-opto-IRE1a", "ER" = "ER-opto-IRE1a", "gap" = ""),
      labels_gp = gpar(fontsize = 5),
      title_gp  = gpar(fontsize = 6)
    )
  )
)

# ── 9. Column labels (every 60 min per condition) ─────────────────────────────
show_times  <- seq(0, 345, by = 60)
cyto_labels <- ifelse(timepoints %in% show_times, paste0(timepoints, "m"), "")
er_labels   <- ifelse(timepoints %in% show_times, paste0(timepoints, "m"), "")
col_labels  <- c(cyto_labels, "", er_labels)   # length 49

# ── 10. Draw heatmap ──────────────────────────────────────────────────────────
message("Drawing heatmap...")
ht <- Heatmap(
  combined_lfc,
  name                 = "log2FC",
  col                  = col_fun,
  na_col               = "white",
  cluster_rows         = FALSE,
  cluster_columns      = FALSE,
  split                = factor(as.character(cluster_vec), levels = as.character(0:5)),
  row_title_gp         = gpar(fontsize = 6, fontface = "bold"),
  row_title_rot        = 0,
  show_row_names       = FALSE,
  column_labels        = col_labels,
  column_names_gp      = gpar(fontsize = 5),
  column_names_rot     = 45,
  top_annotation       = col_ann,
  left_annotation      = row_ann,
  heatmap_legend_param = list(
    title         = "log2FC",
    title_gp      = gpar(fontsize = 6),
    labels_gp     = gpar(fontsize = 5),
    legend_height = unit(2, "cm"),
    grid_width    = unit(0.3, "cm")
  ),
  use_raster     = TRUE,
  raster_quality = 5
)

# ── 11. Save ──────────────────────────────────────────────────────────────────
tmp_path <- file.path(tempdir(), "Fig3B_heatmap.pdf")
pdf(tmp_path, width = 7, height = 8, useDingbats = FALSE)
draw(ht, merge_legend = TRUE)
dev.off()

out_path <- file.path(fig_dir, "Fig3B_heatmap.pdf")
file.copy(tmp_path, out_path, overwrite = TRUE)
message(sprintf("Saved: %s", out_path))
message("Script 03 complete.")
