# =============================================================================
# Script 02: UMAP of All 48 RNA-seq Samples (Fig. S4A–B)
# =============================================================================
# Manuscript: "A topology-dependent IRE1α signaling axis generates a cryptic
#              XBP1 isoform present in human cancers"
# Journal: Nature Chemical Biology
#
# Purpose:
#   Project all 48 RNA-seq samples (24 ER + 24 Cyto, 0–345 min) into 2D UMAP
#   space using VST-normalized expression of the top 500 most variable genes.
#   Two panels:
#     Fig. S4A — samples colored by construct identity (ER vs. Cyto)
#     Fig. S4B — samples colored by time point (continuous gradient)
#
# Input:
#   data/vst_combined.csv   — VST-normalized expression (output of Script 01)
#
# Output:
#   figures/FigS4AB_UMAP.pdf
#
# Dependencies: Script 01 must be run first to generate vst_combined.csv
# =============================================================================

suppressPackageStartupMessages({
  library(uwot)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)

# ── 0. Paths ──────────────────────────────────────────────────────────────────
data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# ── 1. Load VST matrix ────────────────────────────────────────────────────────
message("Loading VST-normalized expression matrix...")
vst_df <- read.csv(file.path(data_dir, "vst_combined.csv"), check.names = FALSE)
rownames(vst_df) <- vst_df$gene
vst_df$gene <- NULL
vst_mat <- as.matrix(vst_df)   # genes × 48 samples
message(sprintf("VST matrix: %d genes x %d samples", nrow(vst_mat), ncol(vst_mat)))

# ── 2. Select top 500 most variable genes ─────────────────────────────────────
gene_vars <- apply(vst_mat, 1, var)
top500    <- names(sort(gene_vars, decreasing = TRUE))[1:500]
vst_t     <- t(vst_mat[top500, ])   # samples × genes (required by uwot)

# ── 3. Run UMAP ───────────────────────────────────────────────────────────────
message("Running UMAP (top 500 variable genes, n_neighbors=15)...")
umap_coords <- umap(vst_t, n_neighbors = 15, min_dist = 0.3,
                    metric = "euclidean", n_epochs = 500, verbose = FALSE)

# ── 4. Build metadata data frame ──────────────────────────────────────────────
timepoints <- seq(0, 345, by = 15)

# Use ASCII-safe labels; Greek alpha rendered via scale label expressions
meta <- data.frame(
  sample    = colnames(vst_mat),
  UMAP1     = umap_coords[, 1],
  UMAP2     = umap_coords[, 2],
  condition = c(rep("Cyto", 24), rep("ER", 24)),
  time_min  = rep(timepoints, 2),
  stringsAsFactors = FALSE
)

# ── 5. Styling ────────────────────────────────────────────────────────────────
# Per user preference: Cyto = #85ab8e, ER = #f6992d
condition_colors <- c("Cyto" = "#85ab8e", "ER" = "#f6992d")

# Use plotmath expressions for Greek alpha in legend labels
condition_labels <- c(
  "Cyto" = expression("Cyto-opto-IRE1" * alpha),
  "ER"   = expression("ER-opto-IRE1" * alpha)
)

base_theme <- theme_classic(base_size = 7) +
  theme(
    text            = element_text(family = "sans", size = 7),
    axis.title      = element_text(size = 7),
    axis.text       = element_text(size = 6),
    legend.title    = element_text(size = 7),
    legend.text     = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    plot.title      = element_text(size = 7, face = "bold"),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# ── 6. Panel A: colored by construct identity ─────────────────────────────────
pA <- ggplot(meta, aes(x = UMAP1, y = UMAP2, color = condition)) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_color_manual(values  = condition_colors,
                     labels  = condition_labels,
                     name    = "Construct") +
  labs(title = "Fig. S4A", x = "UMAP 1", y = "UMAP 2") +
  base_theme

# ── 7. Panel B: colored by time point ────────────────────────────────────────
pB <- ggplot(meta, aes(x = UMAP1, y = UMAP2, color = time_min)) +
  geom_point(size = 1.8, alpha = 0.9) +
  scale_color_viridis_c(name = "Time (min)", option = "plasma") +
  labs(title = "Fig. S4B", x = "UMAP 1", y = "UMAP 2") +
  base_theme

# ── 8. Save combined figure ───────────────────────────────────────────────────
# Write to local workspace first (PDF is a random-access format), then copy
tmp_path <- file.path(tempdir(), "FigS4AB_UMAP.pdf")
pdf(tmp_path, width = 6, height = 2.8, useDingbats = FALSE)
print(pA + pB + plot_layout(ncol = 2))
dev.off()

out_path <- file.path(fig_dir, "FigS4AB_UMAP.pdf")
file.copy(tmp_path, out_path, overwrite = TRUE)
message(sprintf("Saved: %s", out_path))
message("Script 02 complete.")
