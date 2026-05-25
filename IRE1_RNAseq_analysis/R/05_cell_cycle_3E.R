# =============================================================================
# Script 05: Cell-Cycle Gene Set Mean Log2FC (Fig. 3E)
# =============================================================================
# Manuscript: "A topology-dependent IRE1α signaling axis generates a cryptic
#              XBP1 isoform present in human cancers"
# Journal: Nature Chemical Biology
#
# Purpose:
#   Plot mean log2 fold change (relative to t=0) over time for genes
#   associated with cell-cycle phases (G1/S, S, G2/M, M/G1) and cell-cycle
#   arrest in ER- and Cyto-opto-IRE1α cells. Shaded area = SEM across genes.
#
# Input:
#   data/vst_combined.csv      — VST-normalized expression (Script 01 output)
#   data/PanelA_gene_sets.csv  — Cell-cycle gene sets with status annotation
#
# Output:
#   figures/Fig3E_cell_cycle.pdf
#
# Gene set source:
#   PanelA_gene_sets.csv contains 5 gene sets: G1/S, S, G2/M, M/G1, Arrest.
#   Only genes with status == "used" are included (2 genes excluded as
#   low-expression at T0).
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── 0. Paths ──────────────────────────────────────────────────────────────────
data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# ── 1. Load data ──────────────────────────────────────────────────────────────
message("Loading VST matrix and gene sets...")
vst_df <- read.csv(file.path(data_dir, "vst_combined.csv"), check.names = FALSE)
rownames(vst_df) <- vst_df$gene; vst_df$gene <- NULL
vst_mat <- as.matrix(vst_df)

gene_sets_df <- read.csv(file.path(data_dir, "PanelA_gene_sets.csv"),
                         stringsAsFactors = FALSE)

# Keep only "used" genes
gene_sets_df <- gene_sets_df %>% filter(status == "used")
message(sprintf("Gene sets loaded: %d genes across %d sets",
                nrow(gene_sets_df), length(unique(gene_sets_df$gene_set))))

# ── 2. Compute log2FC relative to t=0 ────────────────────────────────────────
timepoints <- seq(0, 345, by = 15)

cyto_cols <- grep("^C", colnames(vst_mat), value = TRUE)
er_cols   <- grep("^E", colnames(vst_mat), value = TRUE)

# ── 3. Compute mean ± SEM per gene set per timepoint ─────────────────────────
compute_set_stats <- function(genes, vst_mat, timepoints, condition_cols, t0_col) {
  # Filter to genes present in VST matrix
  genes_present <- intersect(genes, rownames(vst_mat))
  if (length(genes_present) == 0) return(NULL)

  lfc_mat <- vst_mat[genes_present, condition_cols] -
             vst_mat[genes_present, t0_col]

  # Mean and SEM across genes at each timepoint
  data.frame(
    time_min = timepoints,
    mean_lfc = colMeans(lfc_mat, na.rm = TRUE),
    sem_lfc  = apply(lfc_mat, 2, function(x) sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x))))
  )
}

# Gene set order for faceting
set_order <- c("G1/S", "S", "G2/M", "M/G1", "Arrest")

all_stats <- list()
for (gs in set_order) {
  genes_gs <- gene_sets_df %>% filter(gene_set == gs) %>% pull(gene) %>% unique()

  cyto_stats <- compute_set_stats(genes_gs, vst_mat, timepoints, cyto_cols, "C0")
  er_stats   <- compute_set_stats(genes_gs, vst_mat, timepoints, er_cols,   "E0")

  if (!is.null(cyto_stats)) {
    cyto_stats$condition <- "Cyto-opto-IRE1a"
    cyto_stats$gene_set  <- gs
    cyto_stats$n_genes   <- length(intersect(genes_gs, rownames(vst_mat)))
  }
  if (!is.null(er_stats)) {
    er_stats$condition <- "ER-opto-IRE1a"
    er_stats$gene_set  <- gs
    er_stats$n_genes   <- length(intersect(genes_gs, rownames(vst_mat)))
  }
  all_stats[[gs]] <- bind_rows(cyto_stats, er_stats)
}

plot_df <- bind_rows(all_stats)
plot_df$gene_set <- factor(plot_df$gene_set, levels = set_order)

# ── 4. Plot ───────────────────────────────────────────────────────────────────
condition_colors <- c("Cyto-opto-IRE1a" = "#85ab8e", "ER-opto-IRE1a" = "#f6992d")
condition_labels <- c("Cyto-opto-IRE1a" = "Cyto-opto-IRE1\u03b1",
                      "ER-opto-IRE1a"   = "ER-opto-IRE1\u03b1")

base_theme <- theme_classic(base_size = 7) +
  theme(
    text             = element_text(family = "sans", size = 7),
    axis.title       = element_text(size = 7),
    axis.text        = element_text(size = 6),
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.3, "cm"),
    strip.text       = element_text(size = 7, face = "bold"),
    strip.background = element_blank(),
    plot.title       = element_text(size = 7, face = "bold"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

p <- ggplot(plot_df, aes(x = time_min, y = mean_lfc, color = condition, fill = condition)) +
  geom_ribbon(aes(ymin = mean_lfc - sem_lfc, ymax = mean_lfc + sem_lfc),
              alpha = 0.2, color = NA) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~ gene_set, ncol = 5, scales = "free_y") +
  scale_color_manual(values = condition_colors,
                     labels = condition_labels,
                     name   = "Construct") +
  scale_fill_manual(values  = condition_colors,
                    labels  = condition_labels,
                    name    = "Construct") +
  scale_x_continuous(breaks = seq(0, 345, by = 120),
                     labels = seq(0, 345, by = 120)) +
  labs(
    title = "Fig. 3E: Cell-cycle gene sets",
    x     = "Time (min)",
    y     = expression("Mean log"[2]*"FC (relative to t=0)")
  ) +
  base_theme

# ── 5. Save ───────────────────────────────────────────────────────────────────
tmp_path <- file.path(tempdir(), "Fig3E_cell_cycle.pdf")
pdf(tmp_path, width = 7, height = 2.2, useDingbats = FALSE)
print(p)
dev.off()

out_path <- file.path(fig_dir, "Fig3E_cell_cycle.pdf")
file.copy(tmp_path, out_path, overwrite = TRUE)
message(sprintf("Saved: %s", out_path))
message("Script 05 complete.")
