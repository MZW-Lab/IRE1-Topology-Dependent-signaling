# =============================================================================
# Script 04: RIDD Target Time-Series (Fig. 3C)
# =============================================================================
# Manuscript: "A topology-dependent IRE1alpha signaling axis generates a cryptic
#              XBP1 isoform present in human cancers"
# Journal: Nature Chemical Biology
#
# Purpose:
#   Plot log2 fold change over time for three individually validated canonical
#   RIDD target transcripts (BLOC1S1, BCAM, CD59) in ER- and Cyto-opto-IRE1alpha
#   cells. Values are normalized to total IRE1alpha protein levels to correct for
#   the ~4-fold expression difference between constructs. A 3-point rolling mean
#   is applied to reduce per-timepoint noise inherent to single-sample time-series
#   data; raw values are shown as faint points behind the smoothed line.
#
# Input:
#   data/vst_combined.csv   -- VST-normalized expression (Script 01 output)
#
# Output:
#   figures/Fig3C_RIDD_targets.pdf
#
# IRE1alpha normalization (from manuscript Methods):
#   Western blot quantification of total IRE1alpha normalized to beta-actin:
#     ER fraction mean:   3,144.48 AU
#     Cyto fraction mean: 16,468.07 AU
#     Overall mean:       9,806.28 AU
#   Normalization factors:
#     ER   = 3,144.48 / 9,806.28 = 0.321
#     Cyto = 16,468.07 / 9,806.28 = 1.679
#   Normalized log2FC = observed log2FC / normalization_factor
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# -- 0. Paths ------------------------------------------------------------------
data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# -- 1. IRE1alpha normalization factors (from manuscript Methods) --------------
norm_factor_ER   <- 0.321
norm_factor_Cyto <- 1.679

# -- 2. RIDD target genes ------------------------------------------------------
ridd_genes <- c("BLOC1S1", "BCAM", "CD59")

# -- 3. Load VST matrix --------------------------------------------------------
message("Loading VST-normalized expression matrix...")
vst_df <- read.csv(file.path(data_dir, "vst_combined.csv"), check.names = FALSE)
rownames(vst_df) <- vst_df$gene; vst_df$gene <- NULL
vst_mat <- as.matrix(vst_df)

missing <- setdiff(ridd_genes, rownames(vst_mat))
if (length(missing) > 0) stop(sprintf("Missing genes: %s", paste(missing, collapse = ", ")))

# -- 4. Compute log2FC relative to t=0 ----------------------------------------
timepoints <- seq(0, 345, by = 15)

cyto_cols <- grep("^C", colnames(vst_mat), value = TRUE)
er_cols   <- grep("^E", colnames(vst_mat), value = TRUE)

cyto_lfc <- (vst_mat[ridd_genes, cyto_cols] - vst_mat[ridd_genes, "C0"]) / norm_factor_Cyto
er_lfc   <- (vst_mat[ridd_genes, er_cols]   - vst_mat[ridd_genes, "E0"]) / norm_factor_ER

# -- 5. 3-point rolling mean (uniform boxcar filter) --------------------------
# Applied per gene per condition to reduce per-timepoint noise.
# Edge timepoints use a smaller window (2-point average at t=0 and t=345).
rolling_mean3 <- function(x) {
  n   <- length(x)
  out <- numeric(n)
  for (i in seq_len(n)) {
    idx    <- max(1L, i - 1L):min(n, i + 1L)
    out[i] <- mean(x[idx])
  }
  out
}

# -- 6. Reshape to long format and apply smoothing ----------------------------
make_long <- function(mat, condition) {
  df       <- as.data.frame(mat)
  df$gene  <- rownames(df)
  colnames(df)[seq_along(timepoints)] <- as.character(timepoints)
  df_long  <- pivot_longer(df, cols = as.character(timepoints),
                            names_to = "time_min", values_to = "log2FC_norm")
  df_long$time_min  <- as.numeric(df_long$time_min)
  df_long$condition <- condition
  df_long <- df_long %>%
    group_by(gene, condition) %>%
    arrange(time_min) %>%
    mutate(log2FC_smooth = rolling_mean3(log2FC_norm)) %>%
    ungroup()
  df_long
}

plot_df <- bind_rows(
  make_long(cyto_lfc, "Cyto-opto-IRE1a"),
  make_long(er_lfc,   "ER-opto-IRE1a")
)

# -- 7. Plot -------------------------------------------------------------------
condition_colors <- c("Cyto-opto-IRE1a" = "#85ab8e", "ER-opto-IRE1a" = "#f6992d")

base_theme <- theme_classic(base_size = 7) +
  theme(
    text             = element_text(family = "sans", size = 7),
    axis.title       = element_text(size = 7),
    axis.text        = element_text(size = 6),
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 6),
    legend.key.size  = unit(0.3, "cm"),
    strip.text       = element_text(size = 7, face = "italic"),
    strip.background = element_blank(),
    plot.title       = element_text(size = 7, face = "bold"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

p <- ggplot(plot_df, aes(x = time_min, color = condition)) +
  # Raw values as faint points
  geom_point(aes(y = log2FC_norm), size = 0.5, alpha = 0.3, shape = 16) +
  # 3-point rolling mean as solid line
  geom_line(aes(y = log2FC_smooth), linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~ gene, ncol = 3, scales = "free_y") +
  scale_color_manual(values = condition_colors, name = "Construct") +
  scale_x_continuous(breaks = seq(0, 345, by = 60)) +
  labs(
    title = "Fig. 3C: RIDD target transcripts (3-point rolling mean)",
    x     = "Time (min)",
    y     = "log2FC (normalized to IRE1a expression)"
  ) +
  base_theme

# -- 8. Save -------------------------------------------------------------------
tmp_path <- file.path(tempdir(), "Fig3C_RIDD_targets.pdf")
pdf(tmp_path, width = 5.5, height = 2.2, useDingbats = FALSE)
print(p)
dev.off()

out_path <- file.path(fig_dir, "Fig3C_RIDD_targets.pdf")
file.copy(tmp_path, out_path, overwrite = TRUE)
message(sprintf("Saved: %s", out_path))
message("Script 04 complete.")
