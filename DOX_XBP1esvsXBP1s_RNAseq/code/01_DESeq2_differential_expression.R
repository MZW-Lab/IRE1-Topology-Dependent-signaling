# =============================================================================
# Script 01: DESeq2 Differential Expression Analysis (Run B)
# =============================================================================
# Description:
#   Performs differential expression analysis comparing XBP1-s and XBP1-es
#   isoforms against wild-type (WT) controls, and XBP1-s vs XBP1-es directly.
#
#   Run B uses a pseudo-replicate for XBP1-s replicate 2 (a failed sequencing
#   run), constructed as floor((rep1 + rep3) / 2). See README for details.
#
# Input:
#   data/raw/expression_matrix_raw_counts_cpm.tsv  — raw count matrix (all 9 samples)
#   data/raw/sample_metadata_runB.csv              — sample metadata
#
# Output:
#   data/DEG_results/XBP1s_vs_WT_full.csv
#   data/DEG_results/XBP1s_vs_WT_significant.csv
#   data/DEG_results/XBP1es_vs_WT_full.csv
#   data/DEG_results/XBP1es_vs_WT_significant.csv
#   data/DEG_results/XBP1s_vs_XBP1es_full.csv
#   data/DEG_results/XBP1s_vs_XBP1es_significant.csv
#
# Thresholds:
#   XBP1-s vs WT / XBP1-es vs WT : |log2FC| >= 1.0, padj <= 0.05
#   XBP1-s vs XBP1-es            : |log2FC| >= 0.5, padj <= 0.05
#
# Software: R >= 4.2, DESeq2 >= 1.38, apeglm >= 1.20, ashr >= 2.2
# =============================================================================

library(DESeq2)
library(apeglm)
library(ashr)
library(dplyr)
library(readr)

set.seed(42)

# ── 1. Load raw count matrix ──────────────────────────────────────────────────
mat_raw <- read_tsv("data/raw/expression_matrix_raw_counts_cpm.tsv")

# Extract integer count columns (columns ending in _count)
count_cols <- grep("_count$", colnames(mat_raw), value = TRUE)
counts_raw <- mat_raw[, c("gene_id", "gene_name", count_cols)]

# Rename columns: strip _count suffix and map sample order
# Sample order in matrix: _1 = WT_DOX_1, _2 = WT_DOX_2, _3 = WT_DOX_3,
#                         _4 = XBP1s_DOX_1, _5 = XBP1s_DOX_2 (FAILED),
#                         _6 = XBP1s_DOX_3, _7 = XBP1es_DOX_1,
#                         _8 = XBP1es_DOX_2, _9 = XBP1es_DOX_3
sample_names <- c("WT_DOX_1", "WT_DOX_2", "WT_DOX_3",
                  "XBP1s_DOX_1", "XBP1s_DOX_2_FAILED", "XBP1s_DOX_3",
                  "XBP1es_DOX_1", "XBP1es_DOX_2", "XBP1es_DOX_3")
colnames(counts_raw)[3:ncol(counts_raw)] <- sample_names

# Round to integers (counts may be fractional from upstream pipeline)
count_mat <- round(as.matrix(counts_raw[, sample_names]))
rownames(count_mat) <- counts_raw$gene_id

# ── 2. Construct Run B pseudo-replicate for XBP1s_DOX_2 ──────────────────────
# XBP1s_DOX_2 was a failed sequencing run (~19K total counts vs 5-9M for others).
# A re-sequenced library was provided but remained a transcriptomic outlier
# (PC1 = 25.1 vs 0.0-1.6 for group mates; Pearson r = 0.935-0.943 vs 0.988-0.999).
# Run B replaces it with a pseudo-replicate: floor((rep1 + rep3) / 2).
pseudo_rep2 <- floor((count_mat[, "XBP1s_DOX_1"] + count_mat[, "XBP1s_DOX_3"]) / 2)
count_mat_B <- count_mat[, c("WT_DOX_1", "WT_DOX_2", "WT_DOX_3",
                              "XBP1s_DOX_1", "XBP1s_DOX_3",
                              "XBP1es_DOX_1", "XBP1es_DOX_2", "XBP1es_DOX_3")]
count_mat_B <- cbind(count_mat_B[, 1:4], XBP1s_DOX_2_pseudo = pseudo_rep2,
                     count_mat_B[, 5:8])

# Reorder to match metadata
count_mat_B <- count_mat_B[, c("WT_DOX_1", "WT_DOX_2", "WT_DOX_3",
                                "XBP1s_DOX_1", "XBP1s_DOX_2_pseudo", "XBP1s_DOX_3",
                                "XBP1es_DOX_1", "XBP1es_DOX_2", "XBP1es_DOX_3")]

# ── 3. Sample metadata ────────────────────────────────────────────────────────
meta <- data.frame(
  sample    = colnames(count_mat_B),
  condition = c("WT_DOX", "WT_DOX", "WT_DOX",
                "XBP1s_DOX", "XBP1s_DOX", "XBP1s_DOX",
                "XBP1es_DOX", "XBP1es_DOX", "XBP1es_DOX"),
  row.names = colnames(count_mat_B)
)
meta$condition <- factor(meta$condition,
                         levels = c("WT_DOX", "XBP1s_DOX", "XBP1es_DOX"))

# ── 4. Pre-filtering ──────────────────────────────────────────────────────────
# Keep genes with >= 10 counts in >= 3 samples
keep <- rowSums(count_mat_B >= 10) >= 3
count_mat_B <- count_mat_B[keep, ]
cat(sprintf("Genes after pre-filtering: %d\n", nrow(count_mat_B)))

# ── 5. DESeq2 object ─────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = count_mat_B,
  colData   = meta,
  design    = ~ condition
)
dds <- DESeq(dds)

# ── 6. Helper: extract + shrink results ──────────────────────────────────────
get_results <- function(dds, contrast, coef = NULL,
                        shrink_type = "apeglm",
                        lfc_thresh, padj_thresh,
                        gene_map) {

  if (shrink_type == "apeglm") {
    res_raw <- results(dds, name = coef)
    res     <- lfcShrink(dds, coef = coef, type = "apeglm", res = res_raw)
  } else {
    res_raw <- results(dds, contrast = contrast)
    res     <- lfcShrink(dds, contrast = contrast, type = "ashr", res = res_raw)
  }

  df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    left_join(gene_map, by = "gene_id") %>%
    arrange(padj)

  sig <- df %>%
    filter(!is.na(padj), padj <= padj_thresh,
           abs(log2FoldChange) >= lfc_thresh)

  list(full = df, significant = sig)
}

# Gene ID -> name map
gene_map <- counts_raw %>% select(gene_id, gene_name) %>% distinct()

# ── 7. XBP1-s vs WT (apeglm, |LFC| >= 1, padj <= 0.05) ──────────────────────
res_s_wt <- get_results(
  dds, coef = "condition_XBP1s_DOX_vs_WT_DOX",
  shrink_type = "apeglm", lfc_thresh = 1.0, padj_thresh = 0.05,
  gene_map = gene_map
)
cat(sprintf("XBP1-s vs WT: %d significant DEGs\n", nrow(res_s_wt$significant)))

# ── 8. XBP1-es vs WT (apeglm, |LFC| >= 1, padj <= 0.05) ─────────────────────
res_es_wt <- get_results(
  dds, coef = "condition_XBP1es_DOX_vs_WT_DOX",
  shrink_type = "apeglm", lfc_thresh = 1.0, padj_thresh = 0.05,
  gene_map = gene_map
)
cat(sprintf("XBP1-es vs WT: %d significant DEGs\n", nrow(res_es_wt$significant)))

# ── 9. XBP1-s vs XBP1-es (ashr, |LFC| >= 0.5, padj <= 0.05) ─────────────────
res_s_es <- get_results(
  dds,
  contrast    = c("condition", "XBP1s_DOX", "XBP1es_DOX"),
  shrink_type = "ashr",
  lfc_thresh  = 0.5, padj_thresh = 0.05,
  gene_map    = gene_map
)
cat(sprintf("XBP1-s vs XBP1-es: %d significant DEGs\n", nrow(res_s_es$significant)))

# ── 10. Save outputs ──────────────────────────────────────────────────────────
dir.create("data/DEG_results", recursive = TRUE, showWarnings = FALSE)

write_csv(res_s_wt$full,        "data/DEG_results/XBP1s_vs_WT_full.csv")
write_csv(res_s_wt$significant, "data/DEG_results/XBP1s_vs_WT_significant.csv")
write_csv(res_es_wt$full,       "data/DEG_results/XBP1es_vs_WT_full.csv")
write_csv(res_es_wt$significant,"data/DEG_results/XBP1es_vs_WT_significant.csv")
write_csv(res_s_es$full,        "data/DEG_results/XBP1s_vs_XBP1es_full.csv")
write_csv(res_s_es$significant, "data/DEG_results/XBP1s_vs_XBP1es_significant.csv")

cat("DESeq2 analysis complete. Results saved to data/DEG_results/\n")
