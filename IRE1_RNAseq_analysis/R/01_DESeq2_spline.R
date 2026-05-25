# =============================================================================
# Script 01: DESeq2 Differential Expression with Natural Cubic Spline Design
# =============================================================================
# Manuscript: "A topology-dependent IRE1α signaling axis generates a cryptic
#              XBP1 isoform present in human cancers"
# Journal: Nature Chemical Biology
#
# Purpose:
#   Identify differentially expressed genes (DEGs) across the 345-min time
#   course for ER-opto-IRE1α and Cyto-opto-IRE1α separately, using DESeq2
#   with a natural cubic spline design and likelihood ratio test (LRT).
#
# Input:
#   data/rawcounts_Cyto_symbols.csv  — 28,395 genes × 24 timepoints
#   data/rawcounts_ER_symbols.csv    — 28,395 genes × 24 timepoints
#
# Output:
#   data/DESeq2_Cyto_results.csv     — Full DESeq2 results for Cyto
#   data/DESeq2_ER_results.csv       — Full DESeq2 results for ER
#   data/DESeq2_Cyto_DEGs.csv        — Significant DEGs (padj < 1e-3) for Cyto
#   data/DESeq2_ER_DEGs.csv          — Significant DEGs (padj < 1e-3) for ER
#   data/vst_combined.csv            — VST-normalized expression (all 48 samples)
#
# Method:
#   - Natural cubic spline with df=4 models the time-course shape
#   - LRT compares full model [~ ns(time, df=4)] vs. reduced [~ 1]
#   - Run separately per condition (each condition has 24 timepoints, 1 sample each)
#   - Significance threshold: padj < 1e-3 (as reported in manuscript)
#
# Note on pseudoreplication:
#   This time-series design has one sample per timepoint (no biological
#   replicates per timepoint). DESeq2 with spline LRT is the standard approach
#   for such designs; dispersion is estimated across the time-course.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(splines)
  library(dplyr)
  library(tibble)
})

set.seed(42)

# ── 0. Paths ──────────────────────────────────────────────────────────────────
# When running from the repo root: Rscript R/01_DESeq2_spline.R
data_dir <- "data"

# ── 1. Load raw counts (drop rows with NA gene names) ─────────────────────────
message("Loading raw count matrices...")
cyto_df <- read.csv(file.path(data_dir, "rawcounts_Cyto_symbols.csv"), check.names = FALSE)
er_df   <- read.csv(file.path(data_dir, "rawcounts_ER_symbols.csv"),   check.names = FALSE)

# Remove rows with NA gene symbol (617 unnamed genes in each matrix)
cyto_df <- cyto_df[!is.na(cyto_df$symbol), ]
er_df   <- er_df[!is.na(er_df$symbol), ]

rownames(cyto_df) <- cyto_df$symbol; cyto_df$symbol <- NULL
rownames(er_df)   <- er_df$symbol;   er_df$symbol   <- NULL

cyto_raw <- as.matrix(cyto_df)
er_raw   <- as.matrix(er_df)

timepoints <- seq(0, 345, by = 15)   # 0, 15, 30, ..., 345 min (24 points)
stopifnot(ncol(cyto_raw) == 24, ncol(er_raw) == 24)

message(sprintf("Cyto: %d genes × %d timepoints", nrow(cyto_raw), ncol(cyto_raw)))
message(sprintf("ER:   %d genes × %d timepoints", nrow(er_raw),   ncol(er_raw)))

# ── 2. Helper: run DESeq2 spline LRT for one condition ───────────────────────
run_deseq2_spline <- function(counts_mat, timepoints, df = 4, padj_thresh = 1e-3) {
  # Build colData with spline basis columns
  col_data <- data.frame(
    sample = colnames(counts_mat),
    time   = timepoints,
    row.names = colnames(counts_mat)
  )
  spline_basis <- ns(timepoints, df = df)
  colnames(spline_basis) <- paste0("ns", seq_len(df))
  col_data <- cbind(col_data, spline_basis)

  full_formula    <- as.formula(paste("~", paste(colnames(spline_basis), collapse = " + ")))
  reduced_formula <- as.formula("~ 1")

  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_mat),
    colData   = col_data,
    design    = full_formula
  )

  # Pre-filter: keep genes with >= 10 counts in at least 1 sample
  keep <- rowSums(counts(dds) >= 10) >= 1
  dds  <- dds[keep, ]
  message(sprintf("  Genes after pre-filtering: %d", nrow(dds)))

  # DESeq2 LRT
  dds <- DESeq(dds, test = "LRT", reduced = reduced_formula,
               fitType = "parametric", quiet = TRUE)

  res    <- results(dds, alpha = padj_thresh)
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    arrange(padj)

  vst_mat <- assay(vst(dds, blind = FALSE))

  list(
    results = res_df,
    dds     = dds,
    vst     = vst_mat,
    degs    = res_df %>% filter(!is.na(padj) & padj < padj_thresh)
  )
}

# ── 3. Run DESeq2 for Cyto ────────────────────────────────────────────────────
message("\n=== Running DESeq2 for Cyto-opto-IRE1α ===")
cyto_out <- run_deseq2_spline(cyto_raw, timepoints, df = 4, padj_thresh = 1e-3)
message(sprintf("  Cyto DEGs (padj < 1e-3): %d", nrow(cyto_out$degs)))

# ── 4. Run DESeq2 for ER ──────────────────────────────────────────────────────
message("\n=== Running DESeq2 for ER-opto-IRE1α ===")
er_out <- run_deseq2_spline(er_raw, timepoints, df = 4, padj_thresh = 1e-3)
message(sprintf("  ER DEGs (padj < 1e-3): %d", nrow(er_out$degs)))

# ── 5. Combined VST matrix (all 48 samples, for UMAP and heatmap) ─────────────
message("\nComputing VST on combined count matrix (48 samples)...")
common_genes    <- intersect(rownames(cyto_raw), rownames(er_raw))
combined_counts <- cbind(cyto_raw[common_genes, ], er_raw[common_genes, ])

combined_coldata <- data.frame(
  sample    = colnames(combined_counts),
  condition = c(rep("Cyto", 24), rep("ER", 24)),
  time      = rep(timepoints, 2),
  row.names = colnames(combined_counts)
)

dds_combined <- DESeqDataSetFromMatrix(
  countData = round(combined_counts),
  colData   = combined_coldata,
  design    = ~ condition
)
keep_combined <- rowSums(counts(dds_combined) >= 10) >= 1
dds_combined  <- dds_combined[keep_combined, ]
dds_combined  <- estimateSizeFactors(dds_combined)
vst_combined  <- assay(vst(dds_combined, blind = TRUE))
message(sprintf("  VST matrix: %d genes × %d samples", nrow(vst_combined), ncol(vst_combined)))

# ── 6. Save outputs ───────────────────────────────────────────────────────────
message("\nSaving results...")
write.csv(cyto_out$results, file.path(data_dir, "DESeq2_Cyto_results.csv"), row.names = FALSE)
write.csv(er_out$results,   file.path(data_dir, "DESeq2_ER_results.csv"),   row.names = FALSE)
write.csv(cyto_out$degs,    file.path(data_dir, "DESeq2_Cyto_DEGs.csv"),    row.names = FALSE)
write.csv(er_out$degs,      file.path(data_dir, "DESeq2_ER_DEGs.csv"),      row.names = FALSE)

vst_df       <- as.data.frame(vst_combined)
vst_df$gene  <- rownames(vst_combined)
vst_df       <- vst_df[, c("gene", setdiff(colnames(vst_df), "gene"))]
write.csv(vst_df, file.path(data_dir, "vst_combined.csv"), row.names = FALSE)

message("\nScript 01 complete.")
message(sprintf("  Cyto DEGs: %d | ER DEGs: %d", nrow(cyto_out$degs), nrow(er_out$degs)))
