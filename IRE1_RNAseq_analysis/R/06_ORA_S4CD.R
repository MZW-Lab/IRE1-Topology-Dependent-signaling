# =============================================================================
# Script 06: GO Over-Representation Analysis on DESeq2 DEGs (Fig. S4C–D)
# =============================================================================
# Manuscript: "A topology-dependent IRE1α signaling axis generates a cryptic
#              XBP1 isoform present in human cancers"
# Journal: Nature Chemical Biology
#
# Purpose:
#   Perform GO Biological Process over-representation analysis (ORA) on the
#   top differentially expressed genes from DESeq2 spline analysis for each
#   condition. Both conditions show enrichment for canonical ER stress response
#   terms, validating that each optogenetic construct engages the expected
#   IRE1α signaling program.
#     Fig. S4C — ORA for ER-opto-IRE1α DEGs
#     Fig. S4D — ORA for Cyto-opto-IRE1α DEGs
#
# Input:
#   data/DESeq2_Cyto_DEGs.csv  — Significant DEGs for Cyto (Script 01 output)
#   data/DESeq2_ER_DEGs.csv    — Significant DEGs for ER (Script 01 output)
#
# Output:
#   figures/FigS4CD_ORA.pdf
#
# Method:
#   - clusterProfiler enrichGO(), ontology = "BP"
#   - BH-FDR correction, padj < 0.05
#   - Top 15 terms shown per condition (ranked by gene ratio)
#   - Universe = all expressed genes tested in DESeq2
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ── 0. Paths ──────────────────────────────────────────────────────────────────
data_dir <- "data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

# ── 1. Load DEG lists ─────────────────────────────────────────────────────────
message("Loading DEG lists...")
cyto_degs <- read.csv(file.path(data_dir, "DESeq2_Cyto_DEGs.csv"),
                      stringsAsFactors = FALSE)
er_degs   <- read.csv(file.path(data_dir, "DESeq2_ER_DEGs.csv"),
                      stringsAsFactors = FALSE)

# Load full results for universe definition
cyto_all <- read.csv(file.path(data_dir, "DESeq2_Cyto_results.csv"),
                     stringsAsFactors = FALSE)
er_all   <- read.csv(file.path(data_dir, "DESeq2_ER_results.csv"),
                     stringsAsFactors = FALSE)

message(sprintf("Cyto DEGs: %d | ER DEGs: %d", nrow(cyto_degs), nrow(er_degs)))

# ── 2. Convert gene symbols to Entrez IDs ────────────────────────────────────
message("Converting gene symbols to Entrez IDs...")

convert_to_entrez <- function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
}

cyto_entrez     <- convert_to_entrez(cyto_degs$gene)
er_entrez       <- convert_to_entrez(er_degs$gene)
universe_entrez <- convert_to_entrez(cyto_all$gene)   # all tested genes

message(sprintf("Cyto: %d/%d genes mapped | ER: %d/%d genes mapped",
                length(cyto_entrez), nrow(cyto_degs),
                length(er_entrez),   nrow(er_degs)))

# ── 3. Run GO ORA ─────────────────────────────────────────────────────────────
run_ora <- function(entrez_ids, universe, label) {
  message(sprintf("Running GO ORA for %s...", label))
  go <- enrichGO(
    gene          = entrez_ids,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE,
    minGSSize     = 10,
    maxGSSize     = 500
  )
  # Simplify to reduce redundancy
  go_simplified <- simplify(go, cutoff = 0.7, by = "p.adjust", select_fun = min)
  message(sprintf("  %s: %d significant GO terms (after simplify)", label,
                  nrow(go_simplified@result %>% filter(p.adjust < 0.05))))
  go_simplified
}

cyto_go <- run_ora(cyto_entrez, universe_entrez, "Cyto-opto-IRE1a")
er_go   <- run_ora(er_entrez,   universe_entrez, "ER-opto-IRE1a")

# ── 4. Plot helper ────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = 7) +
  theme(
    text             = element_text(family = "sans", size = 7),
    axis.title       = element_text(size = 7),
    axis.text.x      = element_text(size = 6),
    axis.text.y      = element_text(size = 6),
    legend.title     = element_text(size = 6),
    legend.text      = element_text(size = 5),
    legend.key.size  = unit(0.25, "cm"),
    plot.title       = element_text(size = 7, face = "bold"),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

make_ora_plot <- function(go_obj, title, bar_color, n_terms = 15) {
  df <- go_obj@result %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(n_terms) %>%
    mutate(
      GeneRatio_num = sapply(GeneRatio, function(x) {
        parts <- strsplit(x, "/")[[1]]
        as.numeric(parts[1]) / as.numeric(parts[2])
      }),
      Description = factor(Description, levels = rev(Description))
    )

  ggplot(df, aes(x = GeneRatio_num, y = Description, fill = p.adjust)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = bar_color, high = "grey85",
                        name = "padj", limits = c(0, 0.05)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = title, x = "Gene Ratio", y = NULL) +
    base_theme +
    theme(axis.text.y = element_text(size = 5.5))
}

pC <- make_ora_plot(er_go,   "Fig. S4C: ER-opto-IRE1a DEGs",   "#f6992d")
pD <- make_ora_plot(cyto_go, "Fig. S4D: Cyto-opto-IRE1a DEGs", "#85ab8e")

# ── 5. Save ───────────────────────────────────────────────────────────────────
tmp_path <- file.path(tempdir(), "FigS4CD_ORA.pdf")
pdf(tmp_path, width = 8, height = 5, useDingbats = FALSE)
print(pC + pD + plot_layout(ncol = 2))
dev.off()

out_path <- file.path(fig_dir, "FigS4CD_ORA.pdf")
file.copy(tmp_path, out_path, overwrite = TRUE)
message(sprintf("Saved: %s", out_path))
message("Script 06 complete.")
