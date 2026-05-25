# =============================================================================
# Script 02: GSEA Pathway Analysis
# =============================================================================
# Description:
#   Gene Set Enrichment Analysis (GSEA) using Hallmark and KEGG gene sets
#   for all three comparisons: XBP1-s vs WT, XBP1-es vs WT, XBP1-s vs XBP1-es.
#   Genes are ranked by sign(log2FC) * -log10(pvalue).
#
# Input:
#   data/DEG_results/XBP1s_vs_WT_full.csv
#   data/DEG_results/XBP1es_vs_WT_full.csv
#   data/DEG_results/XBP1s_vs_XBP1es_full.csv
#
# Output:
#   data/GSEA_results/XBP1s_vs_WT_gsea_results.csv
#   data/GSEA_results/XBP1es_vs_WT_gsea_results.csv
#   data/GSEA_results/XBP1s_vs_XBP1es_gsea_results.csv
#   figures/fig4a_gsea_XBP1s_vs_WT.pdf
#   figures/fig4b_gsea_XBP1es_vs_WT.pdf
#   figures/fig4c_gsea_XBP1s_vs_XBP1es.pdf
#
# Software: R >= 4.2, clusterProfiler >= 4.6, msigdbr >= 7.5, ggplot2 >= 3.4
# =============================================================================

library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(readr)

set.seed(42)

# ── Gene sets: Hallmark + KEGG (human) ───────────────────────────────────────
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, gene_symbol)
kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>%
  select(gs_name, gene_symbol)
gene_sets <- bind_rows(hallmark, kegg)

# ── Helper: build ranked gene list ───────────────────────────────────────────
make_ranked_list <- function(df) {
  df <- df %>%
    filter(!is.na(log2FoldChange), !is.na(pvalue), !is.na(gene_name)) %>%
    mutate(rank_stat = sign(log2FoldChange) * -log10(pvalue + 1e-300)) %>%
    arrange(desc(rank_stat)) %>%
    distinct(gene_name, .keep_all = TRUE)
  ranked <- setNames(df$rank_stat, df$gene_name)
  ranked
}

# ── Helper: run GSEA and return tidy results ──────────────────────────────────
run_gsea <- function(ranked, gene_sets, label) {
  res <- GSEA(
    geneList     = ranked,
    TERM2GENE    = gene_sets,
    minGSSize    = 15,
    maxGSSize    = 500,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    eps          = 0,
    seed         = TRUE
  )
  as.data.frame(res) %>%
    mutate(comparison = label) %>%
    arrange(NES)
}

# ── Helper: dotplot ───────────────────────────────────────────────────────────
gsea_dotplot <- function(df, title, top_n = 15) {
  # Select top activated and top suppressed
  top_up   <- df %>% filter(NES > 0) %>% slice_max(NES,  n = top_n)
  top_down <- df %>% filter(NES < 0) %>% slice_min(NES,  n = top_n)
  plot_df  <- bind_rows(top_up, top_down) %>%
    mutate(ID = gsub("HALLMARK_|KEGG_", "", ID),
           ID = gsub("_", " ", ID),
           ID = factor(ID, levels = rev(unique(ID))))

  ggplot(plot_df, aes(x = NES, y = ID, size = setSize,
                      color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#E05C2A", high = "#CCCCCC",
                         name = "padj", limits = c(0, 0.05)) +
    scale_size_continuous(name = "Gene set size", range = c(1.5, 5)) +
    geom_vline(xintercept = 0, lty = 2, lwd = 0.4, color = "grey50") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_classic(base_size = 6) +
    theme(
      text             = element_text(family = "sans", size = 6),
      axis.text        = element_text(size = 5.5),
      axis.title       = element_text(size = 6),
      legend.text      = element_text(size = 5),
      legend.title     = element_text(size = 5.5),
      plot.title       = element_text(size = 6.5, face = "bold"),
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.4)
    )
}

# ── Load DEG results ──────────────────────────────────────────────────────────
s_wt  <- read_csv("data/DEG_results/XBP1s_vs_WT_full.csv")
es_wt <- read_csv("data/DEG_results/XBP1es_vs_WT_full.csv")
s_es  <- read_csv("data/DEG_results/XBP1s_vs_XBP1es_full.csv")

# ── Run GSEA ──────────────────────────────────────────────────────────────────
ranked_s_wt  <- make_ranked_list(s_wt)
ranked_es_wt <- make_ranked_list(es_wt)
ranked_s_es  <- make_ranked_list(s_es)

gsea_s_wt  <- run_gsea(ranked_s_wt,  gene_sets, "XBP1s_vs_WT")
gsea_es_wt <- run_gsea(ranked_es_wt, gene_sets, "XBP1es_vs_WT")
gsea_s_es  <- run_gsea(ranked_s_es,  gene_sets, "XBP1s_vs_XBP1es")

# ── Save CSV results ──────────────────────────────────────────────────────────
dir.create("data/GSEA_results", recursive = TRUE, showWarnings = FALSE)
write_csv(gsea_s_wt,  "data/GSEA_results/XBP1s_vs_WT_gsea_results.csv")
write_csv(gsea_es_wt, "data/GSEA_results/XBP1es_vs_WT_gsea_results.csv")
write_csv(gsea_s_es,  "data/GSEA_results/XBP1s_vs_XBP1es_gsea_results.csv")

# ── Save dotplot figures ──────────────────────────────────────────────────────
dir.create("figures", recursive = TRUE, showWarnings = FALSE)

ggsave("figures/fig4a_gsea_XBP1s_vs_WT.pdf",
       gsea_dotplot(gsea_s_wt,  "XBP1-s vs WT"),
       width = 3.5, height = 4, device = cairo_pdf)

ggsave("figures/fig4b_gsea_XBP1es_vs_WT.pdf",
       gsea_dotplot(gsea_es_wt, "XBP1-es vs WT"),
       width = 3.5, height = 4, device = cairo_pdf)

ggsave("figures/fig4c_gsea_XBP1s_vs_XBP1es.pdf",
       gsea_dotplot(gsea_s_es,  "XBP1-s vs XBP1-es"),
       width = 3.5, height = 4, device = cairo_pdf)

cat("GSEA analysis complete.\n")
