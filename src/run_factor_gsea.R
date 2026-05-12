#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)

L_file <- args[1]
species <- tolower(args[2])
out_prefix <- args[3]

# Read L matrix
L <- read.delim(
  L_file,
  row.names = 1,
  check.names = FALSE,
  sep = "\t"
)

all_genes <- colnames(L)

# Convert SYMBOL to ENTREZID
gene_map <- bitr(
  all_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = if (species == "mouse") {
    org.Mm.eg.db
  } else {
    org.Hs.eg.db
  }
)

gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]

message("Mapped ", nrow(gene_map), " / ", length(all_genes), " gene symbols to Entrez IDs")

mapped_symbols <- gene_map$SYMBOL
L_mapped <- L[, mapped_symbols, drop = FALSE]

colnames(L_mapped) <- gene_map$ENTREZID[
  match(colnames(L_mapped), gene_map$SYMBOL)
]

# Run GSEA per factor
all_results <- list()

for (factor_name in rownames(L_mapped)) {
  scores <- as.numeric(L_mapped[factor_name, ])
  names(scores) <- colnames(L_mapped)

  keep <- is.finite(scores) & !is.na(scores) & scores != 0
  scores <- scores[keep]

  scores_df <- tibble(
    ENTREZID = names(scores),
    score = scores
  ) %>%
    group_by(ENTREZID) %>%
    slice_max(order_by = abs(score), n = 1, with_ties = FALSE) %>%
    ungroup()

  scores <- scores_df$score
  names(scores) <- scores_df$ENTREZID
  scores <- sort(scores, decreasing = TRUE)

  if (length(scores) < 20) {
    warning("Skipping ", factor_name, ": too few mapped genes")
    next
  }

  gsea <- gseGO(
    geneList = scores,
    OrgDb = if (species == "mouse") {
      org.Mm.eg.db
    } else {
      org.Hs.eg.db
    },
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    verbose = FALSE,
    eps = 0
  )

  res <- as.data.frame(gsea)

  if (nrow(res) == 0) {
    next
  }

  res <- res %>%
    mutate(
      factor = factor_name,
      direction = ifelse(NES > 0, "positive_loading", "negative_loading")
    )

  all_results[[factor_name]] <- res

  # Optionally plot factor dotplots independently
  p <- dotplot(
    gsea,
    showCategory = 10,
    x = "GeneRatio",
    color = "p.adjust",
    size = "Count"
  ) +
    ggtitle(factor_name) +
    theme_bw(base_size = 11)

  ggsave(
    paste0(out_prefix, ".", factor_name, ".pdf"),
    p,
    width = 6,
    height = 5
  )
}

combined <- bind_rows(all_results)

write_tsv(combined, paste0(out_prefix, ".tsv"))

# # Select top terms per factor
# top_terms <- combined %>%
#   filter(!is.na(p.adjust)) %>%
#   mutate(
#     GeneRatio_num = abs(NES),
#     Count = setSize
#   ) %>%
#   group_by(factor) %>%
#   arrange(p.adjust, desc(abs(NES)), .by_group = TRUE) %>%
#   slice_head(n = 10) %>%
#   ungroup() %>%
#   mutate(
#     Description = str_wrap(Description, width = 38),
#     factor = factor(factor, levels = rownames(L))
#   )

# # Classic dot plot
# p <- ggplot(
#   top_terms,
#   aes(
#     x = GeneRatio_num,
#     y = fct_reorder(Description, GeneRatio_num),
#     size = Count,
#     color = p.adjust
#   )
# ) +
#   geom_point(alpha = 0.95) +
#   facet_wrap(~ factor, scales = "free_y", ncol = 2) +
#   scale_color_gradient(
#     low = "red",
#     high = "blue",
#     name = "p.adjust"
#   ) +
#   scale_size_continuous(
#     range = c(2, 8),
#     name = "Gene set size"
#   ) +
#   labs(
#     x = "|NES|",
#     y = NULL
#   ) +
#   theme_bw(base_size = 11) +
#   theme(
#     strip.text = element_text(face = "bold", size = 10),
#     axis.text.y = element_text(size = 8),
#     axis.text.x = element_text(size = 9),
#     axis.title.x = element_text(size = 11),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 9),
#     panel.grid.major = element_line(linewidth = 0.3),
#     panel.grid.minor = element_line(linewidth = 0.2),
#     legend.position = "right"
#   )

# ggsave(
#   paste0(out_prefix, ".pdf"),
#   p,
#   width = 11,
#   height = 8.5
# )
