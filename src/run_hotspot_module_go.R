#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

modules_file <- args[1]
autocorrelations_file <- args[2]
species <- tolower(args[3])
out_prefix <- args[4]

org_db <- if (species == "mouse") {
  org.Mm.eg.db
} else if (species == "human") {
  org.Hs.eg.db
} else {
  stop("Species must be 'mouse' or 'human'")
}

modules <- read.delim(
  modules_file,
  row.names = 1,
  check.names = FALSE,
  sep = "\t"
)

autocorrelations <- read.delim(
  autocorrelations_file,
  row.names = 1,
  check.names = FALSE,
  sep = "\t"
)

universe_map <- bitr(
  rownames(autocorrelations),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org_db
) %>%
  distinct(SYMBOL, .keep_all = TRUE)

universe <- unique(universe_map$ENTREZID)
module_ids <- sort(unique(modules$Module[modules$Module != -1]))
all_results <- list()

for (module_id in module_ids) {
  module_symbols <- rownames(modules)[modules$Module == module_id]
  module_genes <- unique(
    universe_map$ENTREZID[universe_map$SYMBOL %in% module_symbols]
  )

  if (length(module_genes) < 10) {
    warning("Skipping module ", module_id, ": fewer than 10 mapped genes")
    next
  }

  enrichment <- enrichGO(
    gene = module_genes,
    universe = universe,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    pAdjustMethod = "BH",
    readable = TRUE
  )

  result <- as.data.frame(enrichment)
  if (nrow(result) == 0) {
    next
  }

  result <- result %>% mutate(module = module_id, .before = 1)
  all_results[[as.character(module_id)]] <- result

  plot <- dotplot(
    enrichment,
    showCategory = 10,
    x = "GeneRatio",
    color = "p.adjust",
    size = "Count"
  ) +
    ggtitle(paste("Module", module_id)) +
    theme_bw(base_size = 11)

  ggsave(
    paste0(out_prefix, ".module_", module_id, ".pdf"),
    plot,
    width = 7,
    height = 5
  )
}

write_tsv(bind_rows(all_results), paste0(out_prefix, ".tsv"))
