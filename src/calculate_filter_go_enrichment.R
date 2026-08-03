#!/usr/bin/env Rscript
# Calculate GO enrichment results for filtered gene sets.

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  stop(
    "Usage: calculate_filter_go_enrichment.R ",
    "<passing|not_passing|both> <out_prefix> <mouse|human> ",
    "<sample_label> <moran.tsv>"
  )
}

status <- args[1]
out_prefix <- args[2]
species <- tolower(args[3])
sample_args <- args[-c(1, 2, 3)]

if (!status %in% c("passing", "not_passing", "both")) {
  stop("Status must be 'passing', 'not_passing', or 'both'")
}
org_db <- if (species == "mouse") {
  org.Mm.eg.db
} else if (species == "human") {
  org.Hs.eg.db
} else {
  stop("Species must be 'mouse' or 'human'")
}

parse_gene_ratio <- function(values) {
  vapply(
    strsplit(as.character(values), "/", fixed = TRUE),
    function(parts) as.numeric(parts[1]) / as.numeric(parts[2]),
    numeric(1)
  )
}

all_results <- list()
sample_labels <- character()
status_values <- if (status == "both") {
  c("passing", "not_passing")
} else {
  status
}

for (i in seq(1, length(sample_args), by = 2)) {
  sample_label <- sample_args[i]
  moran_file <- sample_args[i + 1]
  sample_labels <- c(sample_labels, sample_label)

  moran <- read_tsv(moran_file, show_col_types = FALSE)
  required_columns <- c("gene", "keep")
  missing_columns <- setdiff(required_columns, colnames(moran))
  if (length(missing_columns) > 0) {
    stop(
      moran_file,
      " is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  moran <- moran %>%
    mutate(
      gene = as.character(gene),
      keep = as.logical(keep)
    ) %>%
    filter(!is.na(gene), !is.na(keep))

  duplicate_rows <- sum(duplicated(moran$gene))
  conflicting_symbols <- moran %>%
    group_by(gene) %>%
    summarize(n_keep_values = n_distinct(keep), .groups = "drop") %>%
    filter(n_keep_values > 1) %>%
    nrow()
  if (duplicate_rows > 0) {
    message(
      sample_label,
      ": collapsed ",
      duplicate_rows,
      " duplicate row(s) by gene symbol using keep=any; ",
      conflicting_symbols,
      " symbol(s) had conflicting keep values"
    )
  }
  moran <- moran %>%
    group_by(gene) %>%
    summarize(keep = any(keep), .groups = "drop")

  gene_map <- bitr(
    unique(moran$gene),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org_db
  ) %>%
    distinct(SYMBOL, .keep_all = TRUE)

  universe <- unique(gene_map$ENTREZID)
  for (current_status in status_values) {
    selected_symbols <- if (current_status == "passing") {
      moran$gene[moran$keep]
    } else {
      moran$gene[!moran$keep]
    }
    selected_genes <- unique(
      gene_map$ENTREZID[gene_map$SYMBOL %in% selected_symbols]
    )

    message(
      sample_label,
      " (",
      current_status,
      "): mapped ",
      length(universe),
      " tested genes and ",
      length(selected_genes),
      " selected genes"
    )

    if (length(selected_genes) < 10) {
      warning(
        sample_label,
        " (",
        current_status,
        "): fewer than 10 selected genes mapped to Entrez IDs"
      )
      next
    }

    enrichment <- enrichGO(
      gene = selected_genes,
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
      warning(
        sample_label,
        " (",
        current_status,
        "): GO enrichment returned no terms"
      )
      next
    }

    result_key <- paste(sample_label, current_status, sep = "_")
    all_results[[result_key]] <- result %>%
      mutate(
        sample = sample_label,
        status = current_status,
        GeneRatio_numeric = parse_gene_ratio(GeneRatio),
        .before = 1
      )
  }
}

combined <- bind_rows(all_results)
write_tsv(combined, paste0(out_prefix, ".tsv"))
