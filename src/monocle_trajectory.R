#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(ggplot2)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: run_monocle3.R input_counts.tsv out_prefix [num_dim=50]")
}

input_tsv <- args[1]
out_prefix <- args[2]
num_dim <- ifelse(length(args) >= 3, as.integer(args[3]), 50)

counts_df <- read_tsv(input_tsv, show_col_types = FALSE)

cell_names <- counts_df[[1]]

counts_mat <- Matrix(data.matrix(counts_df[, -1]), sparse = TRUE)
rownames(counts_mat) <- cell_names

# Monocle expects genes x cells
counts_mat <- t(counts_mat)

gene_metadata <- data.frame(
  gene_short_name = rownames(counts_mat),
  row.names = rownames(counts_mat)
)

cell_metadata <- data.frame(
  cell = colnames(counts_mat),
  row.names = colnames(counts_mat)
)

cds <- new_cell_data_set(
  expression_data = counts_mat,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- preprocess_cds(
  cds,
  num_dim = num_dim,
  method = "PCA",
  norm_method = "log"
)

cds <- reduce_dimension(
  cds,
  reduction_method = "UMAP",
  preprocess_method = "PCA"
)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

# Placeholder: use first cell in matrix as root
# TODO: replace with biologically meaningful root selection
root_cell <- colnames(cds)[1]

cds <- order_cells(
  cds,
  root_cells = root_cell
)

saveRDS(
  cds,
  paste0(out_prefix, ".monocle.cds.rds")
)

cell_out <- as.data.frame(colData(cds))
cell_out$cell <- rownames(cell_out)
cell_out$pseudotime <- pseudotime(cds)

write_tsv(
  cell_out,
  paste0(out_prefix, ".monocle.cell_metadata.tsv")
)

theme_publication <- function() {
  theme_classic(base_size = 14) +
    theme(
      legend.title     = element_text(size = 12, face = "bold"),
      legend.text      = element_text(size = 11),
      axis.title       = element_text(size = 13),
      axis.text        = element_text(size = 11),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# Cluster plot
p <- plot_cells(
  cds,
  color_cells_by                = "cluster",
  label_groups_by_cluster       = TRUE,
  label_leaves                  = FALSE,
  label_branch_points           = FALSE,
  cell_size                     = 1.2,
  trajectory_graph_segment_size = 1.2
) + theme_publication()

ggsave(
  paste0(out_prefix, ".monocle.clusters.pdf"),
  p,
  width = 7,
  height = 6
)

# Partition plot
p <- plot_cells(
  cds,
  color_cells_by                = "partition",
  label_groups_by_cluster       = FALSE,
  label_leaves                  = FALSE,
  label_branch_points           = FALSE,
  cell_size                     = 1.2,
  trajectory_graph_segment_size = 1.2
) + theme_publication()

ggsave(
  paste0(out_prefix, ".monocle.partitions.pdf"),
  p,
  width = 7,
  height = 6
)

# Pseudotime plot
p <- plot_cells(
  cds,
  color_cells_by                = "pseudotime",
  label_cell_groups             = FALSE,
  label_leaves                  = TRUE,
  label_branch_points           = TRUE,
  graph_label_size              = 4,
  cell_size                     = 1.2,
  trajectory_graph_segment_size = 1.2
) + theme_publication()

ggsave(
  paste0(out_prefix, ".monocle.pseudotime.pdf"),
  p,
  width = 7,
  height = 6
)

# Trajectory graph plot
p <- plot_cells(
  cds,
  color_cells_by                = "cluster",
  label_cell_groups             = FALSE,
  label_leaves                  = TRUE,
  label_branch_points           = TRUE,
  graph_label_size              = 4,
  cell_size                     = 1.2,
  trajectory_graph_segment_size = 1.2
) + theme_publication()

ggsave(
  paste0(out_prefix, ".monocle.trajectory_graph.pdf"),
  p,
  width = 7,
  height = 6
)

message("Root cell used: ", root_cell)
