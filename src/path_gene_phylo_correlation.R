#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(PATH)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
expr_file <- args[2]
out_file  <- args[3]

# Read tree; supports Newick or Nexus by extension
tree <- if (grepl("\\.(nex|nexus)$", tree_file, ignore.case = TRUE)) {
  ape::read.nexus(tree_file)
} else {
  ape::read.tree(tree_file)
}

# Read expression matrix
# Assumes:
#   - first column is cell ID
#   - remaining columns are genes
#   - rows are cells
expr_df <- read.table(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
)

# Convert to matrix with rownames as cell IDs and columns as genes
cell_ids <- expr_df[[1]]
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
rownames(expr_mat) <- cell_ids
mode(expr_mat) <- "numeric"

# Check that tree tips are present in the expression matrix
missing_in_expr <- setdiff(tree$tip.label, rownames(expr_mat))
if (length(missing_in_expr) > 0) {
  stop(
    sprintf(
      "These tree tips are missing from the expression matrix: %s",
      paste(missing_in_expr, collapse = ", ")
    )
  )
}

# Reorder matrix rows to match tree tip order exactly
expr_mat <- expr_mat[tree$tip.label, , drop = FALSE]

# Build phylogenetic weight matrix.
Winv <- PATH::inv_tree_dist(tree, node = FALSE, norm = FALSE)

# Compute phylogenetic correlations
xc <- PATH::xcor(expr_mat, Winv)

# Extract gene-level auto-correlation statistics from the diagonal
phy_cor <- diag(xc$phy_cor)
z_score <- diag(xc$Z.score)
p_value <- diag(xc$one.sided.pvalue)

# Multiple-testing correction
padj <- p.adjust(p_value, method = "BH")

results <- data.frame(
  gene = colnames(expr_mat),
  phy_cor = phy_cor,
  z_score = z_score,
  p_value = p_value,
  padj = padj,
  pass = (padj < 0.05),
  stringsAsFactors = FALSE
)

write.table(
  results,
  file = out_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
