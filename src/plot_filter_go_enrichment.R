#!/usr/bin/env Rscript
# Plot previously calculated GO enrichment results.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    "Usage: plot_filter_go_enrichment.R ",
    "<enrichment.tsv> <output.pdf> <passing|not_passing> ",
    "<sample_label> [<sample_label> ...]"
  )
}

enrichment_file <- args[1]
output_file <- args[2]
status_to_plot <- args[3]
sample_labels <- args[-c(1, 2, 3)]

if (!status_to_plot %in% c("passing", "not_passing")) {
  stop("Status must be 'passing' or 'not_passing'")
}

combined <- read_tsv(enrichment_file, show_col_types = FALSE)
required_columns <- c(
  "sample",
  "status",
  "Description",
  "GeneRatio_numeric",
  "Count",
  "pvalue",
  "p.adjust"
)
missing_columns <- setdiff(required_columns, colnames(combined))
if (length(missing_columns) > 0) {
  stop(
    enrichment_file,
    " is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

top_results <- combined %>%
  filter(
    status == status_to_plot,
    !is.na(p.adjust)
  ) %>%
  group_by(sample) %>%
  arrange(p.adjust, desc(GeneRatio_numeric), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  mutate(
    Description = str_wrap(Description, width = 42),
    sample = factor(sample, levels = sample_labels)
  )

plot_title <- "GO Biological Process top 20 terms per sample"

max_count <- if (nrow(top_results) > 0) {
  max(top_results$Count, na.rm = TRUE)
} else {
  1
}
adjusted_p_limits <- if (nrow(top_results) > 0) {
  range(top_results$p.adjust, na.rm = TRUE)
} else {
  c(0, 1)
}

make_sample_panel <- function(sample_label) {
  sample_results <- top_results %>%
    filter(as.character(sample) == sample_label)
  legend_position <- if (sample_label == tail(sample_labels, 1)) {
    "right"
  } else {
    "none"
  }

  if (nrow(sample_results) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "None", size = 5) +
        xlim(-1, 1) +
        ylim(-1, 1) +
        labs(title = sample_label, x = NULL, y = NULL) +
        theme_bw(base_size = 13) +
        theme(
          text = element_text(face = "plain", color = "black"),
          plot.title = element_text(
            size = 13,
            face = "plain",
            color = "black",
            hjust = 0.5
          ),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        )
    )
  }

  ggplot(
    sample_results,
    aes(
      x = GeneRatio_numeric,
      y = reorder(Description, GeneRatio_numeric),
      size = Count,
      fill = p.adjust
    )
  ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.95
    ) +
    scale_fill_gradient(
      low = "#333333",
      high = "#D9D9D9",
      trans = "log10",
      limits = adjusted_p_limits,
      name = "Adj P value"
    ) +
    scale_size_continuous(
      range = c(2, 8),
      limits = c(0, max_count),
      name = "Gene count"
    ) +
    guides(
      fill = guide_colorbar(order = 1),
      size = guide_legend(
        order = 2,
        override.aes = list(
          fill = "white",
          color = "black"
        )
      )
    ) +
    labs(title = sample_label, x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 13) +
    theme(
      text = element_text(face = "plain", color = "black"),
      plot.title = element_text(
        size = 13,
        face = "plain",
        color = "black",
        hjust = 0.5
      ),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(
        size = 11,
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 8.5),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      legend.position = legend_position,
      panel.grid = element_blank()
    )
}

sample_plots <- lapply(sample_labels, make_sample_panel)
plot <- wrap_plots(sample_plots, nrow = 1, guides = "keep") +
  plot_annotation(
    title = plot_title,
    theme = theme(
      plot.title = element_text(
        size = 16,
        face = "plain",
        color = "black",
        hjust = 0.5
      )
    )
  )

ggsave(
  output_file,
  plot,
  width = 4 * length(sample_labels) + 4,
  height = 6
)
