.libPaths(c("library", .libPaths()))
library(data.table)
library(ggplot2)
library(cowplot)

data <- readRDS("data/processed/immvar.rds")
metadata <- droplevels(data$metadata[data$metadata$included, ])

combined <- readRDS("results/models/combined.rds")
t_cells <- readRDS("results/models/CD4_T.rds")
monocytes <- readRDS("results/models/Monocyte.rds")

# Summarize the median fraction for each variance component.
summarize_model <- function(model, name) {
  data.frame(
    model = name,
    component = names(model),
    n_probes = nrow(model),
    median_percent = 100 * vapply(model, median, numeric(1))
  )
}

all_summary <- rbind(
  summarize_model(combined, "Combined"),
  summarize_model(t_cells, "CD4_T"),
  summarize_model(monocytes, "Monocyte")
)
write.csv(all_summary, "results/tables/variance_summary.csv", row.names = FALSE)

# Keep a table of the three examples discussed in the report.
marker_ids <- c(TLR4 = "8157524", KDM5D = "8177232", GSTM1 = "7903765")
markers <- data.frame(
  symbol = names(marker_ids),
  probe_id = unname(marker_ids),
  combined[unname(marker_ids), , drop = FALSE],
  check.names = FALSE
)
write.csv(markers, "results/tables/marker_variance.csv", row.names = FALSE)

colors <- c(
  CellType = "#1F78B4", Individual = "#33A02C", Age = "#FF9900",
  Sex = "#E31A1C", Batch = "#6A3D9A", Residuals = "#B5B5B5"
)
component_names <- c(
  CellType = "Cell type", Individual = "Donor", Age = "Age",
  Sex = "Sex", Batch = "Batch", Residuals = "Residual"
)

# Panel A: median fractions across all probes.
combined_summary <- subset(all_summary, model == "Combined")
combined_summary$component <- factor(
  combined_summary$component,
  levels = c("Batch", "Sex", "Age", "Individual", "Residuals", "CellType")
)
p1 <- ggplot(combined_summary,
             aes(median_percent, component, fill = component)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", median_percent)),
            hjust = -0.15, size = 3) +
  scale_fill_manual(values = colors, guide = "none") +
  scale_y_discrete(labels = component_names) +
  scale_x_continuous(limits = c(0, 62), breaks = c(0, 20, 40, 60)) +
  labs(title = "A  Median across 21,462 probes",
       x = "Variance explained (%)", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

# Panel B: examples where different components stand out.
marker_long <- melt(
  as.data.table(markers),
  id.vars = c("symbol", "probe_id"),
  variable.name = "component",
  value.name = "fraction"
)
marker_long$component <- factor(marker_long$component, levels = names(colors))
marker_long$symbol <- factor(marker_long$symbol, levels = c("GSTM1", "KDM5D", "TLR4"))
p2 <- ggplot(marker_long, aes(100 * fraction, symbol, fill = component)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = colors, labels = component_names) +
  scale_x_continuous(breaks = c(0, 50, 100), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 100)) +
  labs(title = "B  Examples by probe", x = "Variance explained (%)",
       y = NULL, fill = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 1))

legend <- get_legend(p2)
figure <- plot_grid(
  plot_grid(p1, p2 + theme(legend.position = "none"),
            nrow = 1, rel_widths = c(1.05, 1)),
  legend, ncol = 1, rel_heights = c(1, 0.16)
)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
ggsave("results/figures/submission_summary.png", figure,
       width = 7, height = 2.2, dpi = 300, bg = "white")
ggsave("results/figures/submission_summary.pdf", figure,
       width = 7, height = 2.2, bg = "white")

cell_median <- median(combined$CellType) * 100
donor_median <- median(combined$Individual) * 100
writeLines(c(
  "# ImmVar results",
  "",
  sprintf("The combined model was fitted to %s probes from 817 arrays and 379 donors.",
          format(nrow(combined), big.mark = ",")),
  sprintf("Cell type explained a median %.1f%% of expression variation; donor explained %.1f%%.",
          cell_median, donor_median),
  "",
  "![Variance fractions](../results/figures/submission_summary.png)",
  "",
  "TLR4, KDM5D and GSTM1 illustrate cell-type, sex-associated and donor-associated patterns.",
  "Donor-associated variation is not a heritability estimate. The supplied GEO matrix has very little recorded batch-associated variation, and its complete preprocessing history is unknown."
), "report/summary.md")
