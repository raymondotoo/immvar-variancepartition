.libPaths(c("library", .libPaths()))
library(variancePartition)
library(data.table)
library(ggplot2)

data <- readRDS("data/processed/immvar.rds")
expression <- data$expression[, data$metadata$included]
metadata <- droplevels(data$metadata[data$metadata$included, ])
metadata$Age <- as.numeric(scale(metadata$Age))
models <- list(Combined = readRDS("results/models/combined.rds"),
               CD4_T = readRDS("results/models/CD4_T.rds"),
               Monocyte = readRDS("results/models/Monocyte.rds"))
combined_formula <- ~ Age + (1 | Sex) + (1 | CellType) + (1 | Batch) + (1 | Individual)
within_formula <- ~ Age + (1 | Sex) + (1 | Batch) + (1 | Individual)

# Check that no probes were lost and every decomposition sums to one.
fit_counts <- data.frame()
for (name in names(models)) {
  result <- models[[name]]
  stopifnot(identical(rownames(result), rownames(expression)),
            all(is.finite(as.matrix(result))), all(as.matrix(result) >= 0),
            max(abs(rowSums(result) - 1)) < 1e-6)
  fit_counts <- rbind(fit_counts, data.frame(model = name, attempted = nrow(expression),
                                           fitted = nrow(result), failed = nrow(expression) - nrow(result)))
}
fwrite(fit_counts, "results/tables/model_fit_counts.csv")

# Compare selected probes with independent lmer fits using a different optimizer.
set.seed(20261104)
audit_ids <- unique(c("8177232", "8157524", "7903765", sample(rownames(expression), 50)))
control <- lme4::lmerControl(optimizer = "bobyqa", calc.derivs = FALSE,
                            optCtrl = list(maxfun = 100000))
audit <- list()
residuals <- list()
for (name in names(models)) {
  samples <- if (name == "Combined") rep(TRUE, nrow(metadata)) else metadata$CellType == name
  sample_data <- droplevels(metadata[samples, ])
  formula <- if (name == "Combined") combined_formula else within_formula
  for (id in audit_ids) {
    sample_data$expression <- as.numeric(expression[id, samples])
    # Zero variance components produce expected singular-fit messages.
    fit <- suppressMessages(lme4::lmer(update(formula, expression ~ .), sample_data,
                                        REML = FALSE, control = control))
    independent <- calcVarPart(fit)
    original <- as.numeric(models[[name]][id, names(independent)])
    audit[[paste(name, id)]] <- data.frame(model = name, probe_id = id,
      component = names(independent), pipeline_fraction = original,
      independent_fraction = unname(independent),
      absolute_difference = abs(original - unname(independent)))
    if (name == "Combined" && id %in% c("8177232", "8157524", "7903765")) {
      residuals[[id]] <- data.frame(probe_id = id, fitted = fitted(fit),
                                    standardized_residual = residuals(fit) / sigma(fit))
    }
  }
}
audit <- rbindlist(audit)
fwrite(audit, "results/tables/independent_fit_check.csv")
stopifnot(max(audit$absolute_difference) < .01)

residuals <- rbindlist(residuals)
residuals$gene <- c("8177232" = "KDM5D", "8157524" = "TLR4", "7903765" = "GSTM1")[residuals$probe_id]
p1 <- ggplot(residuals, aes(fitted, standardized_residual)) +
  geom_point(alpha = .3, size = .7) + geom_hline(yintercept = 0) +
  facet_wrap(~ gene, scales = "free_x") + theme_minimal() +
  labs(x = "Fitted expression", y = "Standardized residual")
p2 <- ggplot(residuals, aes(sample = standardized_residual)) +
  stat_qq(alpha = .4, size = .7) + stat_qq_line() + facet_wrap(~ gene) + theme_minimal() +
  labs(x = "Theoretical normal quantile", y = "Observed residual quantile")
plot <- cowplot::plot_grid(p1, p2, ncol = 1)
ggsave("results/figures/model_diagnostics.png", plot, width = 11, height = 7, dpi = 200, bg = "white")
ggsave("results/figures/model_diagnostics.pdf", plot, width = 11, height = 7, bg = "white")

# Use a random subset, rather than selecting only highly variable probes.
set.seed(20261104)
subset_ids <- sample(rownames(expression), 2000)
fwrite(data.frame(probe_id = subset_ids), "results/tables/sensitivity_probes.csv")
fit_sensitivity <- function(expression, metadata, formula) {
  result <- as.data.frame(fitExtractVarPartModel(expression, formula, metadata,
    REML = FALSE, useInitialFit = FALSE, BPPARAM = SerialParam()))
  stopifnot(identical(rownames(result), rownames(expression)),
            all(is.finite(as.matrix(result))), max(abs(rowSums(result) - 1)) < 1e-6)
  result
}

# Include all 388 paired donors, without using the inconsistent age records.
paired <- droplevels(data$metadata[data$metadata$paired, ])
no_age_formula <- ~ (1 | Sex) + (1 | CellType) + (1 | Batch) + (1 | Individual)
all_paired <- fit_sensitivity(data$expression[subset_ids, paired$sample_id], paired, no_age_formula)

# Compare the two-level sex/cell variables as numeric fixed effects instead.
fixed_metadata <- metadata
fixed_metadata$Sex <- as.numeric(fixed_metadata$Sex == "Male")
fixed_metadata$CellType <- as.numeric(fixed_metadata$CellType == "Monocyte")
fixed_formula <- ~ Age + Sex + CellType + (1 | Batch) + (1 | Individual)
binary_fixed <- fit_sensitivity(expression[subset_ids, ], fixed_metadata, fixed_formula)
sensitivities <- list(Primary_subset = models$Combined[subset_ids, ],
                      All_paired_no_age = all_paired, Binary_fixed = binary_fixed)
summaries <- list()
correlations <- list()
for (name in names(sensitivities)) {
  result <- sensitivities[[name]]
  fwrite(data.frame(probe_id = rownames(result), result), paste0("results/tables/", name, ".csv"))
  summaries[[name]] <- data.frame(model = name, component = names(result),
                                  median_percent = 100 * vapply(result, median, numeric(1)), row.names = NULL)
  if (name == "Primary_subset") next
  components <- intersect(names(result), names(sensitivities$Primary_subset))
  correlations[[name]] <- data.frame(model = name, component = components,
    spearman = vapply(components, function(component) {
      cor(result[[component]], sensitivities$Primary_subset[[component]], method = "spearman")
    }, numeric(1)))
}
summaries <- rbindlist(summaries)
fwrite(summaries, "results/tables/sensitivity_summary.csv")
fwrite(rbindlist(correlations), "results/tables/sensitivity_correlations.csv")
p <- ggplot(summaries, aes(component, median_percent, fill = model)) +
  geom_col(position = "dodge") + theme_minimal() + theme(legend.position = "bottom") +
  labs(x = NULL, y = "Median variance explained (%)", fill = NULL)
ggsave("results/figures/sensitivity_checks.png", p, width = 10, height = 5, dpi = 200, bg = "white")
ggsave("results/figures/sensitivity_checks.pdf", p, width = 10, height = 5, bg = "white")

packages <- c("variancePartition", "lme4", "limma", "BiocParallel", "data.table", "ggplot2", "cowplot")
fwrite(data.frame(package = packages,
  version = vapply(packages, function(p) as.character(packageVersion(p)), character(1))),
  "results/package_versions.csv")
capture.output(sessionInfo(), file = "results/sessionInfo.txt")
writeLines(c(
  "All three main models retained every probe; fractions are finite, nonnegative and sum to one.",
  sprintf("Independent optimizer comparison: %d probes per model; maximum difference %.6f percentage points.",
    length(audit_ids), 100 * max(audit$absolute_difference)),
  "Sensitivity analyses: 2,000 fixed random probes; expanded paired cohort without age; numeric fixed sex/cell effects.",
  "Residual diagnostics are saved separately; these variance fractions are descriptive point estimates."
), "results/checks.txt")
message("Model checks and sensitivity analyses finished.")
