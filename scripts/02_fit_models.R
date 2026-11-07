.libPaths(c("library", .libPaths()))
library(variancePartition)
library(data.table)

data <- readRDS("data/processed/immvar.rds")
expression <- data$expression[, data$metadata$included]
metadata <- droplevels(data$metadata[data$metadata$included, ])
metadata$Age <- as.numeric(scale(metadata$Age))
stopifnot(identical(colnames(expression), rownames(metadata)))
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)

# Fit every retained probe. The same donor can contribute several arrays.
combined_formula <- ~ Age + (1 | Sex) + (1 | CellType) + (1 | Batch) + (1 | Individual)
within_formula <- ~ Age + (1 | Sex) + (1 | Batch) + (1 | Individual)

fit_and_save <- function(expression, metadata, formula, name) {
  stopifnot(identical(colnames(expression), rownames(metadata)))
  message("Fitting ", name, ": ", nrow(expression), " probes, ", ncol(expression), " arrays")
  # A fresh fit avoids boundary solutions inherited from the package's initial fit.
  result <- fitExtractVarPartModel(expression, formula, metadata,
    REML = FALSE, useInitialFit = FALSE, BPPARAM = SerialParam())
  result <- as.data.frame(result)
  stopifnot(identical(rownames(result), rownames(expression)),
            all(is.finite(as.matrix(result))), all(as.matrix(result) >= 0),
            max(abs(rowSums(result) - 1)) < 1e-6)
  saveRDS(result, paste0("results/models/", name, ".rds"))
  fwrite(data.frame(probe_id = rownames(result), result), paste0("results/tables/", name, ".csv"))
  message("Finished ", name)
}

# With no arguments, run all three models. A name can be supplied to run just one.
models <- commandArgs(trailingOnly = TRUE)
if (!length(models)) models <- c("combined", "CD4_T", "Monocyte")
stopifnot(all(models %in% c("combined", "CD4_T", "Monocyte")))
if ("combined" %in% models) {
  fit_and_save(expression, metadata, combined_formula, "combined")
}
for (cell in intersect(c("CD4_T", "Monocyte"), models)) {
  samples <- metadata$CellType == cell
  fit_and_save(expression[, samples], droplevels(metadata[samples, ]), within_formula, cell)
}
capture.output(sessionInfo(), file = "results/sessionInfo.txt")
