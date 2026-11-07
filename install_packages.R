# Run this once from the project folder.
dir.create("library", showWarnings = FALSE)
.libPaths(c(normalizePath("library"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = "library")
}

packages <- c("variancePartition", "data.table", "ggplot2", "cowplot")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  BiocManager::install(missing, lib = "library", update = FALSE, ask = FALSE)
}
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
