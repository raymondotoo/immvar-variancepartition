# Install the packages this project needs. Run it once, from this folder,
# before knitting any notebook:
#
#   Rscript install_packages.R        (terminal)
#   source("install_packages.R")      (RStudio console)
#
# Packages go into ./library, so your system-wide R library is not changed.
# This stays a plain R script because knitting needs rmarkdown to be installed first.

dir.create("library", showWarnings = FALSE)
.libPaths(c(normalizePath("library"), .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = "library")
}

packages <- c("variancePartition", "data.table", "ggplot2", "cowplot",
              "rmarkdown", "knitr")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  BiocManager::install(missing, lib = "library", update = FALSE, ask = FALSE)
}

stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
message("All packages are installed.")
