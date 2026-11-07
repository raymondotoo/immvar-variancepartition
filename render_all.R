# Render all four notebooks in order, from this folder:
#
#   Rscript render_all.R
#
# This does the same as opening each notebook in RStudio and clicking Knit.
# Each notebook runs in a fresh environment and writes an .html file next to it.
# Outside RStudio, rmarkdown needs pandoc (https://pandoc.org/installing.html).

.libPaths(c("library", .libPaths()))

# CPU cores for model fitting in steps 2 and 4. More cores = faster, same results.
workers <- 1

rmarkdown::render("01_prepare_data.Rmd",  envir = new.env())
rmarkdown::render("02_fit_models.Rmd",    envir = new.env(), params = list(workers = workers))
rmarkdown::render("03_make_figures.Rmd",  envir = new.env())
rmarkdown::render("04_check_results.Rmd", envir = new.env(), params = list(workers = workers))
