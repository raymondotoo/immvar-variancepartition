# ImmVar variancePartition (R Markdown edition)

This project uses [variancePartition](https://doi.org/10.1186/s12859-016-1323-z) to ask
how much of the variation in gene expression is associated with **cell type, donor,
age, sex and batch**. The data are naïve CD4 T cells and monocytes from the ImmVar
study (GEO [GSE56035](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56035)).

The analysis is written as four R Markdown notebooks. Each one explains what it
does and why, shows its results, and can be knitted on its own.

## Quick start

1. Open `immvar.Rproj` in RStudio.
2. Install the packages once. They go into a local `library/` folder:
   ```r
   source("install_packages.R")
   ```
3. Open the notebooks **in order** and click **Knit** on each:
   `01` → `02` → `03` → `04`.

To run everything from a terminal instead:

```sh
Rscript install_packages.R
Rscript render_all.R
```

## The notebooks

| Notebook | What it does | Main outputs |
|---|---|---|
| `01_prepare_data.Rmd` | Downloads the GEO files and verifies their checksums. Builds the sample table, chooses donors, checks for batch effects, converts to log2 and adds gene annotation. | `data/processed/immvar.rds`, QC tables |
| `02_fit_models.Rmd` | Explains the mixed model and **fits one probe by hand** to show how a variance fraction is calculated. Then fits every probe with three models. | `results/models/*.rds`, `results/tables/{combined,CD4_T,Monocyte}.csv` |
| `03_make_figures.Rmd` | Summarises the fractions, shows the full distributions and top probes per component, and makes the summary figure. | `variance_summary.csv`, `marker_variance.csv`, `submission_summary.png/pdf`, `report/summary.md` |
| `04_check_results.Rmd` | Checks that every probe fitted, refits with an independent optimiser, plots residual diagnostics and runs two sensitivity analyses. | Check tables, diagnostic figures, `results/checks.txt` |

Each notebook reads only files written by earlier notebooks, so you can rerun
any single step without rerunning the others.

## Saving time: `refit` and `workers`

Fitting every probe (Step 2) is by far the slowest part. On an Apple-silicon
Mac it took about 30 minutes on one CPU core and about 5 minutes on eight. To
avoid repeating it:

- Steps 2 and 4 **reuse saved fits** in `results/models/` if they exist. To force
  a fresh fit, knit with `refit: true`.
- To use several CPU cores, set `workers`, for example to `4`. The results are
  identical; only the speed changes.

In RStudio, set both with **Knit ▸ Knit with Parameters…**. From a terminal,
edit `workers` at the top of `render_all.R`, or run:

```r
rmarkdown::render("02_fit_models.Rmd", params = list(workers = 4, refit = TRUE))
```

## Data and models

GEO reports RMA expression on a linear scale. The starting dataset has 984 arrays
and 21,727 probes. Nine paired donors have inconsistent recorded ages and are
excluded, so the analysis includes 379 paired donors and 817 arrays. Probes with
nonpositive values are removed before the log2 transformation, which leaves
21,462 probes. No additional normalisation or batch correction is applied.

| Model | Formula |
|---|---|
| Combined (both cell types) | `~ Age + (1 | Sex) + (1 | CellType) + (1 | Batch) + (1 | Individual)` |
| Within each cell type | `~ Age + (1 | Sex) + (1 | Batch) + (1 | Individual)` |

Age is standardised and treated as a fixed effect. The other terms are random
effects. The models use maximum likelihood (`REML = FALSE`), and each probe is
fitted from scratch (`useInitialFit = FALSE`). Each probe's fractions add up to 1.

## Results at a glance

![Variance components across probes and examples](results/figures/submission_summary.png)

Across probes, cell type accounts for a median 49.0% of expression variation and
donor for 5.4%. TLR4 illustrates a strong cell-type component, KDM5D a
sex-associated component and GSTM1 a donor-associated component. These are
probe-level estimates, and donor-associated variation is not a heritability
estimate.

The GEO matrix shows very little batch-associated variation. That may reflect
earlier processing, but the full preprocessing history is unknown. Within-cell
donor estimates rely on a small number of repeat arrays, so interpret them
cautiously.

## Folder layout

```
01_prepare_data.Rmd … 04_check_results.Rmd   the analysis, in order
install_packages.R                           one-time package setup
render_all.R                                 knit all notebooks from a terminal
data/source_files.csv                        GEO download links + MD5 checksums
data/raw/, data/processed/                   downloads and prepared data (not in Git)
results/models/                              fitted models (not in Git)
results/tables/, results/figures/            per-probe results, summaries, figures
report/                                      one-page summary (Word, PDF, Markdown)
```

## Where each original script went

| Original script (`immvar-variancepartition/`) | Notebook |
|---|---|
| `install_packages.R` | `install_packages.R` (adds `rmarkdown` and `knitr`) |
| `scripts/01_prepare_data.R` | `01_prepare_data.Rmd` |
| `scripts/02_fit_models.R` | `02_fit_models.Rmd` |
| `scripts/03_make_figures.R` | `03_make_figures.Rmd` |
| `scripts/04_check_results.R` | `04_check_results.Rmd` |
| `run_analysis.R` | `render_all.R` |

The analysis code is unchanged, and the notebooks write the same output files
as the original scripts. The additions are explanations, previews of
intermediate data, the one-probe walkthrough, extra exploratory plots and
tables shown only in the notebooks, and the option to reuse saved fits.

## References

- Hoffman and Schadt (2016). [variancePartition](https://doi.org/10.1186/s12859-016-1323-z).
- Raj et al. (2014). [Original ImmVar study](https://doi.org/10.1126/science.1249547).
- De Jager et al. (2015). [ImmVar overview](https://doi.org/10.1016/j.smim.2015.03.003).
