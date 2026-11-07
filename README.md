# ImmVar variancePartition

This project applies `variancePartition` to gene expression in naïve CD4 T cells and monocytes from ImmVar (GEO GSE56035). The question is how much expression variation is associated with cell type, donor, age, sex and batch.

The one-page summary is in [`report/ImmVar_summary.pdf`](report/ImmVar_summary.pdf) and [`report/ImmVar_summary.docx`](report/ImmVar_summary.docx).

## Run the analysis

Open `immvar.Rproj` in RStudio, or run these commands from this folder:

```sh
Rscript install_packages.R
Rscript run_analysis.R
```

The first script downloads the GEO files. Model fitting takes the most time. To run one model after preparing the data, use `Rscript scripts/02_fit_models.R combined`, `CD4_T` or `Monocyte`.

The scripts are numbered in the order they run:

1. `scripts/01_prepare_data.R` downloads the files, checks their checksums, prepares the expression matrix and sample information.
2. `scripts/02_fit_models.R` fits the combined and within-cell models.
3. `scripts/03_make_figures.R` writes the summary tables, figure and a short Markdown results summary.
4. `scripts/04_check_results.R` checks the fits and runs the sensitivity analyses.

## Data and models

GEO reports RMA expression on a linear scale. The starting dataset has 984 arrays and 21,727 probes. Nine paired donors have inconsistent recorded ages and are excluded; the analysis includes 379 paired donors and 817 arrays. Probes with nonpositive values are removed before log2 transformation, leaving 21,462 probes. No additional normalization or batch correction is applied.

The combined model is:

```r
~ Age + (1 | Sex) + (1 | CellType) + (1 | Batch) + (1 | Individual)
```

The models run separately within each cell type are:

```r
~ Age + (1 | Sex) + (1 | Batch) + (1 | Individual)
```

Age is standardized and treated as a fixed effect. The other terms estimate variance associated with their groups. The models use maximum likelihood (`REML = FALSE`). Each probe is fitted afresh with `useInitialFit = FALSE`. The reported fractions for each probe add up to one.

## Results

![Variance components across probes and examples](results/figures/submission_summary.png)

Across probes, cell type accounts for a median 49.0% of expression variation and donor for 5.4%. TLR4 illustrates a strong cell-type component, KDM5D a sex-associated component and GSTM1 a donor-associated component. These examples are probe-level estimates; donor-associated variation is not a heritability estimate.

The supplied GEO matrix has very little batch-associated variation. This may reflect earlier processing, but the full preprocessing history is unknown. Within-cell donor estimates use a small number of repeated samples, so they should be interpreted cautiously.

The per-probe variance tables and summary figure are included in `results/`. Raw GEO files, the processed expression matrix and fitted model objects are not included in Git. The scripts download the inputs and recreate them. The Word and PDF summary are provided in `report/`.

## References

- Hoffman and Schadt (2016). [variancePartition](https://doi.org/10.1186/s12859-016-1323-z).
- Raj et al. (2014). [Original ImmVar study](https://doi.org/10.1126/science.1249547).
- De Jager et al. (2015). [ImmVar overview](https://doi.org/10.1016/j.smim.2015.03.003).
