# AGENTS.md

Guidance for working in this repository.

## What this is

A `{golem}` Shiny app for visualizing high-throughput proteomics data. The
current dataset is a public PRIDE study (PXD073413) comparing KRAS-mutant vs
wild-type colorectal cancer, but **the app is built to be reusable for other
datasets/studies.**

## Core principle: keep modules and functions generic

The single most important convention: **functions in `R/fct_*.R` and modules in
`R/mod_*.R` must NOT hardcode study-specific information** (gene names, group
labels, condition names, dataset titles, study-specific colors, search terms,
etc.).

- Study-specific values are **parameters** with sensible generic defaults.
- The **app layer is the composition root** — `R/app_ui.R` and `R/app_server.R`
  are where study-specific values are wired in. Hardcoding KRAS / "Mutant" /
  "Wildtype" etc. is acceptable *only* here (and in `data_raw/process_data.R`
  and `inst/app/about/*.md`).

If you add a feature that needs a study-specific string/value, add a parameter
to the function/module and pass the concrete value from `app_server.R` /
`app_ui.R`.

### Examples of this pattern already in place

- `plot_volcano()` — `up_label` / `down_label` / `ns_label` (default
  "Up"/"Down"/"NS"); app passes "KRAS mutant enriched" / "KRAS wild-type
  enriched". Colors are also params.
- `mod_volcano_server()` — `up_label` / `down_label` passed through to
  `plot_volcano()`.
- `mod_protein_detail_server()` — `group_colors` (NULL = default palette) and
  `pubmed_context` (NULL = search gene alone; app passes "KRAS").
- `mod_intro_ui()` — all content via a `content` list (+ markdown files in
  `inst/app/about/`); `groups_value` passed to the server.
- Column names are always params (`fc_col`, `p_col`, `feature_col`,
  `accession_col`) defaulting to this dataset's columns.

## Data

- `data_raw/process_data.R` builds the app data from the raw PRIDE files.
  - `set.seed(42)` before `impute_left_dist()` — imputation is stochastic, so
    the seed keeps results reproducible across rebuilds. Keep it.
  - Outputs: `inst/app/data/cancer_de.rds` (full limma DE table incl. per-sample
    intensity columns named `M#####a`), `inst/app/data/sample_meta.rds`
    (SampleID -> group), and `inst/app/data/cancer_se.rds` (the filtered/imputed
    SummarizedExperiment used by the heatmap).
- SummarizedExperiment layout: assay `intensities`; rowData has the feature
  annotation incl. `Genes` (may be ";"-separated); colData has `mutated`
  (mutant/Wildtype) and `mutation` (fine-grained KRAS variant). The heatmap
  matches genes against a **rowData** column and groups samples by **colData**
  column(s) -- averaging within each group, then z-scoring each row.
- In the DE intensity columns, `0` means undetected/missing. The boxplot in
  `mod_protein_detail` drops zeros and shows per-group detection counts, because
  limma's `logFC` is computed on **imputed** data and can disagree with
  detected-only intensities for high-missingness proteins. Preserve that
  honesty if you touch the boxplot.

## Running / verifying

- Dev: restart R, then source `dev/run_dev.R` (it reloads the package before
  `run_app()`). Plain file edits are NOT picked up by a running app.
- The volcano uses `plotly::toWebGL()` (client-side GPU rendering) and
  `event_register(p, "plotly_click")` so click events fire reliably. Clicks are
  surfaced via `mod_volcano_server()`'s returned `clicked_gene` reactive.

## Provenance

Dataset is CC0 (public domain), PRIDE PXD073413 (Wu H. et al., Soochow
University). The About tab carries an independence disclaimer — this is a
portfolio/demo reanalysis, not affiliated with the original authors.
