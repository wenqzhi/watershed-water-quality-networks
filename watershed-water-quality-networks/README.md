# WRR1 watershed water-quality network analysis

This repository contains the R scripts supplied for Figures 2–6. It is organized for version control and reproducible local use. The scripts retain their analysis logic; the machine-specific data path has been replaced with a single local setting.

## Repository contents

- `scripts/figures/`: figure and analysis scripts, renamed by figure and pipeline step.
- `config.R`: reads the local data-directory setting from `WRR_DATA_ROOT`.
- `docs/RUN_ORDER.md`: suggested execution order and dependencies.
- `.gitignore`: excludes local data, outputs, and geospatial/binary files from Git.

## Configure the data path

The data are not included in this repository. Set the `WRR_DATA_ROOT` environment variable to the existing local folder that contains the source data and intermediate directories used by the scripts. For example, in R before running a script:

```r
Sys.setenv(WRR_DATA_ROOT = "D:/WRR_project/data")
source("scripts/figures/fig02_network_structure.R")
```

`config.R` is sourced by each script, so run scripts with the repository root as the R working directory (or open `WRR1.Rproj` in RStudio first). On Windows, use forward slashes in the path. Do not commit raw data, restricted data, credentials, or API keys.

## Dependencies

Scripts load or install their own R packages. Depending on the figure, packages include `readr`, `dplyr`, `tidyr`, `stringr`, `tibble`, `ggplot2`, `sf`, `igraph`, `patchwork`, `scales`, `lubridate`, and `purrr`. `sf` may also require system geospatial libraries on some computers.

## Figure 4 versions

Two Figure 4 scripts are retained because they implement different comparisons:

- `fig04_yield_castdm.R`: newer yield-oriented CAST-DM comparison; it expects the additional CAST-DM annual export and may use a verified station-area table.
- `fig04_load_comparison_legacy.R`: earlier load-oriented comparison.

Use the version matching the manuscript and reported units. Do not treat the two outputs as interchangeable.

## Important

This is a code organization package, not a fully self-contained reproducibility archive. The scripts depend on local input data and, in several cases, on intermediate results produced by earlier analyses that are not among the supplied files. Review `docs/RUN_ORDER.md` and confirm each required input exists before running.

## Expanded analysis collection

The project now contains 27 R analysis scripts, including network construction, PCS computation, data quality/robustness analyses, and Figures S1–S6. See `docs/SOURCE_MAPPING.md` for the newly added filename mapping and `docs/RUN_ORDER.md` for dependency notes. Additional packages include those declared in the new scripts (e.g. `terra` for flow-direction raster processing).
