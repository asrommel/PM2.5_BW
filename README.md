# Air Pollution and Birthweight in New York City

Code for analyses examining prenatal PM2.5 exposure and birthweight-for-gestational-age z-scores in the Generation C cohort, including effect modification by neighborhood-level structural disadvantage.

## Project structure
- `scripts/`: analysis scripts in execution order
- `R/`: helper functions
- `data/raw/`: restricted raw input data (not shared)
- `data/derived/`: derived analytic files (not shared)
- `outputs/`: generated tables and figures

## Reproducibility
This repository uses `renv` for package version management.

```r
install.packages("renv")
renv::restore()
