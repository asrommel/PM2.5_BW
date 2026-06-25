# Air Pollution and Birthweight in New York City

Code for analyses examining prenatal PM2.5 exposure and birthweight-for-gestational-age z-scores in the Generation C cohort, including effect modification by neighborhood-level structural disadvantage.

## Citation

Tubassum R, Li M, Boychuk N,  Milazzo F, Lee A, Janevic T, Rommel A-S. (_submitted_). Prenatal PM2.5 Exposure and Birthweight in the Context of Neighborhood-Level Structural Disadvantage: A Cohort Study in New York City. _Environmental Research_

---

## Data Availability

**Raw participant-level data are not publicly available** due to data protection requirements (Generation C cohort, IRB-20–00425). Researchers interested in access may contact anna.rommel@mssm.edu.

This repository contains the complete analysis code. To reproduce analyses with authorized access to the data, follow the instructions below.

---

## Project Structure

```
PM2.5_BW/
├── scripts/              # Analysis scripts (execute in numbered order)
├── R/                    # Helper functions
├── data/
│   ├── raw/             # Restricted raw input data (not shared)
│   └── derived/         # Derived analytic files (not shared)
├── outputs/             # Generated tables and figures
├── renv.lock            # R environment snapshot (for reproducibility)
├── .gitignore           # Excludes data and sensitive files
└── README.md
```

---

## System Requirements

- **R version:** 4.1.0 or later
- **Package management:** [renv](https://rstudio.github.io/renv/) (recommended)

### Key Dependencies

- `bdlim` – distributed lag interaction models
- `growthstandards` – INTERGROWTH-21st z-score calculations  
- `mice` – multiple imputation
- Data handling: `dplyr`, `tidyr`, `readr`, `fst`
- Visualization: `ggplot2`, `flextable`, `gt`

### Installing Required Packages

1. **Install from CRAN:**
```r
   install.packages(c("dplyr", "tidyr", "readr", "lubridate", "fst", 
                      "stringr", "mice", "ggplot2", "janitor", "forcats", 
                      "zoo", "scales", "flextable", "officer", "ggtext", "gt"))
```

2. **Install from GitHub:**
```r
   remotes::install_github("AnderWilson/bdlim")
   remotes::install_github("ki-tools/growthstandards")
```

3. **Or use renv to restore the exact environment:**
```r
   renv::restore()
```

---

## Running the Analysis

### Prerequisites
- Authorized access to Generation C cohort data
- Input data files must be placed in `data/raw/` with expected names:
  - `predictions_2025-10-10.fst` – PM2.5 predictions
  - `merged_address_data_2025-10-22.csv` – Address history
  - `birth_data_2025-10-14.csv` – Birth outcomes & covariates
  - `hypertensive_disorders_2024-10-16.csv` – Pregnancy complications
  - `mean_srei_genc.csv` – Structural racism effect index

### Execution
```r
source("scripts/00_setup.R")      # Load libraries, set paths
source("scripts/01_data_processing.R")  # Prepare analytic dataset
source("scripts/02_analysis.R")    # Primary BDLIM analyses
source("scripts/03_sensitivity.R") # Sensitivity analyses
source("scripts/04_tables_figures.R") # Generate outputs
```

**Estimated runtime:** ~1–2 hours (BDLIM sampling is intensive)

---

## Analysis Overview

### Study Population
- Generation C cohort (New York City-based birth cohort)
- Term births (≥259 days gestation) with complete PM2.5 exposure data
- n = 2,194

### Exposure
- Prenatal PM2.5 (µg/m³) from a satellite-based prediction model, XGBoost-IDW Synthesis (XIS-PM2.5, Just et al., 2025) during weeks 1–37 of gestation

### Outcome
- Birthweight-for-gestational-age z-scores (INTERGROWTH-21st standard)

### Key Analysis
- **Primary:** Bayesian distributed lag interaction model (BDLIM) with effect modification by neighborhood structural disadvantage (SREI)
- **Adjusted for:** Pre-pregnancy BMI, maternal age, parity, marital status, insurance, season of conception, year of birth
- **Sensitivity:** Pregnancy-average PM2.5 models; analyses across multiple imputed datasets

### Output Files
- `outputs/tables/Table1_demographics.docx` – Cohort characteristics
- `outputs/tables/Supplementary_Table_Imputed_Results_Full.docx` – Results across imputations
- `outputs/figures/Flowchart_exclusions.pdf` – Participant selection flowchart
- BDLIM summary plots (per model)

---

## Code Quality & Reproducibility

- **Random seed:** Set to 500 for consistency
- **Multiple imputation:** 5 datasets generated via MICE; primary analysis uses draw 1
- **Environment snapshot:** `renv.lock` captures exact package versions
- **Intermediate outputs:** `data/derived/` contains processed datasets for troubleshooting

---

## Modifications for Local Use

If adapting this code for similar analyses:

1. Update data input paths in `scripts/00_setup.R`
2. Modify covariate selection in BDLIM model call
3. Adjust outcome definition (z-scores vs. raw birthweight)
4. Update visualization color schemes in `ggplot2` calls as needed

---

## Contact & Contributing

**Questions about the analysis:** anna.rommel@mssm.edu

**Data access requests:** Contact anna.rommel@mssm.edu

---

## License

CC-BY-4.0

---

## Acknowledgments

**Funding**
This work was supported by the Simons Foundation [PI Rommel: 866027] and the National Institute of Child Health and Human Development (NICHD) [PI Rommel: R01HD109613]. Additional support was provided through the computational and data resources and staff expertise of Scientific Computing and Data at the Icahn School of Medicine at Mount Sinai, made possible by the Clinical and Translational Science Award (CTSA) grant [UL1TR004419] from the National Center for Advancing Translational Sciences. The funders had no role in design, analysis, or the decision to publish. The findings and conclusions in this report are those of the authors and do not necessarily represent the position of the funding agencies. 

**Acknowledgements**
We are deeply grateful to the Generation C participants who generously contributed their time and effort to this research. We also thank Allan C. Just, who led the development of the PM2.5 exposure model used in this study, and Kodi B. Arfer for his invaluable assistance in providing, preparing, and implementing the exposure data, as well as for his technical guidance throughout the project. 
