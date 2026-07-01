# Air Pollution and Birthweight in New York City

Code for analyses examining prenatal PM2.5 exposure and birthweight-for-gestational-age z-scores in the Generation C cohort, including effect modification by neighborhood-level structural disadvantage.

## Citation

Tubassum R, Li M, Boychuk N,  Milazzo F, Lee A, Janevic T, Rommel A-S. (_submitted_). Prenatal PM2.5 Exposure and Birthweight in the Context of Neighborhood-Level Structural Disadvantage: A Cohort Study in New York City. _Paediatric and Perinatal Epidemiology_

---

## Data Availability

**Raw participant-level data are not publicly available** due to data protection requirements (Generation C cohort, IRB-20–00425). Researchers interested in access may contact anna.rommel@mssm.edu.

This repository contains the complete analysis code. To reproduce analyses with authorized access to the data, follow the instructions below.

---

## Project Structure

```
PM2.5_BW/
├── scripts/
│   └── analysis.R        # Complete analysis: data build, primary models,
│                          # all sensitivity analyses, Table 1, flowchart
├── data/                  # Restricted raw input data (not shared; see below)
├── outputs/
│   ├── tables/            # Generated tables (Table 1, eTables 1-4)
│   └── figures/           # Generated figures (Figure 1 flowchart)
├── renv.lock              # R environment snapshot (for reproducibility)
├── .gitignore             # Excludes data and sensitive files
└── README.md
```

> Note: all helper functions, data processing, primary models, and sensitivity
> analyses are contained in the single script `scripts/analysis.R`, run
> top to bottom. There is no separate `R/` helper directory or numbered
> script sequence in the current version of this repository.

---

## System Requirements

- **R version:** 4.1.0 or later
- **Package management:** [renv](https://rstudio.github.io/renv/) (recommended)

### Key Dependencies

- `bdlim` – distributed lag interaction models
- `growthstandards` – INTERGROWTH-21st z-score calculations  
- `mice` – multiple imputation
- `here` – reproducible relative file paths
- Data handling: `dplyr`, `tidyr`, `readr`, `fst`
- Visualization: `ggplot2`, `flextable`, `gt`
- `EValue`, `tableone`, `WeightIt` – sensitivity analyses (E-values, selection bias comparison, inverse probability weighting)

### Installing Required Packages

1. **Install from CRAN:**
```r
   install.packages(c("dplyr", "tidyr", "readr", "lubridate", "fst", 
                   "stringr", "mice", "ggplot2", "janitor", "forcats", 
                   "zoo", "scales", "flextable", "officer", "ggtext", "gt",
                   "here", "EValue", "tableone", "WeightIt"))
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
- Input data files must be placed in `data/` with expected names:
  - `predictions_2025-10-10.fst` – PM2.5 predictions
  - `merged_address_data_2025-10-22.csv` – Address history
  - `birth_data_2025-10-14.csv` – Birth outcomes & covariates
  - `hypertensive_disorders_2024-10-16.csv` – Pregnancy complications
  - `mean_srei_genc.csv` – Structural racism effect index

### Execution

The complete analysis is contained in a single R script:

```r
source("scripts/analysis.R")
```

Or open `scripts/analysis.R` in RStudio and run it section by section. Sections are clearly marked with `########################################` comments.

**Estimated runtime:** several hours. The script fits 11 BDLIM models in total (2 primary, 5 across imputed datasets, 2 restricted to 2021–2022 births, 2 restricted to complete-case exposure series), most at 50,000 iterations. Actual runtime will depend on your hardware; time a full run once and update this estimate accordingly.

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

### Sensitivity Analyses
- Linear regression with pregnancy-average PM2.5, including a PM2.5 × SREI interaction term
- BDLIM restricted to 2021–2022 births (excluding pandemic-onset 2020)
- BDLIM restricted to participants with fully observed 37-week exposure series (no interpolated weeks)
- Comparison of included vs. excluded participants, with inverse probability weighting (IPW) to assess selection bias
- E-values to assess sensitivity to unmeasured confounding

### Output Files
- `outputs/tables/Table1_demographics.docx` – Cohort characteristics (Table 1)
- `outputs/figures/Flowchart_exclusions.pdf` / `.png` – Participant selection flowchart (Figure 1)
- `outputs/tables/Supplementary_Table1_Imputed.docx` – Primary BDLIM results across the 5 multiply imputed datasets (eTable 1)
- `outputs/tables/Supplementary_Table2_2021_2022.docx` – Sensitivity analysis restricted to 2021–2022 births (eTable 2)
- `outputs/tables/Supplementary_Table3_CompleteCases.docx` – Sensitivity analysis restricted to fully observed 37-week exposure series (eTable 3)
- `outputs/tables/Supplementary_Table4_SelectionBias.csv` – Included vs. excluded participant comparison with standardised mean differences (eTable 4)
- BDLIM summary plots (rendered to console/plot device for each model; not saved as files unless you add `ggsave()` calls)

---

## Code Quality & Reproducibility

- **Random seed:** Set to 500 for data processing/imputation and 123 for BDLIM model fitting, for consistency
- **Multiple imputation:** 5 datasets generated via MICE; primary analysis uses draw 1, with a full sensitivity check across all 5 draws (eTable 1)
- **Environment snapshot:** `renv.lock` captures exact package versions. Run `renv::snapshot()` after installing/updating any packages the script uses, and commit the updated `renv.lock`, so `renv::restore()` stays accurate for anyone else running this code.

---

## Modifications for Local Use

If adapting this code for similar analyses:

1. Update the `dir_data`, `dir_outputs`, `dir_tables`, `dir_figures`, and `path_*` variables near the top of `scripts/analysis.R` to point to your own data locations and filenames
2. Modify covariate selection in the BDLIM model calls
3. Adjust outcome definition (z-scores vs. raw birthweight)
4. Update visualization color schemes and labels in the `ggplot2`/`flextable`/`gt` calls as needed

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
