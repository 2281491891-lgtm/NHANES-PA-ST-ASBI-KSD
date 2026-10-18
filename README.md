# NHANES physical activity, sedentary time, ASBI and KSD

Core R code for the NHANES analyses used in the manuscript and supplementary information. The entry point is an already cleaned analysis data frame. Raw-file acquisition, merging, eligibility screening and missing-data processing are outside this code package.

## Files

| File | Purpose |
| --- | --- |
| `NHANES_manuscript_analysis.R` | Model functions, survey inference, standalone plots and optional Excel export |
| `run_analysis.R` | Run the included analyses on `analysis_data` and create named plot objects |
| `call_results.R` | Direct commands to display each table and each individual figure |
| `INPUT_VARIABLES.md` | Required input columns, units and factor coding |
| `package_versions.csv` | R and package versions used for verification |

## Run

Set the R working directory to this folder. The cleaned data frame must already exist in the R session as `analysis_data`. No participant-level data are distributed here.

```r
source("run_analysis.R", encoding = "UTF-8")
```

The original complete analysis stores the compatible cleaned data in `Final_results$data`. In an R session containing that object, assign it before running:

```r
analysis_data <- Final_results$data
source("run_analysis.R", encoding = "UTF-8")
```

To load only the functions, without fitting models:

```r
source("NHANES_manuscript_analysis.R", encoding = "UTF-8")
```

Required packages are `dplyr`, `survey`, `Hmisc` and `ggplot2`. `openxlsx` is needed only for the optional table export. Package installation is not performed automatically.

The code checks the fixed analysis sample of 19,977 participants, including 1,713 participants with a history of KSD. It fails if required variables are missing or the supplied sample is inconsistent. It does not silently exclude or impute participants. PA/ST categories, ASBI and component aggregates are defined within the analysis code. The supplied weekly PA component values are converted to daily values without rereading GPAQ files.

## Included analyses

- Survey-weighted baseline descriptions and group comparisons.
- Categorical PA, ST and ASBI logistic regression in Models 1-3.
- Sixteen joint categories of PA and ST, adjusted prevalence, categorical interaction and ST spline-by-PA interaction.
- Fully adjusted PA, ST and ASBI RCS models and exploratory piecewise slope tests.
- ASBI sensitivity curves for constants 0.01, 0.1 and 0.5, in the full sample and after excluding P1&S1.
- Three separate fully adjusted PA decompositions: intensity, domain, and domain-by-intensity. All components within each decomposition enter the same model. Weighted composition plots show their contributions within P1-P4.
- Six separate continuous-PA mediation analyses, including survey jackknife confidence intervals, detailed effect estimates and signed indirect-to-total ratios.

This package excludes the local cohort, simulated samples, LASSO/Boruta, response-only analyses, the minimum-location bootstrap, alternative PA references, additional PA-adjusted ASBI models and component RCS explorations. No screening flowchart is generated because the data-processing stage is outside the package.

## Analysis settings

The survey design uses `SDMVPSU`, `SDMVSTRA`, nesting, and the pooled examination weight `WTMEC12YR = WTMEC2YR / 6`. Model 1 is unadjusted. Model 2 includes sex, age group, race/ethnicity, education, marital status and poverty-to-income ratio. Model 3 additionally includes diabetes, hyperlipidemia, hypertension, BMI, waist circumference, alcohol, smoking, water-intake category and uric-acid category. PA and ST models additionally adjust for each other's category. ASBI models do not additionally adjust for PA or ST.

| Exposure | Categories | Reference |
| --- | --- | --- |
| PA, MET-min/day | P1 <20, P2 20 to <200, P3 200 to <600, P4 ≥600 | P1 |
| ST, min/day | S1 <240, S2 240 to <360, S3 360 to <480, S4 ≥480 | S1 |
| ASBI | Q1 <-2.14, Q2 -2.14 to <-0.47, Q3 -0.47 to <0.91, Q4 ≥0.91 | Q3 |
| Joint PA and ST | P1&S1 through P4&S4 | P4&S1 |

ASBI is `log((PAQ/60 + c)/(sedentary/60 + c))`, with `c = 0.1` in the main analysis. Survey-weighted quartiles (`qrule = "math"`) are rounded to two decimal places before classifying ASBI. The code recalculates and verifies the three manuscript cutpoints.

RCS models use four knots at the unweighted 5th, 35th, 65th and 95th percentiles. The model is survey weighted. The OR reference is the survey-weighted median, and plots display the 1st-99th percentile range. Each ASBI sensitivity sample uses its own weighted median. Its knot locations and plotting range are retained from the full sample for the same constant. Boundary minima are not annotated. Piecewise tests are conditional on the fitted minimum selected in the same sample.

The joint-category forest uses BH-adjusted q values for the 15 comparisons with P4&S1. Its confidence intervals are the covariate-adjusted, pointwise 95% model intervals. BH adjustment applies to the P values and does not turn these into simultaneous confidence intervals. The heatmap displays nominal P values.

The six mediation models use continuous PA and a fixed +100 MET-min/day contrast from its weighted median. Each model adjusts for demographic covariates, smoking, alcohol, water category and ST category. Other candidate mediators are not mutually adjusted. Both the mediator and outcome models are refitted in each stratified delete-one-PSU JKn replicate. This is deterministic jackknife inference and does not require a random seed. Signed indirect-to-total ratios describe the fitted decomposition and are not additive pathway fractions.

## Results and figures

`Manuscript_results` retains fitted models and detailed numeric results. `Main_tables` contains the two NHANES main tables. `Supplement_tables` contains variable definitions, threshold tests, joint-category estimates, PA component models and detailed mediation estimates. Descriptive table names are used because the full manuscript also contains local-cohort tables that are outside this package.

`Figure_index` lists 11 main and 16 supplementary standalone plot objects. Run individual commands in `call_results.R` to display them. The PA/ST category forest is combined as requested. All other plot objects remain individually callable. The existing Arial font and plot layouts are retained. No images are saved automatically and no figure captions are generated.

Optional export creates two Excel workbooks containing only the main and supplementary publication tables:

```r
Excel_files <- export_manuscript_tables(Manuscript_results, directory = "tables")
```

## Data source

The analysis concerns [NHANES 2007-2018](https://wwwn.cdc.gov/nchs/nhanes/Default.aspx). Reproduction from raw public files additionally requires the original data-preparation workflow, which is intentionally omitted from this core-code release.
