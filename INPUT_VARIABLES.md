# Cleaned analysis data

Supply one R `data.frame` named `analysis_data`, with one row per participant. Retain the original complete-case sample, factor levels and pooled survey weights. Raw questionnaire files are not read by this package.

## Numeric variables

| Column | Meaning and unit |
| --- | --- |
| `SEQN` | Unique participant identifier |
| `age` | Age, years, used in baseline descriptions |
| `bmi` | BMI, kg/m² |
| `waist` | Waist circumference, cm |
| `water` | Previously derived daily plain-water intake, g/day |
| `uric_acid` | Serum uric acid, mg/dL |
| `PAQ` | Total PA, MET-min/day |
| `sedentary` | ST, min/day |
| `WTMEC2YR` | Original two-year examination weight |
| `WTMEC12YR` | Pooled examination weight, `WTMEC2YR / 6` |
| `SDMVPSU` | Masked survey PSU identifier |
| `SDMVSTRA` | Masked survey stratum identifier |
| `PAQ_work_vig` | Previously derived vigorous work/household activity, MET-min/week |
| `PAQ_work_mod` | Previously derived moderate work/household activity, MET-min/week |
| `PAQ_transport` | Previously derived transport activity, MET-min/week |
| `PAQ_leisure_vig` | Previously derived vigorous leisure activity, MET-min/week |
| `PAQ_leisure_mod` | Previously derived moderate leisure activity, MET-min/week |

The five supplied component columns are weekly values. Their sum divided by 7 must equal daily `PAQ`. Zero denotes recorded nonparticipation. The analysis checks these identities before fitting and does not reconstruct the raw activity questionnaires.

## Outcome and categorical covariates

`kidney_stone` is coded 0/1, where 1 denotes a self-reported history of KSD. Numeric 0/1 or a factor with labels `0`, `1` is accepted. All covariates below must be unordered R factors. Preserve the listed level order from the existing cleaned data.

| Column | Levels in order |
| --- | --- |
| `sex` | `Female`, `Male` |
| `age_group` | `20-40`, `40-60`, `60+` |
| `race` | `Mexican American`, `Non-Hispanic Black`, `Non-Hispanic White`, `Other Hispanic`, `Other Race` |
| `education` | `9-11th grade`, `College or above`, `High school graduate`, `Less than 9th grade`, `Some college` |
| `marital` | `Divorced`, `Living with partner`, `Married`, `Never married`, `Separated`, `Widowed` |
| `pir` | `<1`, `1-3.5`, `>=3.5` |
| `diabetes` | `0`, `1` |
| `hyperlipid` | `0`, `1` |
| `hypertension` | `0`, `1` |
| `drink` | `0`, `1` |
| `smoke` | `0`, `1` |
| `water_g` | `1`, `2`, `3`, `4` |
| `uric_g` | `1`, `2`, `3`, `4` |

The existing `age_group` labels correspond to ages 20-39, 40-59 and ≥60 years. They are preserved as stored, including the original label spelling. Age enters the association models as this factor, while BMI and waist circumference enter continuously. The original water categories use cutpoints 550, 1000 and 1650 g/day. Uric-acid categories use cutpoints 4.5, 5.4 and 6.3 mg/dL. Binary disease and lifestyle indicators must already follow the original project's definitions.

The analysis code creates `PA_group`, `ST_group`, `ASBI`, `ASBI_group`, `joint_group` and `exclude_P1S1` from the supplied exposure values. Existing versions of those columns are replaced only within the returned analysis data, ensuring the published definitions are used. Input rows are not removed, and the source data frame is not overwritten.
