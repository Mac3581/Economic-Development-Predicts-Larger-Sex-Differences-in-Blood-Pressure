# Economic Development Predicts Larger Sex Differences in Blood Pressure: Longitudinal Evidence from 150 Countries

## Overview

This repository contains an R analysis pipeline investigating the longitudinal association between socioeconomic development and population-level mean systolic blood pressure (SBP), with particular attention to sex differences and their evolution over time.

The primary analysis combines country-level, sex-specific SBP estimates with gross domestic product (GDP) per capita for 1980, 1990, 2000, and 2008. Supplementary analyses consider the Human Development Index (HDI), the Gender Inequality Index (GII), and lagged GDP per capita. A replication analysis uses blood pressure estimates from the NCD Risk Factor Collaboration (NCD-RisC) for 1975–2015.

The study uses longitudinal mixed-effects models, region-clustered robust inference, within–between (Mundlak) decomposition, and two-way fixed-effects regression. The unit of observation is the country–sex–year. The analyses are ecological and observational; estimated associations should not be interpreted as individual-level or causal effects.

## Research objectives

- Examine the association between GDP per capita and mean SBP across countries and over time.
- Investigate whether this association differs between males and females and changes across historical periods using three-way interaction models.
- Distinguish within-country changes in economic development from between-country differences in economic development.
- Assess sensitivity to temporal autocorrelation, regional dependence, country and year fixed effects, and alternative exposure timings.
- Explore whether sex differences in SBP vary with human development and gender equality, including patterns relevant to the gender equality paradox.
- Assess the reproducibility of the main findings using an alternative international SBP dataset.

## Data sources

The script expects the following CSV files in the R working directory. File names are reproduced exactly as they appear in the code.

| Input file | Description | Years used |
| --- | --- | --- |
| `WebTable 4 appendix.csv` | Primary country- and sex-specific mean SBP estimates and reported uncertainty bounds from the original study's supplementary Web Table 4. | 1980, 1990, 2000, 2008 |
| `API_NY.GDP.PCAP.KD_DS2_en_csv_v2_89 (1).csv` | World Bank GDP per capita, constant US dollars (indicator `NY.GDP.PCAP.KD`). | 1980–2008 for the main analysis; additional historical years for lagged and replication analyses |
| `human-development-index (1).csv` | Human Development Index (HDI) data by country and year. | 1990, 2000, 2008 |
| `GII complete.csv` | Gender Inequality Index (GII) data in wide format. | 1990, 2000, 2008 |
| `NCD_RisC_Lancet_2017_BP_age_standardised_countries.csv` | NCD-RisC age-standardized, sex-specific country-level SBP estimates and 95% uncertainty intervals. | 1975, 1985, 1995, 2005, 2015 |

Data access links:

- World Bank GDP per capita: <https://data.worldbank.org/indicator/NY.GDP.PCAP.KD>
- NCD-RisC blood pressure data: <https://ncdrisc.org/data-downloads-blood-pressure.html>
- UNDP Human Development Data Center (reference for HDI and GII definitions and data): <https://hdr.undp.org/data-center>

The exact publication associated with Web Table 4, the provenance and release versions of the HDI and GII CSV files, and the GDP download vintage should be documented alongside the archived input files before public release. The UNDP link is a reference resource, not confirmation of the exact CSV exports used in the script.

## Software requirements

The analysis is written in R and requires the following packages:

| Purpose | Packages |
| --- | --- |
| Data import, manipulation, and reshaping | `readr`, `dplyr`, `tidyr`, `stringr`, `purrr` |
| Visualization | `ggplot2`, `patchwork` |
| Mixed-effects and fixed-effects models | `lme4`, `lmerTest`, `nlme`, `fixest` |
| Robust statistical inference | `clubSandwich` |
| Descriptive statistics and model diagnostics | `DescTools`, `performance`, `insight`, `car` |
| Country and region classification | `countrycode` |
| Spatial dependence diagnostics | `spdep` |

Install the dependencies once in an R session:

```r
install.packages(c(
  "readr", "dplyr", "tidyr", "stringr", "purrr",
  "ggplot2", "patchwork", "lme4", "lmerTest", "nlme",
  "fixest", "clubSandwich", "DescTools", "performance",
  "insight", "car", "countrycode", "spdep"
))
```

No specific R or package versions are pinned in the supplied code. Recording `sessionInfo()` and using a dependency lockfile (e.g., with `renv`) is recommended for reproducibility.

## Running the analysis

A suggested repository layout is:

```text
project/
├── README.md
├── analysis.R
├── WebTable 4 appendix.csv
├── API_NY.GDP.PCAP.KD_DS2_en_csv_v2_89 (1).csv
├── human-development-index (1).csv
├── GII complete.csv
└── NCD_RisC_Lancet_2017_BP_age_standardised_countries.csv
```

Here, `analysis.R` is a suggested name for the supplied R script; substitute its actual filename if different.

1. Save the analysis script and the five CSV inputs in the same directory.
2. Install the required R packages.
3. Open the project directory as the working directory in R or RStudio.
4. Run the script sequentially, or execute it with:

```r
source("analysis.R")
```

For interactive inspection, run the script section by section in RStudio. Several sections use `View()` and display plots interactively. Model fitting and multi-optimizer checks may be computationally intensive. The supplied script does not automatically export results, figures, or fitted model objects.

## Data preparation

The primary SBP data are imported from a wide-format supplementary table and reshaped into a long-format dataset. Mean SBP and the accompanying lower and upper uncertainty bounds are extracted from text fields. GDP data are similarly reshaped into country–year format and merged with SBP estimates after harmonizing country names.

The primary sample is restricted to countries with available mean SBP and log GDP per capita in 1980. Subsequent waves are retained for these countries, subject to data availability. This baseline eligibility rule does not require a fully balanced panel. The script excludes Cook Islands, the Democratic People's Republic of Korea, and Taiwan from the primary and lagged GDP merges.

Countries are assigned to United Nations geographical groupings using `countrycode`. Intermediate regions are used where available, with subregions used otherwise. This classification provides the clusters for robust inference and the regional proximity structure for Moran's I. The exact number of represented regions depends on the analytical sample.

### Key variables

| Variable | Definition |
| --- | --- |
| `Country` | Country identifier and grouping factor. |
| `Year` | Calendar year of the SBP estimate. |
| `Mean_SBP` / `mean_SBP` | Country- and sex-specific mean systolic blood pressure (mmHg). |
| `Lower_CI`, `Upper_CI` | Reported SBP uncertainty bounds, retained but not used as model weights. |
| `Sex` | Effect-coded as Female = −0.5 and Male = +0.5. |
| `Year_c` | Years since the analysis-specific reference year (1980, 1990, or 1975). |
| `GDP_per_capita` | GDP per capita in constant US dollars. |
| `log_GDP` | Natural logarithm of GDP per capita. |
| `log_GDP_z` | Standardized log GDP per capita (z-score). |
| `UN_regions_pure22` | Combined UN intermediate-region/subregion classification. |
| `gdp_within` | Deviation of standardized log GDP from the country's observed-period mean. |
| `gdp_between` | Country-specific observed-period mean of standardized log GDP. |
| `HDI_z` | Standardized Human Development Index. |
| `gii_inverted` | Reversed Gender Inequality Index, calculated as `1 - gii`. Higher values indicate lower gender inequality. |
| `gii_z` | Standardized reversed GII. Higher values represent greater gender equality as operationalized by this index. |
| `log_GDP_5lag`, `log_GDP_10lag` | Natural log GDP per capita measured five and ten years before the corresponding SBP observation. |

Effect coding makes the intercept in the interaction models the midpoint of the predicted female and male SBP values at the reference year and zero on the standardized predictor. The sex coefficient represents the male–female contrast at those reference values. Standardized log GDP equal to zero corresponds to the mean of the log-transformed GDP values, not the arithmetic mean of GDP in dollars.

## Statistical analysis

### 1. Primary longitudinal mixed-effects models (1980–2008)

The primary models are fitted with `lme4::lmer()` using maximum likelihood (ML) and the `bobyqa` optimizer.

| Model | Fixed effects | Country-level random effects |
| --- | --- | --- |
| `M0` | Intercept only | Intercept |
| `M1` | Sex + standardized log GDP + centered year | Intercept and slopes for sex and year |
| `M2` | Sex × standardized log GDP × centered year, including all lower-order terms | Intercept and slopes for sex and year |

The full interaction model is specified as:

```r
Mean_SBP ~ Sex * log_GDP_z * Year_c +
  (1 + Sex + Year_c | Country)
```

The sex × GDP interaction tests whether the GDP–SBP association differs by sex at the reference year. The three-way interaction assesses whether that sex difference in the GDP–SBP association changes linearly over time. Equivalently, it assesses whether the relationship between GDP and the male–female SBP gap changes with time.

The script computes variance components, intraclass correlation coefficients (ICCs), likelihood-based model comparisons, variance inflation factors (VIFs), and optimizer sensitivity checks where specified.

### 2. Region-clustered robust inference

Cluster-robust variance–covariance matrices are calculated with `clubSandwich::vcovCR()` using the CR2 adjustment and the UN geographical grouping as the clustering variable. Coefficient tests are obtained with `coef_test()`, which uses Satterthwaite degrees of freedom by default. This approach addresses possible dependence among observations within the same geographical cluster while retaining the mixed-effects point estimates.

The geographical clustering is intended to accommodate potential regional dependence. It does not directly model geographic distance or eliminate unmeasured regional confounding. The ordinary likelihood-ratio comparisons in `anova()` are not cluster-robust tests.

### 3. Model-based predictions

Fixed-effect predictions are generated over a grid of standardized development values, separately for females and males and for each observation year. Prediction standard errors are computed using the relevant robust variance–covariance matrix:

`SE = sqrt(diag(X V Xᵀ))`

where `X` is the fixed-effect design matrix and `V` is the coefficient covariance matrix. The plots include pointwise 95% confidence bands calculated as the fitted value ± 1.96 standard errors.

These bands use a normal approximation, rather than the small-sample degrees-of-freedom adjustments used for coefficient tests. Predictions exclude country-specific random effects and describe fitted population-level mean relationships, not prediction intervals for individual countries or people.

### 4. Temporal autocorrelation

Sensitivity analyses use `nlme::lme()` with first-order autoregressive residual correlation (AR(1)). The combined models specify residual dependence within country–sex series; separate male and female models specify it within countries. Random intercepts and year slopes are included in the sex-specific models.

### 5. Regional autocorrelation of residuals

Global Moran's I is calculated on country-averaged conditional residuals from the full mixed-effects model. Countries are defined as neighbors when they share the same UN geographical cluster, and the resulting weights matrix is row-standardized.

This diagnostic assesses residual similarity under a regional-membership proximity definition. It is not a distance-based or border-contiguity spatial analysis.

### 6. Within–between (Mundlak) decomposition of GDP

To distinguish longitudinal within-country associations from cross-country associations, standardized log GDP is decomposed using one record per country–year:

- **Between-country component:** the country-specific mean of standardized log GDP across observed periods.
- **Within-country component:** the deviation from that country's mean at a given observation year.

The full within–between interaction model is:

```r
Mean_SBP ~ Sex * Year_c * (gdp_within + gdp_between) +
  (1 + Sex + Year_c | Country)
```

The model includes separate sex and time interactions for the within- and between-country GDP components. The decomposition retains the original standardized log GDP scale rather than independently standardizing the two components.

The script estimates null, main-effects, and full-interaction versions of the model, applies CR2 inference, and checks optimizer agreement. For the primary dataset, separate figures show within-country and between-country associations. The within-country predictions are restricted to the marginal observed range in each year, with the between-country component held at zero. Between-country predictions hold the within-country component at zero. These conditional predictions do not necessarily correspond to jointly observed predictor combinations.

### 7. Two-way fixed-effects regression

The fixed-effects sensitivity analysis is estimated with `fixest::feols()`:

```r
Mean_SBP ~ Sex * log_GDP_z | Country + Year_c
```

Country and year fixed effects absorb time-invariant country-specific factors and common year-specific shifts in SBP, respectively. Standard errors are clustered by UN geographical region.

The primary fixed-effects figure displays the fitted contribution of the sex, GDP, and interaction terms, excluding the absorbed country and year effects. It is not an absolute predicted SBP level. The fixed-effects specification does not include the three-way interaction with time and therefore addresses a related, but different, estimand.

## Supplementary analyses

### Human Development Index (HDI)

The primary SBP data are combined with HDI values for 1990, 2000, and 2008. HDI is standardized and year is centered at 1990. Null, main-effects, and full-interaction mixed-effects models are estimated with country-level random intercepts. The full model is:

```r
Mean_SBP ~ Sex * HDI_z * Year_c + (1 | Country)
```

The script produces CR2 coefficient tests and predicted SBP curves by HDI, sex, and year. Unlike the primary GDP model, the HDI specification does not include country-level random slopes.

### Lagged GDP per capita

GDP per capita measured five and ten years before each SBP observation is used as an alternative economic exposure. For example, the SBP observation in 1980 is matched to GDP in 1975 and 1970, respectively. Countries must have both lagged GDP measures and an SBP estimate available at the 1980 baseline to enter the lagged analysis sample.

Separate mixed-effects models (`M25` and `M210`) include the three-way interaction between sex, standardized lagged log GDP, and centered year, as well as country-level random intercepts and slopes for sex and year. Additional two-way fixed-effects models are estimated for both lag lengths.

These analyses examine sensitivity to exposure timing. Lagging GDP does not, by itself, establish causality or eliminate reverse causation and confounding.

### Gender Inequality Index (GII)

GII data are harmonized and merged with the primary SBP sample for 1990, 2000, and 2008. The index is reversed (`1 - GII`) and standardized so that higher scores correspond to lower gender inequality.

The full random-intercept model is:

```r
Mean_SBP ~ Sex * gii_z * Year_c.y + (1 | Country)
```

Here, `Year_c.y` is years since 1990, resulting from the merge of the SBP and GII datasets. The model investigates whether the sex-specific association between gender equality and mean SBP varies over time. The analysis includes null and main-effects models, CR2 inference, optimizer checks, and predicted SBP curves.

### Correlations among development indicators

The script computes Pearson correlations among GDP per capita (untransformed), HDI, and reversed GII using one selected record per country and pairwise complete observations. It also reports the number of available observations for each correlation. As currently implemented, the records for the three indicators may come from different years; the resulting correlations should not be described as contemporaneous associations unless the merge is revised to align calendar years.

## Replication analysis: NCD-RisC (1975–2015)

The replication analysis uses age-standardized SBP estimates from NCD-RisC at five ten-year intervals (1975, 1985, 1995, 2005, and 2015). Countries are retained only if SBP and GDP data are available for both sexes at all five selected years. Year is centered at 1975, and log GDP is standardized within the replication dataset.

The following analyses are repeated using the alternative SBP estimates:

- Null, main-effects, and three-way interaction mixed-effects models (`M0A`, `M1A`, `M2A`).
- CR2 region-clustered robust inference and model diagnostics.
- Combined and sex-specific AR(1) sensitivity models.
- Global Moran's I using regional membership weights.
- Within–between GDP decomposition and interaction models.
- Two-way country and year fixed-effects regression (`mod_fe_basencd`).

This is a replication using an alternative SBP estimation series, time window, and sample-selection rule. It should not automatically be characterized as fully independent, because underlying source studies may overlap between the primary and NCD-RisC datasets.

## Outputs

The script generates statistical results in the R console and visualizations in the active graphics device.

| Output | Description |
| --- | --- |
| Descriptive summaries | Sample composition, missing GDP observations, country and region counts, and variable distributions. |
| Model results | Fixed-effect estimates, random-effect estimates, robust standard errors, coefficient tests, and selected confidence intervals. |
| Model diagnostics | ICCs, variance components, VIFs, optimizer checks, likelihood comparisons, and Moran's I. |
| Primary GDP figure | Predicted mean SBP by standardized log GDP, sex, and year. |
| Within–between GDP figures | Separate conditional predictions for within-country and between-country GDP components. |
| Fixed-effects figure | Fitted sex- and GDP-related SBP component, excluding absorbed fixed effects. |
| HDI and GII figures | Sex- and year-specific predicted mean SBP by the respective standardized indicator. |
| Lagged GDP figures | Separate and combined plots for five- and ten-year lagged GDP associations. |
| Correlation results | Pearson correlation matrix and pairwise sample sizes for GDP, HDI, and reversed GII. |
| Replication results | NCD-RisC model estimates and robustness diagnostics. |

No numerical results or substantive empirical conclusions are reported in this README because the analysis outputs have not been supplied or reproduced here.

## Reproducibility and interpretation notes

The following checks are recommended before running the complete pipeline or using its results in a publication. They refer to the supplied script, which has not been modified by this README.

- **Country harmonization:** verify all manual recodes and unmatched countries, preferably using validated ISO codes. In particular, the script maps the U.S. Virgin Islands to the British Virgin Islands and Sint Maarten to the Netherlands Antilles; these are not interchangeable geographical entities. Check country-name consistency across GDP, HDI, GII, and NCD-RisC inputs, and enforce unique country–year or country–sex–year keys as appropriate.
- **Comparable model samples:** baseline eligibility does not guarantee complete observations at later years. Use identical complete-case samples for nested model comparisons, verify positive GDP before logarithmic transformation, and confirm that region identifiers are nonmissing in every fitted sample. The HDI full join may also introduce HDI-only countries into the standardization sample.
- **AR(1) implementation:** replace `lmerControl()` with `nlme::lmeControl()` in calls to `nlme::lme()` and specify missing-data handling consistently. The primary observation intervals are unequal (10, 10, and 8 years); clarify whether the correlation parameter should represent one calendar year or one measurement wave. Consider a continuous-time AR(1) structure if that matches the intended specification.
- **Ten-year lag predictions:** the current plotting block uses a GDP range standardized on a different dataset, labels facets with `Year_c + 1970` rather than the SBP observation year (`Year_c + 1980`), and does not explicitly align design-matrix and covariance-matrix names with the fitted coefficients. Correct these points before interpreting or exporting the ten-year lag figure.
- **Indicator labels and correlations:** label `gii_z` as standardized reversed GII (higher values indicate greater equality), and align GDP, HDI, and GII by a common year or a clearly defined country-level averaging period for correlation analysis.
- **Model stability and inference:** inspect convergence and singular-fit warnings, effective degrees of freedom, cluster sizes, and the appropriateness of random-effects structures. The small number of geographical clusters may limit precision. Consider small-sample-adjusted confidence intervals for plotted contrasts or predictions where appropriate.
- **Time specification and prediction support:** the mixed models assume a linear time trend and linear moderation by time. Some prediction grids span the overall observed exposure range rather than year-specific ranges, potentially extrapolating beyond observed values. Within-country GDP deviations may also be highly correlated with year.
- **Uncertainty and weighting:** the script retains the uncertainty bounds of the published SBP estimates but does not propagate their estimation uncertainty through the regressions or use population or precision weights. The fitted relationships describe country-level observations, not population-weighted global estimates.
- **Input schemas and environment:** several HDI, GII, and NCD-RisC preparation steps select columns by numeric position. Validate the CSV schemas and prefer named columns for future versions. Record the R and package versions, source data versions, and final analysis-specific sample sizes.

## Scope and limitations

The analysis evaluates associations between country-level socioeconomic indicators and country-level mean SBP. The observational and ecological design is subject to residual confounding, differences in data coverage, potential measurement error, and limits to generalizability. The within–between decomposition and fixed-effects models provide complementary perspectives on the data but do not establish causal effects. The HDI and GII analyses use different periods, available samples, and random-effects structures from the primary GDP analysis; comparisons across these models should account for those differences.

## Citation and licensing

Please cite the original SBP publication, the World Bank World Development Indicators, the exact HDI and GII data releases, and the NCD-RisC publication associated with the downloaded blood pressure estimates when using this analysis. The original study citation, repository authorship, code license, and data redistribution terms were not specified in the supplied script and should be added by the repository maintainer before public distribution.
