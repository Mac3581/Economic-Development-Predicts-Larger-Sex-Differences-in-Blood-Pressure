# Economic Development Predicts Larger Sex Differences in Blood Pressure: Longitudinal Evidence from 150 Countries

## Overview

This repository contains the R script and data for "Economic Development Predicts Larger Sex Differences in Blood Pressure: Longitudinal Evidence from 150 Countries" (under review).

The primary analysis combines country-level, sex-specific SBP estimates with gross domestic product (GDP) per capita for 1980, 1990, 2000, and 2008. Supplementary analyses consider the Human Development Index (HDI), the Gender Inequality Index (GII), and lagged GDP per capita. A replication analysis uses blood pressure estimates from the NCD Risk Factor Collaboration (NCD-RisC) for 1975–2015.

The study uses longitudinal mixed-effects models, region-clustered robust inference, within–between (Mundlak) decomposition, and two-way fixed-effects regression. The unit of observation is the country–sex–year. 

## Data sources

Data access links:

- World Bank GDP per capita: <https://data.worldbank.org/indicator/NY.GDP.PCAP.KD>
- NCD-RisC blood pressure data: <https://ncdrisc.org/data-downloads-blood-pressure.html>
- UNDP Human Development Data Center (reference for HDI and GII definitions and data): <https://hdr.undp.org/data-center>

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
