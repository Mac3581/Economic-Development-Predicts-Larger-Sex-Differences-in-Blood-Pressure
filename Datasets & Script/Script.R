# ---------------------------
# Data preparation
# ---------------------------

# Packages
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(nlme)
library(fixest)#for fixed effect regression
library(purrr)
library(DescTools)
library(performance)
library(clubSandwich)
library(countrycode)
library(patchwork)

# Import and reshape systolic blood pressure data

raw_sbp <- read_csv("WebTable 4 appendix.csv", 
                                col_names = FALSE,
                                show_col_types = FALSE)


######
sbp_analysis <- raw_sbp %>%
  select(1:9) %>%
  rename(
    Country = X1,
    male_1980 = X2,
    male_1990 = X3,
    male_2000 = X4,
    male_2008 = X5,
    female_1980 = X6,
    female_1990 = X7,
    female_2000 = X8,
    female_2008 = X9
  ) %>%
  mutate(
    Region = if_else(
      is.na(Country) & !is.na(male_1980),
      male_1980,
      NA_character_
    )
  ) %>%
  fill(Region, .direction = "down") %>%
  filter(
    !is.na(Country),
    Country != "Country",
    !is.na(male_1980)
  ) %>%
  pivot_longer(
    cols = matches("^(male|female)_\\d{4}$"),
    names_to = c("Sex", "Year"),
    names_pattern = "(male|female)_(\\d{4})",
    values_to = "SBP_values"
  ) %>%
  mutate(
    Sex = recode(Sex, male = "Male", female = "Female"),
    Year = as.integer(Year),
    Mean_SBP = as.numeric(str_extract(SBP_values, "^\\d+\\.\\d+")),
    Lower_CI = as.numeric(str_match(SBP_values, "\\((\\d+\\.\\d+)-")[, 2]),
    Upper_CI = as.numeric(str_match(SBP_values, "-(\\d+\\.\\d+)\\)")[, 2])
  ) %>%
  select(Region, Country, Year, Sex, Mean_SBP, Lower_CI, Upper_CI)


# Import and reshape GDP per capita data

gdp_raw <- read_csv("API_NY.GDP.PCAP.KD_DS2_en_csv_v2_89 (1).csv", 
                                                   skip = 4,
                                                   show_col_types = FALSE )

gdp_data <- gdp_raw %>%
  select(
    Country = `Country Name`,
    ISO = `Country Code`,
    starts_with("19"),
    starts_with("20")
  ) %>%
  pivot_longer(
    cols = starts_with(c("19", "20")),
    names_to = "Year",
    values_to = "GDP_per_capita"
  ) %>%
  mutate(
    Year = as.integer(Year),
    GDP_per_capita = as.numeric(GDP_per_capita)
  ) %>%
  filter(!is.na(GDP_per_capita))

# Harmonize country names across SBP and GDP data sets
setdiff(sort(unique(sbp_analysis$Country)),
  sort(unique(gdp_data$Country)))

sort(unique(gdp_data$Country))

gdp_data <- gdp_data %>%
  mutate(Country = recode(Country,
                          "Bahamas, The" = "Bahamas",
                          "Virgin Islands (U.S.)" = "British Virgin Islands",
                          "Cote d'Ivoire" = "Cote d'Ivoire",
                          "Cabo Verde" = "Cape Verde",
                          "Hong Kong SAR, China" = "China (Hong Kong SAR)",
                          "Macao SAR, China" = "China (Macao SAR)",
                          "Congo, Rep." = "Congo",
                          "Czechia" = "Czech Republic",
                          "Congo, Dem. Rep." = "Democratic Republic of the Congo",
                          "Egypt, Arab Rep." = "Egypt",
                          "Gambia, The" = "Gambia",
                          "Iran, Islamic Rep." = "Iran ",
                          "Kyrgyz Republic" = "Kyrgyzstan",
                          "Lao PDR" = "Lao People's Democratic Republic",
                          "Libya" = "Libya",
                          "North Macedonia" = "North Macedonia",
                          "Micronesia, Fed. Sts." = "Micronesia",
                          "Sint Maarten (Dutch part)" = "Netherlands Antilles",
                          "West Bank and Gaza" = "Palestine",
                          "Puerto Rico (US)" = "Puerto Rico",
                          "Korea, Rep." = "Republic of Korea",
                          "Sao Tome and Principe" = "Sao Tomè and Principe",
                          "St. Kitts and Nevis" = "Saint Kitts and Nevis",
                          "St. Lucia" = "Saint Lucia",
                          "St. Vincent and the Grenadines" = "Saint Vincent and the Grenadines",
                          "Slovak Republic" = "Slovakia",
                          "Somalia, Fed. Rep." = "Somalia",
                          "Eswatini" = "Swaziland",
                          "Turkiye" = "Turkey",
                          "Tanzania" = "Tanzania",
                          "United States" = "United States of America",
                          "Venezuela, RB" = "Venezuela",
                          "Yemen, Rep." = "Yemen"
                          ))

# ---------------------------
# Merge SBP and GDP data
# ---------------------------

# Merge data and compute log-transformed GDP
sbp_gdp <- sbp_analysis %>%
  filter(!Country %in% c(
    "Cook Islands",
    "Democratic People's Republic of Korea",
    "Taiwan"
  )) %>%
  left_join(
    gdp_data,
    by = c("Country", "Year")
  ) %>%
  mutate(
    log_GDP = log(GDP_per_capita)
  )


# Restrict analyses to countries with available baseline data in 1980
countries_1980 <- sbp_gdp %>%
  filter(Year == 1980, !is.na(Mean_SBP), !is.na(log_GDP)) %>%
  pull(Country) %>%
  unique()

# Prepare variables for analysis
sbp_gdp_1980sample <- sbp_gdp %>%
  filter(Country %in% countries_1980) %>%
  mutate(
    Sex = factor(Sex, levels = c("Female", "Male")),
    Year_c = Year - 1980
  )

#Check
Desc(sbp_gdp_1980sample$Country, maxrows = Inf, plotit = F)

# Effect coding: Female = -0.5, Male = 0.5
contrasts(sbp_gdp_1980sample$Sex) <- c(-0.5, 0.5)

# Check sample composition by year
sbp_gdp_1980sample %>%
  group_by(Year) %>%
  summarise(
    n_countries = n_distinct(Country),
    n_rows = n(),
    missing_gdp = sum(is.na(GDP_per_capita)),
    .groups = "drop"
  )


# Create UN_regions directly in sbp_gdp_1980sample using UN major regions.

sbp_gdp_1980sample <- sbp_gdp_1980sample %>%
  mutate(
    UN_regions = countrycode(
      Country,
      origin = "country.name",
      destination = "un.regionsub.name"
    ),
    
    UN_intermediate = countrycode(
      Country,
      origin = "country.name",
      destination = "un.regionintermediate.name"
    ),
    
    UN_regions_pure22 = coalesce(
      UN_intermediate,
      UN_regions
    )
  )


table(sbp_gdp_1980sample$UN_regions_pure22, useNA = "ifany")#21 regions
table(sbp_gdp_1980sample$UN_regions, useNA= "ifany")


# ---------------------------
# Main mixed-effects analyses
# ---------------------------

#log_gdp should be converted in z_scores (scaled: mean= 0, sd=1) to make the intercept more interpretable.
#If log_gdp is not scaled, the intercept would correspond to: female, 1980 (or zero), when log_gdp= 0, but log_gdp ranges from 5 to 11.
range(sbp_gdp_1980sample$log_GDP, na.rm = TRUE)#from 170 to 129.525 USD
mean(sbp_gdp_1980sample$log_GDP)#4315 USD

#If we scaled log_gdp, the intercept would correspond to: female, 1980(or zero), when log_g= 0 (which is the average log_gdp = exp(8.37)= GDP_PER_CAPITA= 4315 USD) 

#Creation of log_GDP_z
sbp_gdp_1980sample$log_GDP_z<- scale(sbp_gdp_1980sample$log_GDP, center = T, scale = T)
range(sbp_gdp_1980sample$log_GDP_z, na.rm = TRUE)


#Step 0: Null Model

M0<- lmer(Mean_SBP ~ 1 + (1 | Country),
          data = sbp_gdp_1980sample, 
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000)),
          REML= FALSE)

summary(M0)
insight::get_variance(M0)
performance::icc(M0)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row0 <- as.numeric(rownames(model.frame(M0)))

vcov_region0 <- vcovCR(M0, 
                       cluster = sbp_gdp_1980sample$UN_regions_pure22[row0], 
                       type = "CR2")
coef_test(M0, vcov = vcov_region0)


#Step 1: Fixed Main Effect effects + Random Intercepts and Random Slopes
M1<- lmer(Mean_SBP ~ Sex + log_GDP_z + Year_c + (1 + Sex + Year_c | Country),
          data = sbp_gdp_1980sample, 
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000)),
          REML= FALSE)

summary(M1)
insight::get_variance(M1)
performance::icc(M1)




#Controlling for potential geographical and cultural non-independence= robust standard errors
rowa <- as.numeric(rownames(model.frame(M1)))

vcov_regiona <- vcovCR(M1, 
                      cluster = sbp_gdp_1980sample$UN_regions_pure22[rowa], 
                      type = "CR2")


coef_test(M1, vcov = vcov_regiona)


#Confidence Interval 95%
IC <- conf_int(
  M1,
  vcov = vcov_regiona,
  coefs = c("log_GDP_z"), 
  level = 0.95)
print(IC)


#Step 3: Full Interactions + Random Slopes and Random Intercepts
M2<- lmer(Mean_SBP ~ Sex * log_GDP_z * Year_c + (1 + Sex + Year_c | Country),
           data = sbp_gdp_1980sample,
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000)),
           REML= FALSE)
summary(M2)
insight::get_variance(M2)
performance::icc(M2)

#Optimizer control
model_list1<- allFit(M2)
summary(model_list1)




#Model comparison
anova(M0, M1, M2)


# #Controlling for potential geographical and cultural non-independence= robust standard errors
row <- as.numeric(rownames(model.frame(M2)))
vcov_region <- vcovCR(M2,
                       cluster = sbp_gdp_1980sample$UN_regions_pure22[row],
                      type = "CR2")
# 
# 
coef_test(M2, vcov = vcov_region)
car::vif(M2)

# ---------------------------
# Model-based predictions
# ---------------------------

# Plot predicted SBP by GDP, sex, and year

pred_data_1980 <- expand.grid(
  log_GDP_z = seq(
    min(sbp_gdp_1980sample$log_GDP_z, na.rm = TRUE),
    max(sbp_gdp_1980sample$log_GDP_z, na.rm = TRUE),
    length.out = 100
  ),
  Year_c = c(0, 10, 20, 28),
  Sex = c("Female", "Male")
)

pred_data_1980 <- pred_data_1980 %>%
  mutate(
    Year = Year_c + 1980,
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

# Apply effect coding: Female = -0.5, Male = +0.5

contrasts(pred_data_1980$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

# Design matrix with effect coding

X <- model.matrix(
  ~ Sex * log_GDP_z * Year_c,
  data = pred_data_1980
)

# Extract fixed-effect coefficients

beta <- fixef(M2)

# Ensure matching coefficient order

X <- X[, names(beta), drop = FALSE]

# Predicted values

pred_data_1980$fit <- as.vector(X %*% beta)

# Cluster-robust variance-covariance matrix

V <- as.matrix(vcov_region)
V <- V[names(beta), names(beta), drop = FALSE]

pred_data_1980$se <- sqrt(diag(X %*% V %*% t(X)))

# 95% confidence intervals (normal approximation)

pred_data_1980 <- pred_data_1980 %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

# Plot

ggplot(pred_data_1980,
       aes(x = log_GDP_z,
           y = fit,
           color = Sex,
           fill = Sex)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.20,
              color = NA) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~Year) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Log GDP per capita",
    y = "Predicted SBP (mmHg)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )
# ---------------------------
# Robustness analyses
# ---------------------------

# AR(1) model for temporal autocorrelation
sbp_gdp_1980sample <- sbp_gdp_1980sample %>%
  mutate(
    Country = factor(Country)
  )

model_ar1_1980sample <- lme(
  Mean_SBP ~ Sex * log_GDP_z * Year_c,
  random = ~ 1 + Sex + Year_c | Country,
  correlation = corAR1(form = ~ Year_c | Country/Sex),
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  data = sbp_gdp_1980sample,
  method = "ML"
)
summary(model_ar1_1980sample)

# Sex-specific AR(1) robustness checks
sbp_male <- sbp_gdp_1980sample %>%
  filter(Sex == "Male") %>%
  droplevels()

sbp_female <- sbp_gdp_1980sample %>%
  filter(Sex == "Female") %>%
  droplevels()

# Male model with AR(1)
model_ar1_male <- lme(
  Mean_SBP ~ log_GDP_z * Year_c,
  random = ~ 1  + Year_c  | Country,
  correlation = corAR1(form = ~ Year_c | Country),
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  data = sbp_male,
  method = "ML",
  na.action = na.omit
)

summary(model_ar1_male)

# Female model with AR(1)
model_ar1_female <- lme(
  Mean_SBP ~ log_GDP_z * Year_c,
  random = ~ 1 + Year_c | Country,
  correlation = corAR1(form = ~ Year_c | Country),
  data = sbp_female,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  method = "ML",
  na.action = na.omit
)

summary(model_ar1_female)


#--------------------------
# Spatial Autocorrelation of Residuals: Moran's I
#--------------------------
library(dplyr)
library(spdep)
library(lme4)

# 1. M2 residuals
dati_modello <- model.frame(M2)

dati_modello$residuo <- residuals(M2)



lookup_regionis <- sbp_gdp_1980sample %>%
  select(Country, UN_regions_pure22) %>%
  mutate(Country = as.character(Country)) %>%
  distinct()


stopifnot(!anyDuplicated(lookup_regionis$Country))

dati_modello <- dati_modello %>%
  mutate(Country = as.character(Country)) %>%
  left_join(lookup_regionis, by = "Country")

stopifnot(!anyNA(dati_modello$UN_regions_pure22))

# 3. average residual for each country 
residui_country <- dati_modello %>%
  group_by(Country, UN_regions_pure22) %>%
  summarise(
    residuo = mean(residuo, na.rm = TRUE),
    .groups = "drop"
  )

# 4. UN macro-regions cluster matrix
cluster <- as.character(
  residui_country$UN_regions_pure22
)


W <- outer(cluster, cluster, "==") * 1


diag(W) <- 0


listw <- spdep::mat2listw(
  W,
  style = "W",
  zero.policy = TRUE
)

# 6. Global Moran's I
moran_result <- spdep::moran.test(
  residui_country$residuo,
  listw,
  alternative = "two.sided",
  zero.policy = TRUE
)

print(moran_result)

#--------------------------
# Longitudinal MLM 3 way interaction + variance decomposition LogGDP (Mundlak decomposition)
#--------------------------
#Let us create a new data-frame where Log_gdp_z is decomposed within and between country


# One observation per country-year
gdp_country_year <- sbp_gdp_1980sample %>%
  distinct(Country, Year_c, log_GDP_z)

# Ensure GDP is unique for each country-year
stopifnot(
  nrow(gdp_country_year) ==
    nrow(distinct(gdp_country_year, Country, Year_c))
)

# Decompose into within and between components
gdp_country_year <- gdp_country_year %>%
  group_by(Country) %>%
  mutate(
    gdp_between = mean(log_GDP_z, na.rm = TRUE),
    gdp_within  = log_GDP_z - gdp_between
  ) %>%
  ungroup()

# Merge into the original dataset
sbp_gdp_1980sample <- sbp_gdp_1980sample %>%
  select(-any_of(c("gdp_within", "gdp_between",
                   "gdp_within_std", "gdp_between_std"))) %>%
  left_join(
    gdp_country_year %>%
      select(Country, Year_c, gdp_within, gdp_between),
    by = c("Country", "Year_c")
  )
#--------------------------------
#Longitudinal Multilevel Model 3 way interaction with log_gdp_c variance decomposition + Random Intercept and Random Slopes
#--------------------------------
#Null Model
mod_mundlak_0 <- lmer(
  Mean_SBP ~ 1 +
    (1 |  Country),
  data = sbp_gdp_1980sample,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)

summary(mod_mundlak_0)
performance::icc(mod_mundlak_0)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row2ax <- as.numeric(rownames(model.frame(mod_mundlak_0)))

vcov_region2ax <- vcovCR(mod_mundlak_0, 
                        cluster = sbp_gdp_1980sample$UN_regions_pure22[row2ax], 
                        type = "CR2")

coef_test(mod_mundlak_0, vcov = vcov_region2ax)

#Main Fixed Effect Model
mod_mundlak_1 <- lmer(
  Mean_SBP ~ Sex + Year_c + (gdp_within + gdp_between) +
    (1  + Sex + Year_c | Country),
  data = sbp_gdp_1980sample,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)



summary(mod_mundlak_1)
performance::icc(mod_mundlak_1)


#Controlling for potential geographical and cultural non-independence= robust standard errors
row2a <- as.numeric(rownames(model.frame(mod_mundlak_1)))

vcov_region2a <- vcovCR(mod_mundlak_1, 
                       cluster = sbp_gdp_1980sample$UN_regions_pure22[row2a], 
                       type = "CR2")

coef_test(mod_mundlak_1, vcov = vcov_region2a)

#Optimizer Check
lismod1<-allFit(mod_mundlak_1)
summary(lismod1)

#Full Interaction Effect Model
mod_mundlak_2 <- lmer(
  Mean_SBP ~ Sex * Year_c * (gdp_within + gdp_between) +
    (1 + Sex + Year_c  | Country),
  data = sbp_gdp_1980sample,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)

summary(mod_mundlak_2)
performance::icc(mod_mundlak_2)
insight::get_variance(mod_mundlak_2)

#Optimizer Check
lismod2<-allFit(mod_mundlak_2)
summary(lismod2)
#---------------------------------
# Mundlak Models Comparison
#---------------------------------
anova(mod_mundlak_0, mod_mundlak_1, mod_mundlak_2)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row2 <- as.numeric(rownames(model.frame(mod_mundlak_2)))

vcov_region2 <- vcovCR(mod_mundlak_2, 
                      cluster = sbp_gdp_1980sample$UN_regions_pure22[row2], 
                      type = "CR2")

coef_test(mod_mundlak_2, vcov = vcov_region2)

car::vif(mod_mundlak_2)

cor(sbp_gdp_1980sample$gdp_within, sbp_gdp_1980sample$Year_c)
#------------------------------------------------
# Three way interaction: WITHIN-COUNTRY GDP PLOT
# Sex × Year × Within-country GDP interaction
# Predictions restricted to the observed GDP range
# within each year
#------------------------------------------------



#---------------------------------------------------------
# 1. Observed range of within-country GDP in each year
#---------------------------------------------------------

year_ranges <- sbp_gdp_1980sample %>%
  distinct(Country, Year_c, gdp_within) %>%
  filter(
    is.finite(Year_c),
    is.finite(gdp_within)
  ) %>%
  group_by(Year_c) %>%
  summarise(
    min_gdp = min(gdp_within, na.rm = TRUE),
    max_gdp = max(gdp_within, na.rm = TRUE),
    .groups = "drop"
  )

print(year_ranges)


#---------------------------------------------------------
# 2. Create prediction dataset
#---------------------------------------------------------

new_within <- year_ranges %>%
  mutate(
    gdp_within = map2(
      min_gdp,
      max_gdp,
      ~ seq(.x, .y, length.out = 100)
    )
  ) %>%
  select(Year_c, gdp_within) %>%
  unnest(gdp_within) %>%
  crossing(
    Sex = c("Female", "Male")
  ) %>%
  mutate(
    
    # Hold between-country GDP at zero
    gdp_between = 0,
    
    # Calendar year
    Year = 1980 + Year_c,
    
    # Correct effect coding
    # Female = -0.5; Male = +0.5
    Sex1 = ifelse(
      Sex == "Male",
      0.5,
      -0.5
    ),
    
    # Labels for plot
    Sex_label = factor(
      Sex1,
      levels = c(-0.5, 0.5),
      labels = c("Female", "Male")
    )
  )


#---------------------------------------------------------
# 3. Extract fixed-effect coefficients
#---------------------------------------------------------

betaw <- fixef(mod_mundlak_2)

print(betaw)


#---------------------------------------------------------
# 4. Build design matrix manually
#    Using correct effect coding
#---------------------------------------------------------

Xw <- cbind(
  
  "(Intercept)" = 1,
  
  "Sex1" =
    new_within$Sex1,
  
  "Year_c" =
    new_within$Year_c,
  
  "gdp_within" =
    new_within$gdp_within,
  
  "gdp_between" =
    new_within$gdp_between,
  
  "Sex1:Year_c" =
    new_within$Sex1 *
    new_within$Year_c,
  
  "Sex1:gdp_within" =
    new_within$Sex1 *
    new_within$gdp_within,
  
  "Sex1:gdp_between" =
    new_within$Sex1 *
    new_within$gdp_between,
  
  "Year_c:gdp_within" =
    new_within$Year_c *
    new_within$gdp_within,
  
  "Year_c:gdp_between" =
    new_within$Year_c *
    new_within$gdp_between,
  
  "Sex1:Year_c:gdp_within" =
    new_within$Sex1 *
    new_within$Year_c *
    new_within$gdp_within,
  
  "Sex1:Year_c:gdp_between" =
    new_within$Sex1 *
    new_within$Year_c *
    new_within$gdp_between
)


#---------------------------------------------------------
# 5. Match matrix order to model coefficients
#---------------------------------------------------------

Xw <- Xw[
  ,
  names(betaw),
  drop = FALSE
]


#---------------------------------------------------------
# 6. Variance-covariance matrix
#---------------------------------------------------------

Vw <- as.matrix(vcov_region2)

Vw <- Vw[
  names(betaw),
  names(betaw),
  drop = FALSE
]


#---------------------------------------------------------
# 7. Fixed-effect predictions
#---------------------------------------------------------

new_within$predicted_SBP <- as.numeric(
  Xw %*% betaw
)


#---------------------------------------------------------
# 8. Standard errors of predictions
#---------------------------------------------------------

new_within$SE <- sqrt(
  rowSums(
    (Xw %*% Vw) * Xw
  )
)


#---------------------------------------------------------
# 9. 95% confidence intervals
#---------------------------------------------------------

new_within <- new_within %>%
  mutate(
    lower_CI = predicted_SBP - 1.96 * SE,
    upper_CI = predicted_SBP + 1.96 * SE
  )


#---------------------------------------------------------
# 10. Prepare observed data for plotting
#---------------------------------------------------------

observed_within <- sbp_gdp_1980sample %>%
  filter(is.finite(Year_c)) %>%
  mutate(
    Year = 1980 + Year_c,
    
    Sex_label = factor(
      Sex,
      levels = c("Female", "Male")
    )
  )


#---------------------------------------------------------
# 11. Global x-axis range
#---------------------------------------------------------

x_limits <- range(
  sbp_gdp_1980sample$gdp_within,
  na.rm = TRUE
)


#---------------------------------------------------------
# 12. Plot
#---------------------------------------------------------

ggplot() +
  
  geom_point(
    data = observed_within,
    aes(
      x = gdp_within,
      y = Mean_SBP,
      color = Sex_label
    ),
    alpha = 0.12,
    size = 1
  ) +
  
  geom_ribbon(
    data = new_within,
    aes(
      x = gdp_within,
      ymin = lower_CI,
      ymax = upper_CI,
      fill = Sex_label,
      group = Sex_label
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  geom_line(
    data = new_within,
    aes(
      x = gdp_within,
      y = predicted_SBP,
      color = Sex_label,
      group = Sex_label
    ),
    linewidth = 1.25
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  facet_wrap(
    ~ Year,
    nrow = 1
  ) +
  
  scale_x_continuous(
    limits = x_limits,
    breaks = seq(
      floor(x_limits[1]),
      ceiling(x_limits[2]),
      by = 1
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  labs(
    x = "Within-country standardized Log GDP",
    y = "Predicted SBP (mmHg)",
    color = "Sex",
    fill = "Sex"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "plain"
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    )
  )



#----------------------------------------
# Three-way interaction effect plot:
# Between-country GDP
#----------------------------------------

#---------------------------------------------------------
# 1. Create prediction dataset
#---------------------------------------------------------

new_between <- expand_grid(
  
  gdp_between = seq(
    min(sbp_gdp_1980sample$gdp_between, na.rm = TRUE),
    max(sbp_gdp_1980sample$gdp_between, na.rm = TRUE),
    length.out = 100
  ),
  
  Year_c = c(0, 10, 20, 28),
  
  Sex = c("Female", "Male"),
  
  # Hold within-country GDP constant at zero
  gdp_within = 0
)


#---------------------------------------------------------
# 2. Create labels
#---------------------------------------------------------

new_between <- new_between %>%
  mutate(
    
    Year = 1980 + Year_c,
    
    Sex = factor(
      Sex,
      levels = c("Female", "Male")
    ),
    
    Sex_label = factor(
      Sex,
      levels = c("Female", "Male")
    )
  )


#---------------------------------------------------------
# 3. Apply correct effect coding
#    Female = -0.5; Male = +0.5
#---------------------------------------------------------

contrasts(new_between$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(
    c("Female", "Male"),
    "1"
  )
)


#---------------------------------------------------------
# 4. Design matrix for fixed effects
#---------------------------------------------------------

Xb <- model.matrix(
  ~ Sex * Year_c * (gdp_within + gdp_between),
  data = new_between
)


#---------------------------------------------------------
# 5. Fixed-effect coefficients
#---------------------------------------------------------

betab <- fixef(mod_mundlak_2)

# Match design matrix to model coefficients

Xb <- Xb[
  ,
  names(betab),
  drop = FALSE
]


#---------------------------------------------------------
# 6. Variance-covariance matrix
#---------------------------------------------------------

Vb <- as.matrix(vcov_region2)

Vb <- Vb[
  names(betab),
  names(betab),
  drop = FALSE
]


#---------------------------------------------------------
# 7. Predicted SBP
#---------------------------------------------------------

new_between$predicted_SBP <- as.numeric(
  Xb %*% betab
)


#---------------------------------------------------------
# 8. Standard errors
#---------------------------------------------------------

new_between$SE <- sqrt(
  rowSums(
    (Xb %*% Vb) * Xb
  )
)


#---------------------------------------------------------
# 9. 95% confidence intervals
#---------------------------------------------------------

new_between <- new_between %>%
  mutate(
    lower_CI = predicted_SBP - 1.96 * SE,
    upper_CI = predicted_SBP + 1.96 * SE
  )


#---------------------------------------------------------
# 10. Plot
#---------------------------------------------------------

ggplot(
  new_between,
  aes(
    x = gdp_between,
    y = predicted_SBP,
    color = Sex_label,
    fill = Sex_label
  )
) +
  
  geom_ribbon(
    aes(
      ymin = lower_CI,
      ymax = upper_CI
    ),
    alpha = 0.20,
    color = NA
  ) +
  
  geom_line(linewidth = 1.3) +
  
  facet_wrap(~ Year) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = "Between-country standardized Log GDP",
    y = "Predicted SBP (mmHg)"
  ) +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )




#---------------------------------
# Fixed-effect regression (fixed effects: country + time) 

#---------------------------------

mod_fe_base <- fixest::feols(
  Mean_SBP ~ Sex*log_GDP_z |
    Country + Year_c,
  data = sbp_gdp_1980sample,
  cluster =  ~ UN_regions_pure22, 
)

summary(mod_fe_base)


#Plot



#--------------------------------------------------
# 1. Prediction grid
#--------------------------------------------------

plot_fe <- expand_grid(
  
  log_GDP_z = seq(
    min(sbp_gdp_1980sample$log_GDP_z, na.rm = TRUE),
    max(sbp_gdp_1980sample$log_GDP_z, na.rm = TRUE),
    length.out = 100
  ),
  
  Sex = levels(sbp_gdp_1980sample$Sex)
)

#--------------------------------------------------
# 2. Create sex indicator
#    Female = -0.5, Male = +0.5
#--------------------------------------------------

plot_fe <- plot_fe %>%
  mutate(
    Sex1 = ifelse(Sex == "Male", 0.5, -0.5),
    
    Sex_label = factor(
      Sex1,
      levels = c(-0.5, 0.5),
      labels = c("Female", "Male")
    )
  )

#--------------------------------------------------
# 3. Extract coefficients
#--------------------------------------------------

betafe <- coef(mod_fe_base)

betafe

#--------------------------------------------------
# 4. Design matrix
#--------------------------------------------------

Xfe <- cbind(
  
  Sex1 = plot_fe$Sex1,
  
  log_GDP_z = plot_fe$log_GDP_z,
  
  "Sex1:log_GDP_z" =
    plot_fe$Sex1 * plot_fe$log_GDP_z
)

# Check
colnames(Xfe)
names(betafe)

# Ensure same order
Xfe <- Xfe[
  ,
  names(betafe),
  drop = FALSE
]

#--------------------------------------------------
# 5. Variance-covariance matrix
#--------------------------------------------------

Vfe <- vcov(mod_fe_base)

# Ensure same order
Vfe <- Vfe[
  names(betafe),
  names(betafe),
  drop = FALSE
]

#--------------------------------------------------
# 6. Adjusted fitted component
#--------------------------------------------------

plot_fe$fit <- as.numeric(
  Xfe %*% betafe
)

#--------------------------------------------------
# 7. Standard errors
#--------------------------------------------------

plot_fe$se <- sqrt(
  rowSums((Xfe %*% Vfe) * Xfe)
)

#--------------------------------------------------
# 8. 95% confidence intervals
#--------------------------------------------------

plot_fe <- plot_fe %>%
  mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )

#--------------------------------------------------
# 9. Plot
#--------------------------------------------------

ggplot(
  plot_fe,
  aes(
    x = log_GDP_z,
    y = fit,
    color = Sex_label,
    fill = Sex_label
  )
) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.15,
    color = NA
  ) +
  
  geom_line(
    linewidth = 1.3
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  labs(
    x = "Standardized Log GDP",
    y = "Predicted difference in mean SBP (mmHg)",
    color = "Sex",
    fill = "Sex"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15
    ),
    axis.title = element_text(
      face = "plain"
    ),
    legend.position = "bottom",
    legend.title = element_text(
      face = "bold"
    )
  )








#-----------------------------------------------------
#-----------------------------------------------------
# HDI 
#-----------------------------------------------------
#-----------------------------------------------------
library(readr)
HDI <- read_csv("human-development-index (1).csv")

HDI_2<- HDI %>% rename(Country=Entity)

library(DescTools)

Desc(HDI_2$Country, plotit = F, maxrows = Inf)

HDI_3<- HDI_2 %>% filter (!Country %in% c("World", "North America","Arab States (UNDP)","Africa", "South America", "Asia", "European Union (27)" ,"East Asia and the Pacific (UNDP)", "Sub-Saharan Africa (UNDP)" , "Europe and Central Asia (UNDP)", " European Union (27)", "High-income countries", "High human development (UNDP)", "Latin America and the Caribbean (UNDP)", "Low human development (UNDP)", "Medium human development (UNDP)", "South Asia (UNDP)", " Sub-Saharan Africa (UNDP)", "Very high human development (UNDP)", " Lower-middle-income countries", "Micronesia (country)", " North America", "Upper-middle-income countries", "Lower-middle-income countries", "Low-income countries"))


Desc(HDI_3$Country, plotit = F, maxrows = Inf)


HDI_4<- HDI_3[, -c(2,5)]


HDI_5 <- subset(HDI_4, HDI_4$Year %in% c(1990, 2000, 2008))


Desc(sbp_gdp_1980sample$Country, maxrows = Inf, plotit = F)
Desc(HDI_5$Country, maxrows = Inf, plotit = F)

library(dplyr)

HDI_5$Country <- recode(
  HDI_5$Country,
  "Brunei"      = "Brunei Darussalam",
  "Hong Kong"   = "China (Hong Kong SAR)",
  "South Korea" = "Republic of Korea",
  "Eswatini"    = "Swaziland"
)

sbp_gdp_1980sampleHDI<- full_join(sbp_gdp_1980sample, HDI_5, by= c("Country", "Year"))

sbp_gdp_1980sample_2 <- sbp_gdp_1980sampleHDI %>%
  filter(Year %in% c(1990, 2000, 2008))

View(sbp_gdp_1980sample_2)
sbp_gdp_1980sample_2<- sbp_gdp_1980sample_2[, c(2, 3, 4, 5,6, 7,9, 10,11, 12,14, 16, 17, 18)]


sbp_gdp_1980sample_2<- sbp_gdp_1980sample_2 %>% rename(HDI =`Human Development Index`)
sbp_gdp_1980sample_2$Year_c<- (sbp_gdp_1980sample_2$Year - 1990)


sbp_gdp_1980sample_2$HDI_z<- scale(sbp_gdp_1980sample_2$HDI)

View(sbp_gdp_1980sample_2)

# Effect coding: Female = -0.5; Male = +0.5
sbp_gdp_1980sample_2$Sex <- factor(
  sbp_gdp_1980sample_2$Sex,
  levels = c("Female", "Male")
)

contrasts(sbp_gdp_1980sample_2$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

#----------------------------
# Main mixed effect model
#----------------------------
M0_HDI <- lmer(Mean_SBP~1 + (1|Country), 
               data= sbp_gdp_1980sample_2)
summary(M0_HDI)
performance::icc(M0_HDI)


#Controlling for potential geographical and cultural non-independence= robust standard errors
row3as <- as.numeric(rownames(model.frame(M0_HDI)))

vcov_region3s <- vcovCR(M0_HDI, 
                        cluster = sbp_gdp_1980sample_2$UN_regions_pure22[row3as], 
                        type = "CR2")

coef_test(M0_HDI, vcov = vcov_region3s)

#Main effect model with random intercepts

M1_HDI <- lmer(Mean_SBP ~ Sex + HDI_z + Year_c + (1 |Country),
               data= sbp_gdp_1980sample_2,
               control = lmerControl(
                 optimizer = "bobyqa",
                 optCtrl = list(maxfun = 200000)), 
               REML = FALSE)

summary(M1_HDI)

performance::icc(M1_HDI)

#optimizer check
lism1hdi<- allFit(M1_HDI)
summary(lism1hdi)



#Controlling for potential geographical and cultural non-independence= robust standard errors
row3a <- as.numeric(rownames(model.frame(M1_HDI)))

vcov_region3 <- vcovCR(M1_HDI, 
                       cluster = sbp_gdp_1980sample_2$UN_regions_pure22[row3a], 
                       type = "CR2")

coef_test(M1_HDI, vcov = vcov_region3)

#Full Interaction Effect Model with Random Intercepts

M2_HDI <- lmer(
  Mean_SBP ~ Sex * HDI_z * Year_c + (1 | Country),
  data = sbp_gdp_1980sample_2, 
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)

summary(M2_HDI)
performance::icc(M2_HDI)


#optimizer check
lism2hdi<- allFit(M2_HDI)
summary(lism2hdi)

#Model comparison

anova(M1_HDI, M2_HDI)



# #Controlling for potential geographical and cultural non-independence= robust standard errors
row3 <- as.numeric(rownames(model.frame(M2_HDI)))
# 
vcov_region3 <- vcovCR(M2_HDI, 
                       cluster = sbp_gdp_1980sample_2$UN_regions_pure22[row3], 
                       type = "CR2")#CR2 because of the small number of clusters

coef_test(M2_HDI, vcov = vcov_region3)




#-------------------------------
# plot
#-------------------------------


# Plot predicted SBP by HDI_z, sex, and year
pred_data_1980_HDI <- expand.grid(
  HDI_z = seq(
    min(sbp_gdp_1980sample_2$HDI_z, na.rm = TRUE),
    max(sbp_gdp_1980sample_2$HDI_z, na.rm = TRUE),
    length.out = 100
  ),
  Year_c = c(0, 10, 18),
  Sex = c("Female", "Male")
)

pred_data_1980_HDI <- pred_data_1980_HDI %>%
  mutate(
    Year = Year_c + 1990,
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

# Same effect coding as in the fitted models
contrasts(pred_data_1980_HDI$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

X_HDI <- model.matrix(
  ~ Sex * HDI_z * Year_c,
  data = pred_data_1980_HDI
)

beta_HDI <- fixef(M2_HDI)
X_HDI <- X_HDI[, names(beta_HDI), drop = FALSE]
pred_data_1980_HDI$fit <- X_HDI %*% beta_HDI
V_HDI <- as.matrix(vcov_region3)
V_HDI <- V_HDI[names(beta_HDI), names(beta_HDI), drop = FALSE]
pred_data_1980_HDI$se <- sqrt(diag(X_HDI %*% V_HDI %*% t(X_HDI)))

pred_data_1980_HDI <- pred_data_1980_HDI %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

ggplot(pred_data_1980_HDI,
       aes(x = HDI_z,
           y = fit,
           color = Sex,
           fill = Sex)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.20,
              color = NA) +
  
  geom_line(linewidth = 1.3) +
  facet_wrap(~Year) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Human Development Index (z_scores)",
    y = "Predicted SBP (mmHg)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )


#=============================================
# gdp 5 and 10 years lag: SUPPLEMENTARY ANALYSES
#=============================================
sbp_analysis_10lag <- raw_sbp %>%
  select(1:9) %>%
  rename(
    Country = X1,
    male_1980 = X2, male_1990 = X3, male_2000 = X4, male_2008 = X5,
    female_1980 = X6, female_1990 = X7, female_2000 = X8, female_2008 = X9
  ) %>%
  mutate(
    Region = if_else(
      is.na(Country) & !is.na(male_1980),
      male_1980,
      NA_character_
    )
  ) %>%
  fill(Region, .direction = "down") %>%
  filter(
    !is.na(Country),
    Country != "Country",
    !is.na(male_1980)
  ) %>%
  pivot_longer(
    cols = matches("^(male|female)_\\d{4}$"),
    names_to = c("Sex", "Year"),
    names_pattern = "(male|female)_(\\d{4})",
    values_to = "SBP_values"
  ) %>%
  mutate(
    Sex = recode(Sex, male = "Male", female = "Female"),
    Year = as.integer(Year),
    Mean_SBP = as.numeric(str_extract(SBP_values, "^\\d+\\.\\d+")),
    Lower_CI = as.numeric(str_match(SBP_values, "\\((\\d+\\.\\d+)-")[, 2]),
    Upper_CI = as.numeric(str_match(SBP_values, "-(\\d+\\.\\d+)\\)")[, 2])
  ) %>%
  select(Region, Country, Year, Sex, Mean_SBP, Lower_CI, Upper_CI)

# 
gdp_base <- gdp_raw %>%
  select(
    Country = `Country Name`,
    ISO = `Country Code`,
    starts_with("19"),
    starts_with("20")
  ) %>%
  pivot_longer(
    cols = starts_with(c("19", "20")),
    names_to = "Year",
    values_to = "GDP_per_capita"
  ) %>%
  mutate(
    GDP_per_capita = as.numeric(GDP_per_capita),
    Year = as.numeric(Year)
  ) %>%
  filter(!is.na(GDP_per_capita))


gdp_base <- gdp_base %>%
  mutate(Country = recode(Country,
                          "Bahamas, The" = "Bahamas",
                          "Virgin Islands (U.S.)" = "British Virgin Islands",
                          "Cote d'Ivoire" = "Cote d'Ivoire",
                          "Cabo Verde" = "Cape Verde",
                          "Hong Kong SAR, China" = "China (Hong Kong SAR)",
                          "Macao SAR, China" = "China (Macao SAR)",
                          "Congo, Rep." = "Congo",
                          "Czechia" = "Czech Republic",
                          "Congo, Dem. Rep." = "Democratic Republic of the Congo",
                          "Egypt, Arab Rep." = "Egypt",
                          "Gambia, The" = "Gambia",
                          "Iran, Islamic Rep." = "Iran ",
                          "Kyrgyz Republic" = "Kyrgyzstan",
                          "Lao PDR" = "Lao People's Democratic Republic",
                          "Libya" = "Libya",
                          "North Macedonia" = "North Macedonia",
                          "Micronesia, Fed. Sts." = "Micronesia",
                          "Sint Maarten (Dutch part)" = "Netherlands Antilles",
                          "West Bank and Gaza" = "Palestine",
                          "Puerto Rico (US)" = "Puerto Rico",
                          "Korea, Rep." = "Republic of Korea",
                          "Sao Tome and Principe" = "Sao Tomè and Principe",
                          "St. Kitts and Nevis" = "Saint Kitts and Nevis",
                          "St. Lucia" = "Saint Lucia",
                          "St. Vincent and the Grenadines" = "Saint Vincent and the Grenadines",
                          "Slovak Republic" = "Slovakia",
                          "Somalia, Fed. Rep." = "Somalia",
                          "Eswatini" = "Swaziland",
                          "Turkiye" = "Turkey",
                          "Tanzania" = "Tanzania",
                          "United States" = "United States of America",
                          "Venezuela, RB" = "Venezuela",
                          "Yemen, Rep." = "Yemen"
  ))



gdp_10lag <- gdp_base %>%
  mutate(Year = Year + 10) %>%
  rename(GDP_10lag = GDP_per_capita)


gdp_5lag <- gdp_base %>%
  mutate(Year = Year + 5) %>%
  rename(GDP_5lag = GDP_per_capita)


sbp_gdp_combined <- sbp_analysis_10lag %>% # 
  filter(!Country %in% c(
    "Cook Islands",
    "Democratic People's Republic of Korea",
    "Taiwan"
  )) %>%
  
  left_join(
    gdp_10lag, 
    by = c("Country", "Year")
  ) %>%
  
  left_join(
    gdp_5lag %>% select(-ISO), 
    by = c("Country", "Year")
  ) %>%
  
  mutate(
    log_GDP_10lag = log(GDP_10lag),
    log_GDP_5lag = log(GDP_5lag)
  )


countries_valid <- sbp_gdp_combined %>%
  filter(
    Year == 1980, 
    !is.na(Mean_SBP), 
    !is.na(log_GDP_10lag), 
    !is.na(log_GDP_5lag)
  ) %>%
  pull(Country) %>%
  unique()

# Final Dataset
sbp_gdp_sample <- sbp_gdp_combined %>%
  filter(Country %in% countries_valid) %>%
  mutate(
    Sex = factor(Sex, levels = c("Female", "Male")),
    Year_c = Year - 1980
  )

# Create UN_regions directly in sbp_gdp_sample using UN major regions.

sbp_gdp_sample <- sbp_gdp_sample %>%
  mutate(
    UN_regions = countrycode(
      Country,
      origin = "country.name",
      destination = "un.regionsub.name"
    ),
    
    UN_intermediate = countrycode(
      Country,
      origin = "country.name",
      destination = "un.regionintermediate.name"
    ),
    
    UN_regions_pure22 = coalesce(
      UN_intermediate,
      UN_regions
    )
  )


# Effect coding: Female = -0.5; Male = +0.5
contrasts(sbp_gdp_sample$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

table(sbp_gdp_sample$UN_regions_pure22, useNA = "ifany")#21 regions
table(sbp_gdp_sample$UN_regions, useNA= "ifany")




#=================================================
# Model LOG_GDP 5 years
#=================================================
#log_gdp_5: Full Interactions + Random Slopes and Random Intercepts
M25<- lmer(Mean_SBP ~ Sex * scale(log_GDP_5lag) * Year_c + (1 + Sex + Year_c | Country),
           data = sbp_gdp_sample,
           control = lmerControl(
             optimizer = "bobyqa",
             optCtrl = list(maxfun = 200000)),
           REML= FALSE)
summary(M25)
insight::get_variance(M25)
performance::icc(M25)


lisM25<- allFit(M25)
summary(lisM25)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row25 <- as.numeric(rownames(model.frame(M25)))

vcov_region25 <- vcovCR(M25, 
                        cluster = sbp_gdp_sample$UN_regions_pure22[row25], 
                        type = "CR2")

coef_test(M25, vcov = vcov_region25)



#=====================
# Plot 5 years lag
#=====================

pred_data_1980_lag5 <- expand.grid(
  log_GDP_z_5 = seq(
    min(model.frame(M25)[["scale(log_GDP_5lag)"]],
        na.rm = TRUE),
    max(model.frame(M25)[["scale(log_GDP_5lag)"]],
        na.rm = TRUE),
    length.out = 100
  ),
  Year_c = c(0, 10, 20, 28),
  Sex = c("Female", "Male")
)

pred_data_1980_lag5 <- pred_data_1980_lag5 %>%
  mutate(
    Year = Year_c + 1980,
    Sex = factor(
      Sex,
      levels = c("Female", "Male")
    )
  )

# Correct effect coding

contrasts(pred_data_1980_lag5$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

# Design matrix

X5 <- model.matrix(
  ~ Sex * log_GDP_z_5 * Year_c,
  data = pred_data_1980_lag5
)

# Fixed-effect coefficients

beta5 <- fixef(M25)

colnames(X5) <- gsub("log_GDP_z_5", "scale(log_GDP_5lag)",
                     colnames(X5), fixed = TRUE)
X5 <- X5[, names(beta5), drop = FALSE]

# Cluster-robust variance-covariance matrix

V5 <- as.matrix(vcov_region25)

V5 <- V5[names(beta5), names(beta5), drop = FALSE]

# Predicted values

pred_data_1980_lag5$fit <- as.numeric(
  X5 %*% beta5
)

# Standard errors

pred_data_1980_lag5$se <- sqrt(
  rowSums((X5 %*% V5) * X5)
)

# Confidence intervals

pred_data_1980_lag5 <- pred_data_1980_lag5 %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

# Plot

p5 <- ggplot(
  pred_data_1980_lag5,
  aes(
    x = log_GDP_z_5,
    y = fit,
    color = Sex,
    fill = Sex
  )
) +
  geom_ribbon(
    aes(ymin = lwr, ymax = upr),
    alpha = 0.20,
    color = NA
  ) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~Year) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Log GDP per capita (5-year lag)",
    y = "Predicted SBP (mmHg)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )

print(p5)



#=================================================
# Model LOG_GDP 10 years
#=================================================
#log_gdp_10: Full Interactions + Random Slopes and Random Intercepts
M210<- lmer(Mean_SBP ~ Sex * scale(log_GDP_10lag) * Year_c + (1 + Sex + Year_c | Country),
            data = sbp_gdp_sample,
            control = lmerControl(
              optimizer = "bobyqa",
              optCtrl = list(maxfun = 200000)),
            REML= FALSE)
summary(M210)
insight::get_variance(M210)
performance::icc(M210)

lisM210<- allFit(M210)
summary(lisM210)


#Controlling for potential geographical and cultural non-independence= robust standard errors
row210 <- as.numeric(rownames(model.frame(M210)))

vcov_region210 <- vcovCR(M210, 
                         cluster = sbp_gdp_sample$UN_regions_pure22[row210], 
                         type = "CR2")

coef_test(M210, vcov = vcov_region210)

#=====================
#Plot 10 years lag
#=====================
# Plot predicted SBP by GDP, sex, and year
pred_data_1980_lag10 <- expand.grid(
  log_GDP_z_10 = seq(
    min(scale(sbp_gdp_combined$log_GDP_10lag), na.rm = TRUE),
    max(scale(sbp_gdp_combined$log_GDP_10lag), na.rm = TRUE),
    length.out = 100
  ),
  Year_c = c(0, 10, 20, 28),
  Sex = c("Female", "Male")
)

pred_data_1980_lag10 <- pred_data_1980_lag10 %>%
  mutate(
    Year = Year_c + 1970,
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

contrasts(pred_data_1980_lag10$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

X10 <- model.matrix(
  ~ Sex * log_GDP_z_10 * Year_c,
  data = pred_data_1980_lag10
)

beta10 <- fixef(M210)
pred_data_1980_lag10$fit <- X10 %*% beta10
V10 <- as.matrix(vcov_region210)
pred_data_1980_lag10$se <- sqrt(diag(X10 %*% V10 %*% t(X10)))

pred_data_1980_lag10 <- pred_data_1980_lag10 %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

p10<-ggplot(pred_data_1980_lag10,
            aes(x = log_GDP_z_10,
                y = fit,
                color = Sex,
                fill = Sex)) +
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.20,
              color = NA) +
  
  geom_line(linewidth = 1.3) +
  facet_wrap(~Year) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Standardized Log GDP per capita 10 years lagged",
    y = "Predicted SBP (mmHg)"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )

#Two graphs in one
p10/p5




#==================================================
# Two way fixed effect with log gdp of 5 and 10 years before
#==================================================

mod_fe_base5 <- fixest::feols(
  Mean_SBP ~ Sex * scale(log_GDP_5lag) |
    Country + Year_c,
  data = sbp_gdp_sample,
  cluster = ~ UN_regions_pure22
)

summary(mod_fe_base5)


mod_fe_base10 <- fixest::feols(
  Mean_SBP ~ Sex * scale(log_GDP_10lag) |
    Country + Year_c,
  data = sbp_gdp_sample,
  cluster = ~ UN_regions_pure22
)

summary(mod_fe_base10)

#================================================
#Gender Inequality Paradox
#================================================

library(readr)
GII <- read_csv("GII complete.csv")

GII_2<- GII %>% rename(Country= country)


table(GII_2$Country)
setdiff(sbp_gdp_1980sample$Country, GII_2$Country)



GII_2 <- GII_2 %>%
  mutate(Country = case_when(
    Country == "Korea (Republic of)" ~ "Republic of Korea",
    Country == "Hong Kong, China (SAR)" ~ "China (Hong Kong SAR)",
    Country == "Bolivia (Plurinational State of)" ~ "Bolivia",
    Country == "T\xfcrkiye" | Country == "Türkiye" ~ "Turkey",
    Country == "United States" ~ "United States of America",
    Country == "Congo (Democratic Republic of the)" ~ "Democratic Republic of the Congo",
    Country == "Eswatini (Kingdom of)" ~ "Swaziland",
    Country == "Cabo Verde" ~ "Cape Verde",
    Country == "Timor-Leste" ~ "East Timor",
    Country == "Lao People's Democratic Republic" ~ "Laos",
    Country == "Russian Federation" ~ "Russia",
    Country == "Tanzania (United Republic of)" ~ "Tanzania",
    Country == "Viet Nam" ~ "Vietnam",
    TRUE ~ Country 
  ))


GII_3 <- GII_2 %>%
  semi_join(sbp_gdp_1980sample, by = "Country")

table(GII_3$Country)
View(GII_3)

GII_4<- GII_3[, c(2,6, 16, 24)]
head(GII_4)


library(tidyr)
library(dplyr)

GII_5 <- GII_4 %>%
  pivot_longer(
    cols = starts_with("gii_"),   
    names_to = "Year",         
    names_prefix = "gii_",       
    values_to = "gii"         
  ) %>%
  mutate(Year = as.numeric(Year))  
head(GII_5)

GII_5$Year_c <- GII_5$Year - 1990



sbp_gdp_1980sampleGII<- left_join(sbp_gdp_1980sample, GII_5, by= c("Country", "Year"))



sbp_gdp_1980sample_2GII <- sbp_gdp_1980sampleGII %>%
  filter(Year %in% c(1990, 2000, 2008))
View(sbp_gdp_1980sample_2GII)
sbp_gdp_1980sample_2_smallgii<- sbp_gdp_1980sample_2GII[, c(2, 3, 4, 5, 6, 7, 12, 14, 15, 16, 17, 18, 19)]
View(sbp_gdp_1980sample_2_smallgii)

sbp_gdp_1980sample_2_smallgii$gii_inverted<- 1- sbp_gdp_1980sample_2_smallgii$gii


sbp_gdp_1980sample_2_smallgii$gii_z<- scale(sbp_gdp_1980sample_2_smallgii$gii_inverted)

# Effect coding: Female = -0.5; Male = +0.5
sbp_gdp_1980sample_2_smallgii$Sex <- factor(
  sbp_gdp_1980sample_2_smallgii$Sex,
  levels = c("Female", "Male")
)

contrasts(sbp_gdp_1980sample_2_smallgii$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)



#----------------------------
# Main mixed effect model
#----------------------------
M0_gii <- lmer(Mean_SBP~1 + (1|Country), 
               data= sbp_gdp_1980sample_2_smallgii)
summary(M0_gii)
performance::icc(M0_gii)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row7as <- as.numeric(rownames(model.frame(M0_gii)))

vcov_region7as <- vcovCR(M0_gii, 
                         cluster = sbp_gdp_1980sample_2_smallgii$UN_regions_pure22[row7as], 
                         type = "CR2")
coef_test(M0_gii, vcov = vcov_region7as)

#Main effect model with random intercepts

M1_gii <- lmer(Mean_SBP ~ Sex + gii_z + Year_c.y + (1 |Country),
               data= sbp_gdp_1980sample_2_smallgii,
               control = lmerControl(
                 optimizer = "bobyqa",
                 optCtrl = list(maxfun = 200000)), 
               REML = FALSE)

summary(M1_gii)

performance::icc(M1_gii)

#optimizer check
lism1gii<- allFit(M1_gii)
summary(lism1gii)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row7a <- as.numeric(rownames(model.frame(M1_gii)))

vcov_region7a <- vcovCR(M1_gii, 
                        cluster = sbp_gdp_1980sample_2_smallgii$UN_regions_pure22[row7a], 
                        type = "CR2")
coef_test(M1_gii, vcov = vcov_region7a)

#Full Interaction Effect Model with Random Intercepts

M2_gii <- lmer(
  Mean_SBP ~ Sex * gii_z * Year_c.y + (1 | Country),
  data = sbp_gdp_1980sample_2_smallgii, 
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)

summary(M2_gii)
performance::icc(M2_gii)


#optimizer check
lism2gii<- allFit(M2_gii)
summary(lism2gii)

#Model comparison

anova(M1_gii, M2_gii)



#Controlling for potential geographical and cultural non-independence= robust standard errors
row7 <- as.numeric(rownames(model.frame(M2_gii)))

vcov_region7 <- vcovCR(M2_gii, 
                       cluster = sbp_gdp_1980sample_2_smallgii$UN_regions_pure22[row7], 
                       type = "CR2")
coef_test(M2_gii, vcov = vcov_region7)




#-------------------------------
# plot
#-------------------------------

# Plot predicted SBP by gii_z, sex, and year

pred_data_1990_gii <- expand.grid(
  gii_z = seq(
    min(sbp_gdp_1980sample_2_smallgii$gii_z, na.rm = TRUE),
    max(sbp_gdp_1980sample_2_smallgii$gii_z, na.rm = TRUE),
    length.out = 100
  ),
  Year_c = c(0, 10, 18),
  Sex = c("Female", "Male")
)

pred_data_1990_gii <- pred_data_1990_gii %>%
  mutate(
    Year = Year_c + 1990,
    Year_c.y = Year_c,
    Sex = factor(Sex, levels = c("Female", "Male"))
  )

# Effect coding: Female = -0.5; Male = +0.5

contrasts(pred_data_1990_gii$Sex) <- matrix(
  c(-0.5, 0.5),
  ncol = 1,
  dimnames = list(c("Female", "Male"), "1")
)

X_gii <- model.matrix(
  ~ Sex * gii_z * Year_c.y,
  data = pred_data_1990_gii
)

beta_gii <- fixef(M2_gii)

X_gii <- X_gii[, names(beta_gii), drop = FALSE]

pred_data_1990_gii$fit <- X_gii %*% beta_gii

V_gii <- as.matrix(vcov_region7)

V_gii <- V_gii[names(beta_gii), names(beta_gii), drop = FALSE]

pred_data_1990_gii$se <- sqrt(
  diag(X_gii %*% V_gii %*% t(X_gii))
)

pred_data_1990_gii <- pred_data_1990_gii %>%
  mutate(
    lwr = fit - 1.96 * se,
    upr = fit + 1.96 * se
  )

ggplot(pred_data_1990_gii,
       aes(x = gii_z,
           y = fit,
           color = Sex,
           fill = Sex)) +
  
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),
              alpha = 0.20,
              color = NA) +
  
  geom_line(linewidth = 1.3) +
  
  facet_wrap(~Year) +
  
  theme_minimal(base_size = 14) +
  
  labs(
    x = "Gender Inequality Index (z_scores)",
    y = "Predicted SBP (mmHg)"
  ) +
  
  theme(
    strip.text = element_text(face = "bold"),
    legend.title = element_blank()
  )



#==================================
#Correlation Matrix gii, hdi, gdp
#==================================
library(dplyr)

# ============================================================
# 1. CLEAN DATASETS AND REMOVE DUPLICATES
# ============================================================

# Dataset 1: GDP
ds1_clean <- sbp_gdp_1980sample %>%
  arrange(is.na(GDP_per_capita)) %>%
  distinct(Country, .keep_all = TRUE)

# Dataset 2: HDI
# Prioritize observations with non-missing HDI values
ds2_clean <- sbp_gdp_1980sampleHDI %>%
  arrange(is.na(`Human Development Index`)) %>%
  distinct(Country, .keep_all = TRUE)

# Dataset 3: GII
# Prioritize observations with non-missing GII values
ds3_clean <- sbp_gdp_1980sample_2_smallgii %>%
  arrange(is.na(gii_inverted)) %>%
  distinct(Country, .keep_all = TRUE)



# ============================================================
# 2. MERGE DATASETS BY COUNTRY
# ============================================================

df_joined <- ds1_clean %>%
  inner_join(ds2_clean, by = "Country") %>%
  inner_join(ds3_clean, by = "Country")


# ============================================================
# 3. CHECK MERGED DATA
# ============================================================

# Display column names
colnames(df_joined)

# Number of countries in the merged dataset
n_distinct(df_joined$Country)

# Check HDI distribution
summary(df_joined$`Human Development Index`)

# Count missing HDI values
sum(is.na(df_joined$`Human Development Index`))


# ============================================================
# 4. SELECT VARIABLES FOR CORRELATION ANALYSIS
# ============================================================

# Note: Check column names after joining.
# GDP_per_capita.x may have a different suffix.

df_cor_data <- df_joined %>%
  select(
    GDP  = GDP_per_capita.x,
    HDI  = `Human Development Index`,
    GII  = gii_inverted
  )


# ============================================================
# 5. CHECK MISSING VALUES
# ============================================================

# Count missing values for each variable
colSums(is.na(df_cor_data))

# Summary statistics
summary(df_cor_data)


# ============================================================
# 6. COMPUTE PEARSON CORRELATION MATRIX
# ============================================================

# Pairwise complete observations
cor_matrix <- cor(
  df_cor_data,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# Display correlation matrix
print(round(cor_matrix, 3))


# ============================================================
# 7. COMPUTE INDIVIDUAL CORRELATIONS
# ============================================================

# GDP and HDI
cor(
  df_cor_data$GDP,
  df_cor_data$HDI,
  use = "pairwise.complete.obs",
  method = "pearson"
)

# GDP and GII
cor(
  df_cor_data$GDP,
  df_cor_data$GII,
  use = "pairwise.complete.obs",
  method = "pearson"
)


# HDI and GII
cor(
  df_cor_data$HDI,
  df_cor_data$GII,
  use = "pairwise.complete.obs",
  method = "pearson"
)





# ============================================================
# 8. CHECK NUMBER OF OBSERVATIONS PER CORRELATION
# ============================================================

# Count complete observations for each pair of variables
n_matrix <- outer(
  names(df_cor_data),
  names(df_cor_data),
  Vectorize(function(x, y) {
    sum(complete.cases(df_cor_data[, c(x, y)]))
  })
)

dimnames(n_matrix) <- list(
  names(df_cor_data),
  names(df_cor_data)
)

print(n_matrix)





##########################################################################
#===============================
# Replication Analysis with the NCD.Risc Dataset: https://ncdrisc.org/data-downloads-blood-pressure.html
#===============================
NCD <- read_csv("NCD_RisC_Lancet_2017_BP_age_standardised_countries.csv")

NCD1<- NCD [, c(1,3, 4, 5, 6, 7)]

NCD2 <- NCD1 %>% 
  rename(Country = `Country/Region/World`, 
         Sex= Sex, 
         Year= Year, 
         mean_SBP= `Mean systolic blood pressure (mmHg)`, 
         Lower_CI= `Mean systolic blood pressure lower 95% uncertainty interval (mmHg)`, 
         Upper_CI= `Mean systolic blood pressure upper 95% uncertainty interval (mmHg)`) 

NCD2$Sex <- ifelse(NCD2$Sex=="Men", "Male", "Female")
NCD2$Year<- as.integer(NCD2$Year)

#Merge GDP data with NCD2 data

GDP_data_2<- gdp_data [, c(1, 3, 4)]

NCD_GDP<- NCD2 %>% full_join(GDP_data_2, by= c("Country", "Year"))

NCD_GDP$log_GDP<- log(NCD_GDP$GDP_per_capita)

head(NCD_GDP)

# Restrict analyses to countries with available baseline data in 1975

library(dplyr)



NCD_GDP_selected <- NCD_GDP %>%
  filter(Year %in% c(1975, 1985, 1995, 2005, 2015))

countries_selected <- NCD_GDP_selected %>%
  filter(
    !is.na(mean_SBP),
    !is.na(log_GDP),
    Sex %in% c("Male", "Female")
  ) %>%
  group_by(Country) %>%
  filter(
    n_distinct(Year[Sex == "Male"]) == 5,
    n_distinct(Year[Sex == "Female"]) == 5
  ) %>%
  ungroup()


# Prepare variables for analysis
sbp_gdp_1975sample <- countries_selected %>%
  mutate(
    Year_c = Year - 1975
  )

#Check
Desc(sbp_gdp_1975sample$Country, maxrows = Inf, plotit = F)

# Effect coding: Female = -0.5, Male = 0.5
sbp_gdp_1975sample$Sex<- as.factor(sbp_gdp_1975sample$Sex)
contrasts(sbp_gdp_1975sample$Sex) <- c(-0.5, 0.5)

# Check sample composition by year
sbp_gdp_1975sample %>%
  group_by(Year) %>%
  summarise(
    n_countries = n_distinct(Country),
    n_rows = n(),
    missing_gdp = sum(is.na(GDP_per_capita)),
    .groups = "drop"
  )

# Create UN_regions directly in sbp_gdp_1975sample using UN major regions.




sbp_gdp_1975sample<- sbp_gdp_1975sample %>%
  mutate(
    UN_regions = countrycode(
      Country,
      origin = "country.name",
      destination = "un.regionsub.name"
    ),
    
    UN_intermediate = countrycode(
      Country,
      origin = "country.name",
      destination = "un.regionintermediate.name"
    ),
    
    UN_regions_pure22 = coalesce(
      UN_intermediate,
      UN_regions
    )
  )


table(sbp_gdp_1975sample$UN_regions_pure22, useNA = "ifany")#21 regions



# ---------------------------
# Main mixed-effects analyses
# ---------------------------

#log_gdp should be converted in z_scores (scaled: mean= 0, sd=1) to make the intercept more interpretable.
#If log_gdp is not scaled, the intercept would correspond to: female, 1975 (or zero), when log_gdp= 0, but log_gdp ranges from 5 to 11.
range(sbp_gdp_1975sample$log_GDP, na.rm = TRUE)
mean(sbp_gdp_1975sample$log_GDP)

#If we scaled log_gdp, the intercept would correspond to: female, 1980(or zero), when log_g= 0 (which is the average log_gdp = exp(8.37)= GDP_PER_CAPITA= 4315 USD) 

#Creation of log_GDP_z
sbp_gdp_1975sample$log_GDP_z<- scale(sbp_gdp_1975sample$log_GDP, center = T, scale = T)
range(sbp_gdp_1975sample$log_GDP_z, na.rm = TRUE)


#Step 0: Null Model

M0A<- lmer(mean_SBP ~ 1 + (1 | Country),
          data = sbp_gdp_1975sample, 
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000)),
          REML= FALSE)

summary(M0A)
insight::get_variance(M0A)
performance::icc(M0A)


#Controlling for potential geographical and cultural non-independence= robust standard errors
rowa1x <- as.numeric(rownames(model.frame(M0A)))

vcov_regiona1x <- vcovCR(M0A, 
                        cluster = sbp_gdp_1975sample$UN_regions_pure22[rowa1x], 
                        type = "CR2")


coef_test(M0A, vcov = vcov_regiona1x)


#Step 1: Fixed Main Effect effects + Random Intercepts and Random Slopes
M1A<- lmer(mean_SBP ~ Sex + log_GDP_z + Year_c + (1 + Sex + Year_c | Country),
          data = sbp_gdp_1975sample, 
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000)),
          REML= FALSE)

summary(M1A)



#Controlling for potential geographical and cultural non-independence= robust standard errors
rowa1 <- as.numeric(rownames(model.frame(M1A)))

vcov_regiona1 <- vcovCR(M1A, 
                       cluster = sbp_gdp_1975sample$UN_regions_pure22[rowa1], 
                       type = "CR2")


coef_test(M1A, vcov = vcov_regiona1)


#Step 3: Full Interactions + Random Slopes and Random Intercepts
M2A<- lmer(mean_SBP ~ Sex * log_GDP_z * Year_c + (1 + Sex + Year_c | Country),
          data = sbp_gdp_1975sample,
          control = lmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 200000)),
          REML= FALSE)
summary(M2A)


#Optimizer control
model_list1a<- allFit(M2A)
summary(model_list1a)




#Model comparison
anova(M0A, M1A, M2A)


#Controlling for potential geographical and cultural non-independence= robust standard errors
rowA <- as.numeric(rownames(model.frame(M2A)))

vcov_regionA <- vcovCR(M2A, 
                      cluster = sbp_gdp_1975sample$UN_regions_pure22[rowA], 
                      type = "CR2")


coef_test(M2A, vcov = vcov_regionA)
car::vif(M2A)



# ---------------------------
# Robustness analyses
# ---------------------------

# AR(1) model for temporal autocorrelation
sbp_gdp_1975sample <- sbp_gdp_1975sample %>%
  mutate(
    Country = factor(Country)
  )

model_ar1_1975sample <- lme(
  mean_SBP ~ Sex * log_GDP_z * Year_c,
  random = ~ 1 + Sex + Year_c | Country,
  correlation = corAR1(form = ~ Year_c | Country/Sex),
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  data = sbp_gdp_1975sample,
  method = "ML"
)
summary(model_ar1_1975sample)

# Sex-specific AR(1) robustness checks
sbp_male1975 <- sbp_gdp_1975sample %>%
  filter(Sex == "Male") %>%
  droplevels()

sbp_female1975 <- sbp_gdp_1975sample %>%
  filter(Sex == "Female") %>%
  droplevels()

# Male model with AR(1)
model_ar1_male1975 <- lme(
  mean_SBP ~ log_GDP_z * Year_c,
  random = ~ 1  + Year_c  | Country,
  correlation = corAR1(form = ~ Year_c | Country),
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  data = sbp_male1975,
  method = "ML",
  na.action = na.omit
)

summary(model_ar1_male1975)

# Female model with AR(1)
model_ar1_female1975 <- lme(
  mean_SBP ~ log_GDP_z * Year_c,
  random = ~ 1 + Year_c | Country,
  correlation = corAR1(form = ~ Year_c | Country),
  data = sbp_female1975,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  method = "ML",
  na.action = na.omit
)

summary(model_ar1_female1975)

#--------------------------
# Spatial Autocorrelation of Residuals: Moran's I
#--------------------------

library(dplyr)
library(spdep)
library(lme4)

# 1. M2A residuals
dati_modello2 <- model.frame(M2A)

dati_modello2$residuo2 <- residuals(M2A)



lookup_regioni <- sbp_gdp_1975sample %>%
  select(Country, UN_regions_pure22) %>%
  mutate(Country = as.character(Country)) %>%
  distinct()


stopifnot(!anyDuplicated(lookup_regioni$Country))

dati_modello2 <- dati_modello2 %>%
  mutate(Country = as.character(Country)) %>%
  left_join(lookup_regioni, by = "Country")

stopifnot(!anyNA(dati_modello2$UN_regions_pure22))

# 3. average residual for each country 
residui_country2 <- dati_modello2 %>%
  group_by(Country, UN_regions_pure22) %>%
  summarise(
    residuo = mean(residuo2, na.rm = TRUE),
    .groups = "drop"
  )

# 4. UN macro-regions cluster matrix
cluster2 <- as.character(
  residui_country2$UN_regions_pure22
)


W2 <- outer(cluster2, cluster2, "==") * 1


diag(W2) <- 0


listw2 <- spdep::mat2listw(
  W2,
  style = "W",
  zero.policy = TRUE
)

# 6. Global Moran's I
moran_result2 <- spdep::moran.test(
  residui_country2$residuo,
  listw2,
  alternative = "two.sided",
  zero.policy = TRUE
)

print(moran_result2)




#--------------------------
# Longitudinal MLM 3 way interaction + variance decomposition LogGDP (Mundlak decomposition)
#--------------------------
#Let us create a new data-frame where Log_gdp_z is decomposed within and between country


# One observation per country-year
gdp_country_year1975 <- sbp_gdp_1975sample %>%
  distinct(Country, Year_c, log_GDP_z)

# Ensure GDP is unique for each country-year
stopifnot(
  nrow(gdp_country_year1975) ==
    nrow(distinct(gdp_country_year1975, Country, Year_c))
)

# Decompose into within and between components
gdp_country_year1975 <- gdp_country_year1975 %>%
  group_by(Country) %>%
  mutate(
    gdp_between = mean(log_GDP_z, na.rm = TRUE),
    gdp_within  = log_GDP_z - gdp_between
  ) %>%
  ungroup()

# Merge into the original dataset
sbp_gdp_1975sample <- sbp_gdp_1975sample %>%
  select(-any_of(c("gdp_within", "gdp_between",
                   "gdp_within_std", "gdp_between_std"))) %>%
  left_join(
    gdp_country_year1975 %>%
      select(Country, Year_c, gdp_within, gdp_between),
    by = c("Country", "Year_c")
  )
#--------------------------------
#Longitudinal Multilevel Model 3 way interaction with log_gdp_c variance decomposition + Random Intercept and Random Slopes
#--------------------------------
#Null Model
mod_mundlak_0A <- lmer(
  mean_SBP ~ 1 +
    (1 |  Country),
  data = sbp_gdp_1975sample,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)

summary(mod_mundlak_0A)
performance::icc(mod_mundlak_0A)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row2a1z <- as.numeric(rownames(model.frame(mod_mundlak_0A)))

vcov_region2a1z <- vcovCR(mod_mundlak_0A, 
                         cluster = sbp_gdp_1975sample$UN_regions_pure22[row2a1z], 
                         type = "CR2")

coef_test(mod_mundlak_0A, vcov = vcov_region2a1z)

#Main Fixed Effect Model
mod_mundlak_1A <- lmer(
  mean_SBP ~ Sex + Year_c + (gdp_within + gdp_between) +
    (1  + Sex + Year_c | Country),
  data = sbp_gdp_1975sample,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)



summary(mod_mundlak_1A)



#Controlling for potential geographical and cultural non-independence= robust standard errors
row2a1 <- as.numeric(rownames(model.frame(mod_mundlak_1A)))

vcov_region2a1 <- vcovCR(mod_mundlak_1A, 
                        cluster = sbp_gdp_1975sample$UN_regions_pure22[row2a1], 
                        type = "CR2")

coef_test(mod_mundlak_1A, vcov = vcov_region2a1)

#Optimizer Check
lismod1a<-allFit(mod_mundlak_1A)
summary(lismod1a)

#Full Interaction Effect Model
mod_mundlak_2A <- lmer(
  mean_SBP ~ Sex * Year_c * (gdp_within + gdp_between) +
    (1 + Sex + Year_c  | Country),
  data = sbp_gdp_1975sample,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 200000)),
  REML = FALSE
)

summary(mod_mundlak_2A)


#Optimizer Check
lismod2a<-allFit(mod_mundlak_2A)
summary(lismod2a)
#---------------------------------
# Mundlak Models Comparison
#---------------------------------
anova(mod_mundlak_0A, mod_mundlak_1A, mod_mundlak_2A)

#Controlling for potential geographical and cultural non-independence= robust standard errors
row2aa <- as.numeric(rownames(model.frame(mod_mundlak_2A)))

vcov_region2aa <- vcovCR(mod_mundlak_2A, 
                       cluster = sbp_gdp_1975sample$UN_regions_pure22[row2aa], 
                       type = "CR2")

coef_test(mod_mundlak_2A, vcov = vcov_region2aa)

car::vif(mod_mundlak_2A)

cor(sbp_gdp_1975sample$gdp_within, sbp_gdp_1975sample$Year_c)



#---------------------------------
# Fixed-effect regression (fixed effects: country + time) 
#One way clustered errors: it allows standard error country correlations and heteroskedasticity 
#---------------------------------

mod_fe_basencd <- fixest::feols(
  mean_SBP ~ Sex*log_GDP_z |
    Country + Year_c,
  data = sbp_gdp_1975sample,
  cluster =  ~ UN_regions_pure22, 
)

summary(mod_fe_basencd)



#==================================
# End of the script!
#=================================