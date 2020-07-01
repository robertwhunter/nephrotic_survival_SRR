################################################################
################################################################
####
#### Code used in Kolb et al. (2020)...
#### ...Death and ESKD in adult nephrotic sydrnome
####
#### Survival analysis and competing risk models
#### Estimating mortality in age- and sex-matched population
####
################################################################
################################################################


#### LOAD PACKAGES ----
library(tidyverse)
library(lubridate)
library(stringr)
library(survminer)
library(survival)
library(cmprsk)
library(riskRegression)
library(prodlim)


#### DATA SOURCES ----

# MAIN DATASET = df_all
#
# data for all patients in the study


# EVENT DATA SUBSET = df_events
#
# 

# code first event (ESRF or death)
df_events %>% mutate(
  first_event = case_when(
    outcome == "Death_without_ESRF" ~ "death",
    outcome == "ESRF_then_death" ~ "ESRF",
    outcome == "ESRF_without_death" ~ "ESRF",
    outcome == "Neither_ESRF_nor_death" ~ "censor"
  ),
  time_to_event = case_when(
    first_event == "death" | first_event == "censor" ~ survival_without_death,
    first_event == "ESRF" ~ survival_without_ESRF_days
  )
) -> df_events

# re-level so that censor is top level (pre-requisite for competing risks model)
df_events$first_event %>% 
  as.factor() %>% 
  fct_relevel("censor", "ESRF", "death") -> df_events$first_event

# filter out patients with ESRF at baseline
df_events %>% 
  filter((first_event == "ESRF" & time_to_event <= 0) == FALSE) -> df_events

# select only older patients with primary NS for subgroup analysis
df_events_primold <- df_events %>% filter(classification == "primary", Over_60 == TRUE)



#### SURVIVAL & COXPH ANALYSIS ----

#### Censoring model (i.e not competing risks)
model_all <- coxph(Surv(time = survival_days, event = Died) ~Age + Sex + classification + Urgency + alb_cut + CKD_EPI + Hb + SIMD16_Quintile, 
      data = df_all)

fitmort_all <- survfit(
  Surv(time = survival_days, event = Died) ~classification + Over_60, 
  data = df_all
)


#### COMPETING RISKS ANALYSIS FOR ESRF AND COMPETING RISK OF DEATH (Fine-Gray) ----
# 
# use three alternative methods as a sense-check
# all three methods gave similar results
# results from method 1 in paper
#
#

# Conventional model without competing risk (i.e. censoring for death) - NOT USED
censor_fit_ESRF <- survfit(
  Surv(time = survival_without_ESRF_days, event = ESRF) ~classification + Over_60,
  data = df_events
)


#### CR MODEL - method 1 (finegray function)
model_formula <- Surv(time = time_to_event, event = first_event) ~ .

df_finegray_primold <- finegray(
  formula = model_formula,
  data = df_events_primold,
  etype = "ESRF"
)

coxph(Surv(fgstart, fgstop, fgstatus) ~Age + Sex + Diagnosis + alb_cut + CKD_EPI_0m + Hb0m + SIMD16_Quintile, 
  data = df_finegray_primold,
  weights = fgwt,
  na.action = "na.omit") -> model_CR_primold


#### CR MODEL - alternative method 2 (crr function)
#
#
# first build a model matrix for the co-variates
# use na.pass to keep NAs so that number of rows stays the same - else error in crr()
options(na.action = "na.pass")
model.matrix(~ Age + Sex + Diagnosis + alb_cut + CKD_EPI_0m + Hb0m + SIMD16_Quintile, 
             data = df_events_primold)[,-1] -> cov_primold

# call crr to perform regression, modelling the hazard of the cumulative incidence function (Fine & Gray) 
crr(ftime = df_events_primold$time_to_event,
    fstatus = df_events_primold$first_event,
    cov1 = cov_primold,
    failcode = "ESRF",
    cencode = "censor"
) -> model_CR_primold_cmprsk


#### CR MODEL - alternative method 3 (FGR function)
FGR(Hist(time_to_event, first_event) ~ Age + Sex + Diagnosis + alb_cut + CKD_EPI_0m + Hb0m + SIMD16_Quintile, 
    data = df_events_primold, 
    cause = "ESRF") -> model_CR_primold_FGR


#### Event curves & cumulative incidence functions
finegray(
  Surv(time = time_to_event, event = first_event) ~ . + strata(classification, Over_60),
  data = df_events,
  etype = "ESRF"
) -> FG_data_ESRF

survfit(
  Surv(fgstart, fgstop, fgstatus) ~ classification + Over_60,
  data = FG_data_ESRF, 
  weight = fgwt
) -> FG_fit_ESRF

cuminc(ftime = df_events$time_to_event,
       fstatus = df_events$first_event,
       cencode = "censor",
       group = df_events$classification
) -> CIF_all_classification

cuminc(ftime = df_events$time_to_event,
       fstatus = df_events$first_event,
       cencode = "censor",
       group = df_events$Over_60
) -> CIF_all_age


#### Estimating age- and sex-matched survival in the general population ----

#### get NRS data ----
# download from...
# https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/life-expectancy/life-expectancy-at-scotland-level/scottish-national-life-tables/2015-2017/national-life-tables

t <- read.csv("Data - input/NRS/NRS_by_age.csv", skip = 4, nrows = 101)
t <- t %>% select("Age" = x, "Male_qx" = qx, "Female_qx" = qx.1)


# predicted 3-yr mortality approximated as 1 - [(1 - qx)^3], where qx is the age- and sex-specific one-year probability of death at the time of biopsy.  

# compute 3 yr mortality from 1 yr mortality (assuming no change in qx over the 3 yrs)
x_yr_mort <- function(yrs, qx) {
  p_not_dying_each_yr <- 1-qx
  p_not_dying <- p_not_dying_each_yr ^ yrs
  p_dying <- 1 - p_not_dying
  p_dying %>% return()
}

t %>% mutate(
  `3yr_mort_male` = x_yr_mort(yrs = 3, qx = Male_qx),
  `3yr_mort_female` = x_yr_mort(yrs = 3, qx = Female_qx)
) -> t

df_male <- t %>% select(Age, "3yr_mort" = `3yr_mort_male`)
df_female <- t %>% select(Age, "3yr_mort" = `3yr_mort_female`)
df_male$Sex = c("M")
df_female$Sex = c("F")
df_life_exp <- rbind(df_male, df_female)

# merge into patient data
df_short <- df_all %>% select(Sex, Age, classification, Over_60)
df_short$Age <- round(df_short$Age)
df_short <- left_join(df_short, df_life_exp, by = c("Sex", "Age"))
