# Realistic worked example in a Singaporean emergency department
options(scipen = 100, digits = 5)
library(tidyverse)

# Import data
df_ed_raw <- readRDS("C:/Users/Robin/NUS Dropbox/EDData/ED2A/ED2A_Blythe Robin Daniel/data_all_temp2022.11.18.RDS")

df_ed <- df_ed_raw |>
  filter(REGISTRATION_DATE >= "2017-01-01") |>
  mutate(ID = row_number(),
         outcome_icu = ifelse(n_icu_overall_current > 0, 1, 0),
         outcome_died_30d = ifelse(outcome_mortality_30d == TRUE, 1, 0)) |>
  select(ID, AGE, PULSE, RESPIRATION, BPSYSTOLIC, BPDIASTOLIC, AllCancer, outcome_admit, outcome_icu, outcome_died_30d) |>
  rowwise() |>
  mutate(pred_risk = sum(
    case_when(AGE < 30 ~ 0,
              AGE >= 30 & AGE <= 49 ~ 8,
              AGE >= 50 & AGE <= 79 ~ 14,
              AGE >= 80 ~ 19),
    case_when(PULSE < 60 ~ 1,
              PULSE >= 60 & PULSE <= 69 ~ 0,
              PULSE >= 70 & PULSE <= 94 ~ 2,
              PULSE >= 95 & PULSE <= 109 ~ 6,
              PULSE >= 110 ~ 9),
    case_when(RESPIRATION < 16 ~ 8,
              RESPIRATION >= 16 & RESPIRATION <= 19 ~ 0,
              RESPIRATION >= 20 ~ 6),
    case_when(BPSYSTOLIC < 100 ~ 8,
              BPSYSTOLIC >= 100 & BPSYSTOLIC <= 114 ~ 5,
              BPSYSTOLIC >= 115 & BPSYSTOLIC <= 149 ~ 2,
              BPSYSTOLIC >= 150 ~ 0),
    case_when(BPDIASTOLIC < 50 ~ 3,
              BPDIASTOLIC >= 50 & BPDIASTOLIC <= 94 ~ 0,
              BPDIASTOLIC >= 95 ~ 2),
    case_when(AllCancer == "0no" ~ 0,
              AllCancer == "1local" ~ 6,
              AllCancer == "2metastatic" ~ 14)
    )) |>
  ungroup() |>
  na.omit()

saveRDS(df_ed, file = "C:/Users/Robin/NUS Dropbox/EDData/ED2A/ED2A_Blythe Robin Daniel/dat_scored.RDS")
remove(df_ed_raw)
