# Results: Output results for manuscript 
# Author: Nadine Zakkak

library(dplyr)
library(tidyr)
library(flextable)

source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
path_save <- "./tables_figs/"

# this data will be needed several times - import it once in the beginngin
source("./scripts/00_cancer_stage_dataprep.R") #import data from MySQL


#Get info about age at highest risk by symptom, cancer and sex -----
# Data from 4_bootstrap_analysis.R
# Advanced stage:
advanced_predict <- read.csv(file = "./results/bootstrap_10000/advanced_predict_bs.csv")
early_predict <- read.csv(file = "./results/bootstrap_10000/early_predict_bs.csv")

final_predict <- advanced_predict |> mutate(stage='advanced') |> bind_rows(early_predict|>mutate(stage='early'))

max_prob <- 
  final_predict |>
  group_by(stage, gender, cancer_site, index_labels) |>
  filter(prob == max(prob)) |>
  mutate(prob = 100*prob, conf.low = 100*conf.low, conf.high = 100*conf.high) |>
  arrange(gender, stage, cancer_site, index_labels) |>
  select(gender, stage, cancer_site, index_labels, age, prob, conf.low, conf.high)

max_prob |>
  mutate(across(where(is.numeric), format_numb)) |>
  flextable() |>
  merge_v(j=c(1,2,3)) |>
  hline(i=c(2,4,6,8,10,12,14)) |>
  fontsize(part = "all", size = 9) |>
  bold(part = "header") |>
  font(part = "all", fontname =  "Calibri") |>
  save_as_docx(path = paste0(path_save, "max_abs_risk.docx"))

# Stage cohort proportion with both symptoms -----
# Replica of code in 2b_stage_descriptive.R
prop_symptom <- prop.table(table(df$sx))*100
write.csv(prop_symptom, paste0(path_save, "prop_symptom.csv"), row.names = F)

# Stage cohort proportion with missing stage -----
# Replica of code in 2b_stage_descriptive.R
prop_stage <- df |> group_by(stage_best) |> summarise(n=n()) |> mutate(N=sum(n), prop = format_numb(100*n/N))
write.csv(prop_stage, paste0(path_save, "prop_stage.csv"), row.names = F)

# Stage cohort proportion with missing stage by age (25-69, 70+) -----
# Replica of code in 2b_stage_descriptive.R
prop_stage_age <- 
  df |>
  mutate(age70_bin = case_when(
    age10_cat %in% c("30-39", "40-49", "50-59", "60-69") ~ "<70",
    age10_cat %in% c("70-79", "80-89", "90+") ~ ">=70")) |>
  group_by(age70_bin, stage_best) |> summarise(n=n()) |> mutate(N = sum(n), prop = format_numb(100*n/N))
write.csv(prop_stage_age, paste0(path_save, "prop_stage_by_age.csv"), row.names = F)

# Stage cohort proportion of advanced stage by cancer site -----
prop_stage_site <- df_cc |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup() |>
  group_by(cancer_site, advanced) |>
  summarise(n=n()) |> mutate(prop = 100*n/sum(n)) |>
  filter(advanced)
write.csv(prop_stage_site, paste0(path_save, "prop_stage_by_site.csv"), row.names = F)

# Stage cohort proportion of symptoms by cancer site -----
prop_sx_site <- df_cc |>
  left_join(df_cc |> group_by(cr_eventtype) |> summarise(N = n())) |>
  mutate(Neither = sx_eventtype_14 == 0 & sx_eventtype_4 == 0,
         `Rectal Bleeding` = as.logical(as.numeric(sx_eventtype_14) - 1),
         CIBH = as.logical(as.numeric(sx_eventtype_4) -1)) |>
  pivot_longer(cols = c(`Rectal Bleeding`, CIBH, Neither), names_to = "symptom") |>
  group_by(cr_eventtype, N, symptom) |>
  summarise(n = sum(value)) |>
  mutate(prop = format_numb(100*n/N)) |>
  ungroup() |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup() |>
  select(-cr_eventtype) |>
  select(cancer_site, N, symptom, n, prop)
write.csv(prop_sx_site, paste0(path_save, "prop_sx_by_site.csv"), row.names = F)


# Risk of cancer at 80 years old -----
# Data from 3a_cancer_incidence_analysis.R
risk_cancer_preds <- readRDS(file="./results/cancer_incidence_cohort/allspline_pred.rds")

risk_80_cancer <- 
  risk_cancer_preds |>
  group_by(gender, cancer_site, index_labels) |>
  filter(age>80) |>
  mutate(age=round(age), predicted = 100*predicted, conf.low = 100*conf.low, conf.high = 100*conf.high) |>
  filter(age==80) |>
  select(age, index_labels, cancer_site, gender, predicted, conf.low, conf.high) |>
  arrange(index_labels, cancer_site, gender) 
write.csv(risk_80_cancer, "./results/cancer_incidence_cohort/risk_80.csv")
write.csv(risk_80_cancer, paste0(path_save, "risk_cancer_80.csv"), row.names = F)


# Poportion of advanced stage cancer in cases at 60 years old -----
# Data from 3b_stage_analysis.R
load(file="./results/cancer_stage_cohort/allmodels_spline_pred.RData")
risk_advanced_preds <- allmodels_spline_pred; rm(allmodels_spline_pred)

risk_60_advanced <- 
  risk_advanced_preds |>
  group_by(gender, cancer_site, sx) |>
  filter(age>60) |>
  mutate(age=round(age), predicted=100*predicted, conf.low = 100*conf.low, conf.high = 100*conf.high) |>
  filter(age==60) |>
  select(age, cancer_site, sx, gender, predicted, conf.low, conf.high) |>
  arrange(cancer_site, sx, gender)
write.csv(risk_60_advanced, "./results/cancer_stage_cohort/risk_60.csv")
write.csv(risk_60_advanced, paste0(path_save, "risk_adv_60.csv"), row.names = F)

# Sensitivity analysis (1) excluding screen detected cancers descriptives ------
## number of cases detected via screening ----
prop_route <- df_cc |>
  group_by(final_route) |>
  summarise(n=n()) |>
  mutate(N=sum(n), prop=100*n/N) |>
  ungroup()
write.csv(prop_route, paste0(path_save, "prop_route.csv"), row.names = F)

## Descriptives of patients that remain after excluding screen detected cancers -----
df_noscreen <- df_cc |> mutate(screen = final_route == "Screening") |> filter(!screen)
# Stage cohort proportion of advanced stage by cancer site -----
noscreen_prop_stage_site <- df_noscreen |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup() |>
  group_by(cancer_site, advanced) |>
  summarise(n=n()) |> mutate(N=sum(n), prop = 100*n/sum(n)) |>
  filter(advanced) |>
  select(cancer_site, advanced, N, n, prop)
write.csv(noscreen_prop_stage_site, paste0(path_save, "noscreen_prop_stage_site.csv"), row.names = F)

# Stage cohort proportion of symptoms by cancer site -----
noscreen_prop_sx_site <- df_noscreen |>
  left_join(df_noscreen |> group_by(cr_eventtype) |> summarise(N = n())) |>
  mutate(Neither = sx_eventtype_14 == 0 & sx_eventtype_4 == 0,
         `Rectal Bleeding` = as.logical(as.numeric(sx_eventtype_14) - 1),
         CIBH = as.logical(as.numeric(sx_eventtype_4) -1)) |>
  tidyr::pivot_longer(cols = c(`Rectal Bleeding`, CIBH, Neither), names_to = "symptom") |>
  group_by(cr_eventtype, N, symptom) |>
  summarise(n = sum(value)) |>
  mutate(prop = format_numb(100*n/N)) |>
  ungroup() |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup() |>
  select(-cr_eventtype) |>
  select(cancer_site, N, symptom, n, prop)
write.csv(noscreen_prop_sx_site, paste0(path_save, "noscreen_prop_sx_by_site.csv"), row.names = F)


# Symptom cohort no screen description ----
source("./scripts/00_cancer_indicidence_dataprep.R") #import data from MySQL
df_screen <- df |> filter(first_crc_route == "Screening")
# variables to stratify by 
categ <- c(quo(gender), quo(age10_cat), quo(imd2015_5), quo(index_eventtype))
desc_summary <-
  lapply(categ, FUN = function(x){ # summary per covariate
    df_screen |>
      mutate(!!x := as.character(!!x)) |>
      group_by(!!x, cancer_site) |>
      summarise(n = n()) |>
      mutate(N = sum(n),
             prop = n/N,
             prop100 = format_numb(100*prop),
             lower = calculate_ci(n, N)$lower,
             upper = calculate_ci(n, N)$upper
      ) |>
      mutate(!!x := as.factor(!!x)) |>
      mutate(prop100 = ifelse(n<=6, "-", prop100),
             n = ifelse(n<=6, "n<=6", as.character(n)),
             N = ifelse(N<=6, "N<=6", as.character(N)),
             summary = sprintf("%s (%s%%)", n, prop100)
      ) |>
      ungroup()
  })

desc_summary_df <- bind_rows(desc_summary) |>
  bind_rows( # add info of total patients
    df_screen |> group_by(cancer_site) |>
      summarise(n = n()) |>
      mutate(N = sum(n),
             prop = n/N,
             prop100 = format_numb(100*prop),
             lower = calculate_ci(n, N)$lower,
             upper = calculate_ci(n, N)$upper,
             Total = "Total") |>
      mutate(prop100 = ifelse(n<=6, "-", prop100),
             n = ifelse(n<=6, "n<=6", as.character(n)),
             N = ifelse(N<=6, "N<=6", as.character(N)),
             summary = sprintf("%s (%s%%)", n, prop100)
             
      ) |>
      ungroup()
  ) |>
  mutate(across(!n & !N & !summary, as.character)) |>
  recode_gender() |>
  mutate(imd2015_5 = ifelse(imd2015_5 == 1, "1 - Least", ifelse(imd2015_5 == 5, "5 - Most", imd2015_5))) |>
  rowwise() |>
  mutate(index_eventtype = ifelse(!is.na(index_eventtype), 
                                  names(symptom_codes)[symptom_codes == index_eventtype], NA),
         cancer_site = ifelse(!is.na(cancer_site), names(cancer_codes)[cancer_codes == cancer_site], NA))|>
  ungroup() |>
  rename(Gender = gender,
         Age = age10_cat,
         `IMD (2015)` = imd2015_5,
         `Index Symptom` = index_eventtype) |>
  pivot_longer(cols = c(Age, 
                        # Age_cont, 
                        `IMD (2015)`, `Index Symptom`, Gender, Total), names_to = "covariate") |> 
  filter(!is.na(value)) 

desc_summary_df_wide <- desc_summary_df |>
  pivot_wider(id_cols = c(covariate, value, N), names_from = cancer_site, 
              values_from = c(summary), names_prefix = "summary_")

desc_summary_df_wide |>
  select(covariate, value, N, summary_Colon, summary_Rectal) |>
  rename_with(.fn = ~gsub("summary_", "", .x), cols=starts_with("summary_")) |>
  flextable() |>
  merge_v(j=1) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "desc_cancer_incidence_screen.docx"))

