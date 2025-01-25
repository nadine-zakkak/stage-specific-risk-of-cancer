# Results: Output all required figures and tables in Appendix 
# Author: Nadine Zakkak

library(dplyr)
library(tidyr)
library(flextable)
library(ggplot2)

source("./scripts/00_functions.R")
source("./scripts/00_global_var.R")
path_save <- "./tables_figs/"
save_image_var_one_and_half <- list(device = "jpeg", unit = "mm", width = 140, height = 120, dpi=1000)
save_image_var_full <- list(device = "jpeg", unit = "mm", width = 190, height = 150, dpi=1000)


# Table A.1 ----
# Comparison of characteristics based on cohort origin (symptom cohort vs random sample)
# Data from 2b_stage_descriptive.R
compare_cohort <- read.csv("./results/cancer_stage_cohort/compare_cohort_desc.csv")

compare_cohort |>
  mutate(p.value = round(p.value, 3)) |>
  flextable() |>
  set_header_labels(covariate="", value="",
                    random.sample = "Cases nested within random sample",
                    symptom.cohort = "Cases nested within symptom cohort",
                    p.value = "P value") |>
  footnote(i=1, j=5, ref_symbols = "*", value = as_paragraph("Wilcoxon-rank sum test")) |>
  footnote(i=2:14, j=5, ref_symbols = "+", value = as_paragraph("t-test")) |>
  merge_v(c(1,5)) |>
  hline(i=c(1,3,8,10,12,14)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "TableA1.docx"))

# Figure A.1 ----
# Descriptive of cancer stage cohort incl. missing stages
# Data from 2b_stage_descriptive.R 
stage_desc <- read.csv("./results/cancer_stage_cohort/raw_desc_stage_all.csv")

total_label <-
  stage_desc |>
  group_by(cancer_site, covariate, value) |>
  summarise(N = sum(n)) |>
  mutate(label = paste0("N=",N)) |>
  ungroup()

stage_desc <- 
  stage_desc |>
  mutate(stage_best = case_when(
    stage_best == 1 ~ "I", 
    stage_best == 2 ~ "II",
    stage_best == 3 ~ "III",
    stage_best == 4 ~ "IV",
    stage_best == -99 ~ "Unknown"),
    stage_best = factor(stage_best, levels = rev(c("I", "II", "III", "IV", "Unknown"))),
    value = factor(value, levels = c("<50", "50-59", "60-69", "70-79", "80+",
                                     "Men", "Women",
                                     "1 - Least", "2", "3", "4", "5 - Most",
                                     "Rectal Bleeding", "CIBH", "Neither",
                                     "Total"))) |>
  ungroup() |>
  select(covariate, value, cancer_site, prop, lower, upper, stage_best, N)

figa1 <- 
  ggplot() +
  geom_bar(aes(x = value, y = prop, fill = stage_best), stat = "identity", data = stage_desc) +
  ylim(-13, 101)+
  geom_text(size = 2, hjust = 0, colour = "#5A5A5A", aes(x = value, y = -13, label = label), data = total_label) +
  scale_fill_manual(values = rev(c("#5340ff", "lightblue", "#C47B7B", "darkred", "grey"))) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  facet_grid(covariate~cancer_site, scales = "free_y") +
  labs(x = "",
       y = "Proportion of Patients (%)") +
  guides(fill=guide_legend(reverse=T, title="Stage")) +
  theme_minimal()+
  theme(strip.text = element_text(size = 7, colour = "black"),
        axis.text.x = element_text(size = 6, colour = "black"),
        axis.text.y = element_text(size = 6, colour = "black"),
        legend.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 7, colour = "black"),
        axis.title.x = element_text(size = 7, colour =  "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())
do.call(ggsave, c(list(plot = figa1, filename = paste0(path_save, "FigA1.jpeg")), save_image_var_full))

# Figure A2 -----
# Risk of advanced stage cancer within cancer cases from diff data sources (study extract and CR 2018-2019)
# Replica of script 9c_stage_cprd_cr_analysis.R
# Data from 3b_stage_analysis.R and 9b_analyse_cr_data.R
load(file="./results/cancer_stage_cohort/allmodels_spline_pred.RData")
cr_model_preds <- readRDS("cr_model_preds.rds")

cr_model_preds <-
  cr_model_preds |>
  mutate(data_source = "Cancer registry data (2018-2019)",
         gender = ifelse(male==0, "Women", "Men")) |>
  rename(cancer_site = cancer) |>
  ungroup() |>
  select(cancer_site, gender, age, pr, pr_lb, pr_ub, data_source)

allmodels_spline_pred <- 
  allmodels_spline_pred |>
  mutate(predicted = 100*predicted, conf.low = 100*conf.low, conf.high = 100*conf.high, 
         sx = factor(sx, levels = c("Rectal Bleeding", "CIBH", "Neither")),
         data_source = "Study extract") |>
  select(cancer_site, gender, age, pr=predicted, pr_lb=conf.low, pr_ub=conf.high, sx, data_source)

figa2 <- ggplot() +
  geom_line(aes(x = age, y = pr, color = sx, linetype = data_source), 
            linewidth = linesize, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = pr_lb, ymax = pr_ub, fill = sx),
              alpha = .1, data = allmodels_spline_pred) +
  geom_line(aes(x = age, y = pr, linetype = data_source),
            linewidth = linesize, data = cr_model_preds) +
  geom_ribbon(aes(x = age, ymin = pr_lb, ymax = pr_ub),
              alpha = .1, data = cr_model_preds) +
  facet_grid(gender~cancer_site) +
  scale_fill_manual(
    breaks = c("Rectal Bleeding", "CIBH", "Neither"),
    values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_color_manual(breaks = c("Rectal Bleeding", "CIBH", "Neither"),
                     values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_linetype_manual(values=c("dashed", "solid") ,"Data source") +
  legend +
  theme_light() +
  theme +
  theme(text=element_text(size=7)) +
  labs(y = "Probability of advanced stage cancer in cases (%)",
       x = x_label)
do.call(ggsave, c(list(plot = figa2, filename = paste0(path_save, "FigA2.jpeg")), save_image_var_one_and_half))

# Figure A3 ------
# Sensitivity analysis (1): Risk of advanced cancer in cancer cases predictions excluding screen-detected cancer
# Data from 5_sensitivity_analysis.R
sensitivity_excl_screen_preds <-  read.csv("./results/sensitivity_analysis/excl_screen/allmodels_spline_pred.csv") |>
  rename(index_labels = sx) |>
  mutate(index_labels = factor(index_labels, c("Rectal Bleeding", "CIBH", "Neither")))

figa3 <- plot_preds_data(sensitivity_excl_screen_preds, "Probability of advanced stage cancer in cases (%)")
do.call(ggsave, c(list(plot = figa3, filename = paste0(path_save, "FigA3.jpeg")), save_image_var_one_and_half))


# Figure A4 ------
# Sensitivity analysis (2): Risk of advanced cancer (IV) in cancer cases predictions
# Data from 5_sensitivity_analysis.R
sensitivity_stage4_preds <-  read.csv("./results/sensitivity_analysis/stage/allmodels_spline_pred.csv") |>
  rename(index_labels = sx) |>
  mutate(index_labels = factor(index_labels, c("Rectal Bleeding", "CIBH", "Neither")))

figa4 <- plot_preds_data(sensitivity_stage4_preds, "Probability of advanced stage cancer in cases (%)")
do.call(ggsave, c(list(plot = figa4, filename = paste0(path_save, "FigA4.jpeg")), save_image_var_one_and_half))


# Table A.2 ----
get_pred_at_age <- function(data, target_age, group_cols){
  data |>
    group_by(!!!syms(group_cols)) |>
    filter(age>target_age) |>
    mutate(age=round(age), predicted = 100*predicted, conf.low = 100*conf.low, conf.high = 100*conf.high) |>
    filter(age==target_age) 
}

cancer_preds <- readRDS(file = "./results/cancer_incidence_cohort/allspline_pred.rds") |>
  mutate(index_labels = factor(index_labels, levels = c("Rectal Bleeding", "CIBH")))

group_cols <- c("gender", "cancer_site", "index_labels")
cancer_preds_40 <- get_pred_at_age(cancer_preds, 40, group_cols)
cancer_preds_60 <- get_pred_at_age(cancer_preds, 60, group_cols)  
cancer_preds_80 <- get_pred_at_age(cancer_preds, 80, group_cols)

cancer_preds_40to80 <- rbind(cancer_preds_40, cancer_preds_60, cancer_preds_80) |>
  mutate( predicted = format_numb(predicted), conf.low = format_numb(conf.low), conf.high = format_numb(conf.high),
          pred_CI = sprintf("%s (%s, %s)", predicted, conf.low, conf.high)) |>
  select(age, cancer_site, symptom=index_labels, gender, pred_CI) |>
  pivot_wider(id_cols = c(gender, cancer_site, symptom), 
              names_from = age,
              values_from = pred_CI, names_vary = "slowest", names_prefix = "Age ") |>
  arrange(gender, cancer_site, symptom) 

cancer_preds_40to80 |>
  flextable() |>
  set_header_labels(gender = "", symptom = "", cancer_site="") |>
  add_header_row(values = c("Sex", "Cancer site", "Symptom", "Probability of cancer % (95% CI)"), colwidths = c(1,1,1,3)) |>
  merge_v(j=c(1,2,3)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "TableA2.docx"))


# Table A.3 ----
load(file="./results/cancer_stage_cohort/allmodels_spline_pred.RData")
risk_advanced_preds <- allmodels_spline_pred; rm(allmodels_spline_pred)
risk_advanced_preds <- risk_advanced_preds |> rename(index_labels=sx)

group_cols <- c("gender", "cancer_site", "index_labels")
risk_40_advanced <- get_pred_at_age(risk_advanced_preds, 40, group_cols)
risk_60_advanced <- get_pred_at_age(risk_advanced_preds, 60, group_cols)
risk_80_advanced <- get_pred_at_age(risk_advanced_preds, 80, group_cols)

risk_advanced_40to80 <- rbind(risk_40_advanced, risk_60_advanced, risk_80_advanced) |>
  mutate( predicted = format_numb(predicted), conf.low = format_numb(conf.low), conf.high = format_numb(conf.high),
          pred_CI = sprintf("%s (%s, %s)", predicted, conf.low, conf.high)) |>
  select(age, cancer_site, symptom=index_labels, gender, pred_CI) |>
  pivot_wider(id_cols = c(gender, cancer_site, symptom), 
              names_from = age,
              values_from = pred_CI, names_vary = "slowest", names_prefix = "Age ") |>
  arrange(gender, cancer_site, symptom) 

risk_advanced_40to80 |>
  flextable() |>
  set_header_labels(gender = "", symptom = "", cancer_site="") |>
  add_header_row(values = c("Sex", "Cancer site", "Symptom", "Probability of advanced stage cancer in cases % (95% CI)"), colwidths = c(1,1,1,3)) |>
  merge_v(j=c(1,2,3)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "TableA3.docx"))

# Table A.4 ----
advanced_predict <- read.csv(file = "./results/bootstrap_10000/advanced_predict_bs.csv")
early_predict <- read.csv(file = "./results/bootstrap_10000/early_predict_bs.csv")

final_predict <- advanced_predict |> mutate(stage='advanced') |> bind_rows(early_predict |> mutate(stage='early')) |>
  rename(predicted=prob)

group_cols <- c("stage", "gender", "cancer_site", "index_labels")
final_predict_40 <- get_pred_at_age(final_predict, 40, group_cols)
final_predict_60 <- get_pred_at_age(final_predict, 60, group_cols)
final_predict_80 <- get_pred_at_age(final_predict, 80, group_cols)

final_predict_40to80 <- rbind(final_predict_40, final_predict_60, final_predict_80) |>
  mutate( predicted = format_numb(predicted), conf.low = format_numb(conf.low), conf.high = format_numb(conf.high),
          pred_CI = sprintf("%s (%s, %s)", predicted, conf.low, conf.high)) |>
  select(age, cancer_site, symptom=index_labels, gender, pred_CI) |>
  pivot_wider(id_cols = c(stage, gender, cancer_site, symptom), 
              names_from = age,
              values_from = pred_CI, names_vary = "slowest", names_prefix = "Age ") |>
  arrange(stage, gender, cancer_site, symptom) 

final_predict_40to80 |>
  flextable() |>
  set_header_labels(gender = "", symptom = "", cancer_site="") |>
  add_header_row(values = c("Stage", "Sex", "Cancer site", "Symptom", "Probability of (non-)advanced stage cancer % (95% CI)"), 
                 colwidths = c(1,1,1,1,3)) |>
  merge_v(j=c(1,2,3)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "TableA4.docx"))
