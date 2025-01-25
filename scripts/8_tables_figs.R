# Results: Output all required figures and tables and some results for manuscript 
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


# Table 1 ----
# Descriptive of cancer incidence cohort
# data from 2a_cancer_incidence_descriptive.R
table1 <- read.csv("./results/cancer_incidence_cohort/desc_all.csv")
table1 |>
  select(covariate, value, N, summary_Colon, summary_Rectal) |>
  rename_with(.fn = ~gsub("summary_", "", .x), cols=starts_with("summary_")) |>
  flextable() |>
  set_header_labels(
    covariate = "", value = "", 
    Colon = "Colon cancer within 12 months",
    Rectal = "Rectal cancer within 12 months"
  ) |>
  merge_v(j=1) |>
  add_table_theme() |>
  save_as_docx(path=paste0(path_save, "Table1.docx"))

# Figure 2 ----
# Risk of cancer predictions 
# Data from 3a_cancer_incidence_analysis.R
cancer_preds <- readRDS(file = "./results/cancer_incidence_cohort/allspline_pred.rds") |>
  mutate(index_labels = factor(index_labels, levels = c("Rectal Bleeding", "CIBH")))

fig2 <- plot_preds_data(cancer_preds, "Probability of cancer (%)")
do.call(ggsave, c(list(plot = fig2, filename = paste0(path_save, "Fig2.jpeg")), save_image_var_one_and_half))


# Table 2 -----
# Odds ratios for cancer risk prediction models
# Data from 3a_cancer_incidence_analysis.R
colon_OR <- read.csv(file = "./results/cancer_incidence_cohort/colon_spline.csv")
rectal_OR <- read.csv(file = "./results/cancer_incidence_cohort/rectal_spline.csv")

combine_OR <- bind_rows(colon_OR |> mutate(cancer_site = "Colon"),
                        rectal_OR |> mutate(cancer_site = "Rectal"))

table2 <- combine_OR |>
  filter(grepl("index_eventtype", term)) |> #extract symptom coefficients only
  mutate(symptom = ifelse(term == "index_eventtype4", "Change in Bowel Habit", "Rectal Bleeding"),
         p = ifelse(p < 0.01, "<0.01", format_numb(p)),
         OR = format_numb(OR), lower = format_numb(lower), upper = format_numb(upper),
         OR_CI = paste0(OR, " (", lower, ", ", upper, ")")) |>
  bind_rows(expand.grid(gender = c("Women", "Men"), 
                        symptom = "Rectal Bleeding", 
                        cancer_site = c("Colon", "Rectal"),
                        OR_CI = "Reference",
                        p = "-")) |>
  pivot_wider(id_cols = c(gender, symptom), 
              names_from = cancer_site,
              values_from = c(OR_CI, p), names_vary = "slowest") |>
  arrange(gender, desc(symptom))

table2 |>
  flextable() |>
  set_header_labels(
    gender = "", symptom = "",
    OR_CI_Colon = "OR (95% CI)", p_Colon = "p-value",
    OR_CI_Rectal = "OR (95% CI)", p_Rectal = "p-value"
  ) |>
  add_header_row(values = c("", "Symptom", "Colon", "Rectal"), colwidths = c(1,1,2,2)) |>
  merge_v(j=1) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "Table2.docx"))

# Figure 3 -----
# Descriptive of cancer stage cohort (complete case)
# Data from 2b_stage_descriptive.R 
stage_desc <- read.csv("./results/cancer_stage_cohort/raw_desc_stage_cc.csv")

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


fig3 <-  ggplot() +
  geom_bar(aes(x = value, y = prop, fill = stage_best), stat = "identity", data = stage_desc) +
  ylim(-13, 101)+
  geom_text(size = 2, hjust = 0, colour = "#5A5A5A", aes(x = value, y = -13, label = label), data = total_label) +
  scale_fill_manual(values = rev(c("#5340ff", "lightblue", "#C47B7B", "darkred"))) +
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
do.call(ggsave, c(list(plot = fig3, filename = paste0(path_save, "Fig3.jpeg")), save_image_var_full))

# Figure 4 ----
# Risk of advanced stage in cancer cases predictions 
# Data from 3b_stage_analysis.R
adv_preds <- read.csv(file = "./results/cancer_stage_cohort/allmodels_spline_pred.csv") |> 
  rename(index_labels=sx) |>
  mutate(index_labels = factor(index_labels, levels = c("Rectal Bleeding", "CIBH", "Neither")))

fig4 <- plot_preds_data(adv_preds, "Probability of advanced stage cancer in cases (%)")
do.call(ggsave, c(list(plot = fig4, filename = paste0(path_save, "Fig4.jpeg")), save_image_var_one_and_half))

# Table 3 -----
# Odds ratios for advanced stage in cancer cases models
# Data from 3b_stage_analysis.R
adv_OR <- read.csv(file = "./results/cancer_stage_cohort/allmodels_spline.csv")

table3 <- 
  adv_OR |>
  filter(grepl("sx", term)) |>
  mutate(symptom = ifelse(term == "sx_eventtype_41", "Change in Bowel Habit", "Rectal Bleeding"),
         symptom = factor(symptom, levels = c("Rectal Bleeding", "Change in Bowel Habit"))) |>
  mutate(p = ifelse(p<0.01, "<0.01", format_numb(p)),
         OR = format_numb(OR), lower = format_numb(lower), upper = format_numb(upper),
         OR_CI = paste0(OR, " (", lower, ", ", upper, ")")) |>
  pivot_wider(id_cols = c(gender, symptom), 
              names_from = cancer_site,
              values_from = c(OR_CI, p), names_vary = 'slowest') |>
  arrange(gender, symptom)

table3 |>
  flextable() |>
  set_header_labels(
    gender = "", symptom = "",
    OR_CI_Colon = "OR (95% CI)", p_Colon = "p-value",
    OR_CI_Rectal = "OR (95% CI)", p_Rectal = "p-value"
  ) |>
  add_header_row(values = c("", "Symptom", "Colon", "Rectal"), colwidths = c(1,1,2,2)) |>
  merge_v(j=1) |>
  add_table_theme() |>
  save_as_docx(path = paste0(path_save, "Table3.docx"))

# Figure 5 -----
# Risk of advanced stage cancer predictions
# Data from 4_bootstrap_analysis.R
abs_adv_preds <- read.csv("./results/bootstrap_10000/advanced_predict_bs.csv") |> rename(predicted=prob) |>
  mutate(index_labels = factor(index_labels, levels = c("Rectal Bleeding", "CIBH")))

fig5 <- plot_preds_data(abs_adv_preds, "Probability of advanced stage cancer (%)") + ylim(0,3.5)
do.call(ggsave, c(list(plot = fig5, filename = paste0(path_save, "Fig5.jpeg")), save_image_var_one_and_half))

# Figure 6 -----
# Risk of non-advanced stage cancer predictions
# Data from 4_bootstrap_analysis.R
abs_non_adv_preds <- read.csv("./results/bootstrap_10000/early_predict_bs.csv") |> rename(predicted=prob) |>
  mutate(index_labels = factor(index_labels, levels = c("Rectal Bleeding", "CIBH")))

fig6 <- plot_preds_data(abs_non_adv_preds, "Probability of non-advanced stage cancer (%)") + ylim(0,3.5)
do.call(ggsave, c(list(plot = fig6, filename = paste0(path_save, "Fig6.jpeg")), save_image_var_one_and_half))


