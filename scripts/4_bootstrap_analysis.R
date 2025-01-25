# Purpose: Combine analysis: Find the risk of advanced and non-advanced cancer in patients with rectal bleeding or change in bowel habit
# using bootstrap to produce confidence interval
# Author: Nadine Zakkak

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(stringr, lubridate, forcats, tidyr, dplyr, marginaleffects,
               splines, purrr, ggplot2, flextable)

# global variables ----
source("./scripts/00_functions_bootstrap.R") #load bootstrap functions
source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
seed_script <- 153525
n_bootstrap <- 10000 #set number of bootstrap rep
output_dir <- paste0("./results/bootstrap_", n_bootstrap, "/")
if(!dir.exists(output_dir)) { dir.create(output_dir) }

set.seed(seed_script)
#cancer prediction -----
source("./scripts/00_cancer_indicidence_dataprep.R")

## prepare data for modelling ----
cancer_df <- 
  df |>
  mutate(age_decade = (age_sx - 60)/10) #centre ate at 60yo and present as decades
cancer_df <- cancer_df |> 
  mutate(colon  = !is.na(cancer_site) & cancer_site == cancer_codes[['Colon']],
         rectal = !is.na(cancer_site) & cancer_site == cancer_codes[['Rectal']])
rm(df); gc()

## bootstrap ----
#repeat bootstrap function n times
print("colon cancer")
colon_predict <- bind_rows(lapply(1:n_bootstrap, 
                                  function(i) bootstrap(cancer_df, 
                                                        as.formula(paste("colon~",paste(c("ns(age_decade, knots = c(-2, 0, 1, 2)) ", "index_eventtype"), collapse="+"))),
                                                        groups = syms("gender"),
                                                        stage = F, incidence_predict_data)), .id = "column_label") |>
  mutate(cancer_site = "Colon")

print("rectal cancer")
rectal_predict <- bind_rows(lapply(1:n_bootstrap, 
                                   function(i) bootstrap(cancer_df, 
                                                         as.formula(paste("rectal~",paste(c("ns(age_decade, knots = c(-2, 0, 1, 2)) ", "index_eventtype"), collapse="+"))),
                                                         groups = syms("gender"),
                                                         stage = F, incidence_predict_data)), .id = "column_label") |>
  mutate(cancer_site = "Rectal")

cancer_predict <- bind_rows(colon_predict, rectal_predict)

write.csv(cancer_predict, file = paste0(output_dir, "cancer_predict_bs.csv"), row.names = F)

# to check with other analysis 
cancer_plot <- cancer_predict |>
  group_by(age_decade, index_labels, gender, cancer_site) |>
  summarise(conf.low = quantile(predicted, 0.025),
            conf.high = quantile(predicted, 0.975),
            predicted = mean(predicted)) |>
  mutate(index_labels = factor(index_labels, c('Rectal Bleeding', 'CIBH')),
         cancer_site = factor(cancer_site, c('Colon', 'Rectal'))) |>
  ggplot() +
  geom_line(aes(x = 60+(age_decade*10), y = 100*predicted, color = index_labels), 
            size = linesize) +
  geom_ribbon(aes(x = 60+(age_decade*10), ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
              alpha = .2) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of Cancer",
       caption = "bootstrap",
       y = y_label,
       x = x_label)
ggsave(plot = cancer_plot, filename = paste0(output_dir, "risk_cancer_bs.png"),
       unit = "cm", width = 20, height = 15)

#stage prediction ------
source("./scripts/00_cancer_stage_dataprep.R")

## define outcomes & symptoms ----
stage_df <- df |> 
  filter(stage_bin != "missing") |> # complete case - excl. records with missing stage info  
  mutate(advanced = stage_bin == "advanced")
stage_df <-
  stage_df |>
  mutate(sx_eventtype_4 = factor(sx_eventtype_4, c(0,1)),
         sx_eventtype_14 = factor(sx_eventtype_14, c(0,1)))

#bootstrap ----
print("stage")
stage_predict <- bind_rows(lapply(1:n_bootstrap, function(i) bootstrap(stage_df, 
                                                                       as.formula(paste("advanced~",paste(c("ns(age_decade, knots = c(0))", "sx_eventtype_4", "sx_eventtype_14"), collapse="+"))),
                                                                       groups = syms(c("cr_eventtype","gender")),
                                                                       stage = T, stage_predict_data)), .id = "column_label") 

write.csv(stage_predict, file = paste0(output_dir, "stage_predict_bs.csv"), row.names = F)

# to check with other analysis 
stage_plot <- stage_predict |>
  group_by(age_decade, sx, gender, cancer_site) |>
  summarise(conf.low = quantile(predicted, 0.025),
            conf.high = quantile(predicted, 0.975),
            predicted = mean(predicted)) |>
  mutate(index_labels = factor(sx, c('Rectal Bleeding', 'CIBH', 'Neither')),
         cancer_site = factor(cancer_site, c('Colon', 'Rectal'))) |>
  ggplot() +
  geom_line(aes(x = 60+(age_decade*10), y = 100*predicted, color = index_labels), 
            size = linesize) +
  geom_ribbon(aes(x = 60+(age_decade*10), ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
              alpha = .2) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  ylim(1,3.5) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of advanced cancer in cancer cases",
       caption = "bootstrap",
       y = y_label,
       x = x_label)
ggsave(plot = stage_plot, filename = paste0(output_dir, "stage_cancer_bs.png"),
       unit = "cm", width = 20, height = 15)

# Combine analysis -----
print("combine analysis")
advanced_predict <-
  cancer_predict |> rename(cancer_predict = predicted) |>
  inner_join(stage_predict |> rename(stage_predict = predicted), 
             by = c("age_decade"="age_decade", "index_labels"="sx", "gender"="gender", "cancer_site"="cancer_site",
                    "column_label"="column_label")) |>
  mutate(prob = cancer_predict*stage_predict) |>
  group_by(age_decade, index_labels, gender, cancer_site) |>
  summarise(conf.low  = quantile(prob, 0.025),
            conf.high = quantile(prob, 0.975),
            prob = mean(prob)) |>
  mutate(index_labels = factor(index_labels, c('Rectal Bleeding', 'CIBH')),
         cancer_site  = factor(cancer_site, c('Colon', 'Rectal')),
         age           = 60 + (age_decade*10))

write.csv(advanced_predict, file = paste0(output_dir, "advanced_predict_bs.csv"), row.names = F)

early_predict <-
  cancer_predict |> rename(cancer_predict = predicted) |>
  inner_join(stage_predict |> rename(stage_predict = predicted), 
             by = c("age_decade"="age_decade", "index_labels"="sx", "gender"="gender", "cancer_site"="cancer_site",
                    "column_label"="column_label")) |>
  mutate(prob = cancer_predict*(1-stage_predict)) |>
  group_by(age_decade, index_labels, gender, cancer_site) |>
  summarise(conf.low  = quantile(prob, 0.025),
            conf.high = quantile(prob, 0.975),
            prob = mean(prob)) |>
  mutate(index_labels = factor(index_labels, c('Rectal Bleeding', 'CIBH')),
         cancer_site  = factor(cancer_site, c('Colon', 'Rectal')),
         age          = 60 + (age_decade*10))

write.csv(early_predict, file = paste0(output_dir, "early_predict_bs.csv"), row.names = F)


advanced_plot <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*prob, color = index_labels), 
            size = linesize, data = advanced_predict) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
              alpha = .2, data = advanced_predict) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  ylim(0,3.5) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of advanced stage cancer",
       y = y_label,
       x = x_label)
ggsave(plot = advanced_plot, filename = paste0(output_dir, "risk_advanced_cancer_bs.png"),
       unit = "cm", width = 20, height = 15)

early_plot <-
  ggplot() +
  geom_line(aes(x = age, y = 100*prob, color = index_labels), 
            size = linesize, data = early_predict) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
              alpha = .2, data = early_predict) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  ylim(0,3.5) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of non-advanced stage cancer",
       y = y_label,
       x = x_label)
ggsave(plot = early_plot, filename = paste0(output_dir, "risk_early_cancer_bs.png"),
       unit = "cm", width = 20, height = 15)



# Get info about age at highest risk by symptom, cancer and sex -----
# Advanced stage:
advanced_predict <- read.csv(file = paste0(output_dir, "advanced_predict_bs.csv"))
early_predict <- read.csv(file = paste0(output_dir, "early_predict_bs.csv"))

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
  merge_v(j = c(1,2,3)) |>
  hline(i = c(2,4,6,8,10,12,14)) |>
  fontsize(part = "all", size = 9) |>
  bold(part = "header") |>
  font(part = "all", fontname =  "Calibri") |>
  save_as_docx(path = paste0(output_dir, "max_abs_risk.docx"))
