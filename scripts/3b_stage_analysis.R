## Purpose: Cancer stage analysis, using colon and rectal cancer cohorts
## Author:  Nadine Zakkak 

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(forcats, tidyr, dplyr, purrr, marginaleffects, broom, splines, gt, ggplot2)

source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_stage_dataprep.R") #import data from MySQL

## Output directory ----
output_dir <- "./results/cancer_stage_cohort/"
if(!dir.exists(output_dir)) { dir.create(output_dir) }

df_cc <-
  df_cc |>
  mutate(sx_eventtype_4 = factor(sx_eventtype_4, c(0,1)),
         sx_eventtype_14 = factor(sx_eventtype_14, c(0,1)))

# Predict function ----
model_predict <- function(model, pred_data){
  predictions(model,
              pred_data,
              vcov = "HC3" # robust standard error
  ) |>
    mutate(sx = factor(case_when(
      sx_eventtype_4 == 1 ~ "CIBH",
      sx_eventtype_14 == 1 ~ "Rectal Bleeding",
      sx_eventtype_14 == 0 ~ "Neither"
    ), 
    levels = c("Rectal Bleeding", "CIBH",  "Neither")))
}

# observed proportions ----
observed_prop <-  
  df_cc |>
  mutate(age_decade = floor(age_decade)) |> # groups into 10y age category (30-40, 50-60,...)
  group_by(cr_eventtype, gender, sx, age_decade) |>
  summarise(n = sum(advanced),
            N = n()) |>
  mutate(prop = n / N ,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper
  ) |>
  filter(sx != "Both") |>
  recode_gender() |>
  rowwise() |>
  mutate(cr_eventtype = names(cancer_codes[which(cancer_codes == cr_eventtype)])) |>
  ungroup() |>
  mutate(age = 60+age_decade*10)

# Logistic regression with age as natural cubic spline with knot at 60yo and symptoms
print("Logistic regression with spline")
allmodels_spline <-
  df_cc |>
  group_by(cr_eventtype, gender) |>
  nest() |>
  mutate(glms = map(
    data, ~glm(advanced ~ ns(age_decade, knots = c(0)) + sx_eventtype_4 + sx_eventtype_14,
               family = binomial(link = "logit"),
               data = .x)
  )) |>
  select(cr_eventtype, gender, glms) |>
  recode_gender() |>
  rowwise() |>
  mutate(cancer_site = names(cancer_codes[which(cancer_codes == cr_eventtype)])) |>
  ungroup()

allmodels_spline_tidy <-
  allmodels_spline |>
  mutate(estimates = map(glms, tidy, conf.int = T, exponentiate = TRUE),
         stat = map(glms, glance)) |>
  select(cancer_site, gender, estimates, stat) |>
  unnest(cols = c(estimates, stat)) |>
  select(cancer_site, gender, term = term, OR = estimate, lower = conf.low, upper = conf.high, p = p.value, AIC)

write.csv(allmodels_spline_tidy, file = paste0(output_dir, "allmodels_spline.csv"), row.names = F)

# OR table
allmodels_spline_tidy |>
  filter(grepl("sx", term)) |>
  pivot_wider(id_cols     = c(gender, term), 
              names_from  = cancer_site,
              values_from = c(OR, lower, upper, p)) |>
  select(contains("Colon"), contains("Rectal"), everything()) |>
  mutate(term = ifelse(term == "sx_eventtype_41", "Change in Bowel Habit", "Rectal Bleeding"),
         term = factor(term, levels = c("Rectal Bleeding", "Change in Bowel Habit"))) |>
  arrange(gender, term) |>
  group_by(gender, term) |>
  rename(symptom = term) |>
  mutate(Colon    = paste0(format_numb(OR_Colon), " (", format_numb(lower_Colon), ", ", format_numb(upper_Colon), ")"),
         Rectal   = paste0(format_numb(OR_Rectal), " (", format_numb(lower_Rectal), ", ", format_numb(upper_Rectal), ")"),
         p_Colon  = format_numb(p_Colon),
         p_Rectal = format_numb(p_Rectal)) |>
  select(-c(OR_Colon, lower_Colon, upper_Colon, OR_Rectal, lower_Rectal, upper_Rectal)) |>
  select(gender, symptom, Colon, p_Colon, Rectal, p_Rectal) |>
  flextable() |>
  add_table_theme() |>
  save_as_docx(path=paste0(output_dir, "table_OR.docx"))

allmodels_spline_pred <-
  allmodels_spline |>
  mutate(
    predictions = map(
      .x = glms, 
      ~ model_predict(., stage_predict_data) 
    )
  ) |>
  select(cancer_site, gender, predictions) |>
  unnest(predictions) |>
  mutate(age = 60+age_decade*10)

write.csv(allmodels_spline_pred, file = paste0(output_dir, "allmodels_spline_pred.csv"), row.names = F)
save(allmodels_spline_pred, file=paste0(output_dir, "allmodels_spline_pred.RData"))

# with observed
allmodels_spline_plot_observ <-
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = sx), 
            size = linesize, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = sx), 
              alpha = .2, data = allmodels_spline_pred) +
  geom_point(aes(x = age, y = 100*prop, color = sx),
             data = observed_prop) +
  geom_step(aes(x = age, y = 100*prop, color = sx),
            size = .5, data = observed_prop) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of advanced cancer in cancer cases",
       y = y_label,
       x = x_label)
ggsave(plot = allmodels_spline_plot_observ, filename = paste0(output_dir, "allmodels_spline_plot_observ.png"), 
       unit = "cm", width = 20, height = 15)

# without observed
allmodels_spline_plot <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = sx), 
            size = linesize, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = sx), 
              alpha = .2, data = allmodels_spline_pred) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of advanced cancer in cancer cases",
       y = y_label,
       x = x_label)
ggsave(plot = allmodels_spline_plot, filename = paste0(output_dir, "allmodels_spline_plot.png"), 
       unit = "cm", width = 20, height = 15)


rm(list = ls())
gc()
