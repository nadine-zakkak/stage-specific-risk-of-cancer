## Purpose: Cancer incidence analysis, using rectal bleeding & change in bowel habit cohorts
## Author:  Nadine Zakkak 

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(forcats, tidyr, dplyr, purrr, marginaleffects, broom, splines, gt, ggplot2)

source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_incidence_dataprep.R") #extract data form MySQL
output_dir <- "./results/cancer_incidence_cohort/"
if(!dir.exists(output_dir)) { dir.create(output_dir) }

## prepare data for modelling ----
df_backup <- df
df <- 
  df |>
  mutate(age_decade = (age_sx - 60)/10) #centre ate at 60yo and present as decades

## define outcomes ----
df <- df |> 
  mutate(colon = !is.na(cancer_site) & cancer_site == cancer_codes[['Colon']],
         rectal = !is.na(cancer_site) & cancer_site == cancer_codes[['Rectal']])

# Predict function -----
model_predict <- function(model, pred_data){
  predictions(model,
              newdata = pred_data,
              vcov = "HC3") |> # robust standard error
    rowwise() |>
    mutate(index_labels = factor(names(symptom_codes)[symptom_codes == index_eventtype])) |>
    mutate(index_labels = fct_relevel(index_labels, 'Rectal Bleeding'))
}

# Plot function ----
plot_predict <- function(predict_df, observe_df, colours, legend_title, legend, theme, title, y_label, x_label){
  ggplot() +
    geom_line(aes(x = age, y = 100*predicted, color = index_labels),
              size = linesize, data = predict_df) +
    geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels),
                alpha = .2, data = predict_df) +
    geom_point(aes(x = age, y = 100*prop, color = index_labels),
               data = observe_df) +
    geom_step(aes(x = age, y = 100*prop, color = index_labels),
              size = .5, data = observe_df) +
    scale_fill_manual(values = colours, legend_title) +
    scale_color_manual(values = colours, legend_title) +
    legend +
    facet_grid(~gender)+
    theme_light() +
    theme +
    labs(title = title,
         y = y_label,
         x = x_label)
}

# Observed proportions ----
observed_prop <- df |> 
  mutate(age_decade = floor(age_decade)) |> 
  group_by(index_eventtype, gender, age_decade) |> 
  summarise(Colon = sum(colon),
            Rectal = sum(rectal), 
            N = n()) |>
  pivot_longer(cols = c(Colon, Rectal), names_to = "cancer_site", values_to = "n") |>
  mutate(prop = n / N,
         lower = calculate_ci(n,N)$lower,
         upper = calculate_ci(n,N)$upper) |>
  recode_gender() |>
  rowwise() |>
  mutate(index_labels = factor(names(symptom_codes)[symptom_codes == index_eventtype])) |>
  mutate(index_labels = fct_relevel(index_labels, 'Rectal Bleeding')) |>
  ungroup() |>
  arrange(cancer_site, index_eventtype) |>
  mutate(age = 60+age_decade*10)

# Colon ----

# Logistic regression adjusted for age as natural cubic spline with knot at 40, 60,70 and 80yo and symptoms
colon_spline <-
  df |>
  group_by(gender) |>
  nest() |>
  mutate(glms = map(
    data, ~glm(colon ~ ns(age_decade, knots = c(-2, 0, 1, 2)) + index_eventtype,
               family = binomial(link = "logit"),
               data = .x)
  )) |>
  select(gender, glms) |>
  recode_gender() |>
  ungroup()

colon_spline_tidy <-
  colon_spline |>
  mutate(estimates = map(glms, tidy, conf.int = T, exponentiate = TRUE),
         stat = map(glms, glance)) |>
  select(gender, estimates, stat) |>
  unnest(cols = c(estimates, stat)) |>
  select(gender, term = term, OR = estimate, lower = conf.low, upper = conf.high, p = p.value, AIC)

write.csv(colon_spline_tidy, file = paste0(output_dir, "colon_spline.csv"), row.names = F)

colon_spline_pred <-
  colon_spline |>
  mutate(
    predictions = map(
      .x = glms, 
      ~ model_predict(.,incidence_predict_data) 
    )
  ) |>
  select(gender, predictions) |>
  unnest(predictions) |>
  mutate(age = 60+age_decade*10)

write.csv(colon_spline_pred, file = paste0(output_dir, "colon_spline_pred.csv"), row.names = F)

colon_spline_plot <- plot_predict(colon_spline_pred, observed_prop |> filter(cancer_site == "Colon"),
                                  colours = c(colours$`Rectal Bleeding`, colours$CIBH),
                                  legend_title = legend_title,legend = legend,
                                  theme = theme,
                                  x_label= x_label, y_label =  y_label, "Risk of Colon Cancer")
ggsave(plot = colon_spline_plot, filename = paste0(output_dir, "colon_spline_plot.png"), 
       unit = "cm", width = 20, height = 15)

# Rectum -----
# Logistic regression adjusted for age as natural cubic spline with knot at 40, 60, 70, 80yo and symptoms
rectal_spline <-
  df |>
  group_by(gender) |>
  nest() |>
  mutate(glms = map(
    data, ~glm(rectal ~ ns(age_decade, knots = c(-2, 0, 1, 2)) + index_eventtype,
               family = binomial(link = "logit"),
               data = .x)
  )) |>
  select(gender, glms) |>
  recode_gender() |>
  ungroup()

rectal_spline_tidy <-
  rectal_spline |>
  mutate(estimates = map(glms, tidy, conf.int = T, exponentiate = TRUE),
         stat = map(glms, glance)) |>
  select(gender, estimates, stat) |>
  unnest(cols = c(estimates, stat)) |>
  select(gender, term = term, OR = estimate, lower = conf.low, upper = conf.high, p = p.value, AIC)

write.csv(rectal_spline_tidy, file = paste0(output_dir, "rectal_spline.csv"), row.names = F)

rectal_spline_pred <-
  rectal_spline |>
  mutate(
    predictions = map(
      .x = glms, 
      ~ model_predict(.,incidence_predict_data) 
    )
  ) |>
  select(gender, predictions) |>
  unnest(predictions) |>
  mutate(age = 60+age_decade*10)

write.csv(rectal_spline_pred, file = paste0(output_dir, "rectal_spline_pred.csv"), row.names = F)

rectal_spline_plot <- plot_predict(rectal_spline_pred, observed_prop |> filter(cancer_site == "Rectal"), 
                                   colours = c(colours$`Rectal Bleeding`, colours$CIBH),
                                   legend_title = legend_title,legend = legend,
                                   theme = theme,
                                   x_label= x_label, y_label =  y_label, "Risk of Rectal Cancer")
ggsave(plot = rectal_spline_plot, filename = paste0(output_dir, "rectal_spline_plot.png"), 
       unit = "cm", width = 20, height = 15)

# Combined results ----
allspine_tidy <- bind_rows(colon_spline_tidy |> mutate(cancer_site = "Colon"),
                           rectal_spline_tidy |> mutate(cancer_site = "Rectal"))
allspine_tidy |>
  filter(grepl("index_eventtype", term)) |>
  pivot_wider(id_cols = c(gender, term), 
              names_from = cancer_site,
              values_from = c(OR, lower, upper, p)) |>
  select(contains("Colon"), contains("Rectal"), everything()) |>
  mutate(term = ifelse(term == "index_eventtype4", "Change in Bowel Habit", "Rectal Bleeding")) |>
  arrange(gender, term) |>
  group_by(gender, term) |>
  rename(symptom = term) |>
  gt() |>
  cols_merge(columns = c(OR_Colon, lower_Colon, upper_Colon), pattern="{1} ({2}, {3})") |>
  cols_merge(columns = c(OR_Rectal, lower_Rectal, upper_Rectal), pattern="{1} ({2}, {3})") |>
  fmt_number(columns = c(1,2,3,4,5,6,7,8), decimals = 2) |>
  gtsave(paste0(output_dir, "table_OR.rtf"))

allspline_pred <- bind_rows(colon_spline_pred |> mutate(cancer_site = "Colon"), rectal_spline_pred |> mutate(cancer_site = "Rectal"))
write.csv(allspline_pred, file = paste0(output_dir, "allspline_pred.csv"), row.names = F)
saveRDS(allspline_pred, file = paste0(output_dir, "allspline_pred.rds"))
save(allspline_pred, file=paste0(output_dir, "allspline_pred.RData"))


# with observed
allspline_pred_plot_observ <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = index_labels), 
            size = linesize, data = allspline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
              alpha = .2, data = allspline_pred) +
  geom_point(aes(x = age, y = 100*prop, color = index_labels),
             data = observed_prop ) +
  geom_step(aes(x = age, y = 100*prop, color = index_labels),
            size = .5, data = observed_prop) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  legend + 
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme + 
  labs(title = "Risk of cancer",
       y = y_label,
       x = x_label)
ggsave(plot = allspline_pred_plot_observ, filename = paste0(output_dir, "allspline_pred_plot_observ.png"), 
       unit = "cm", width = 20, height = 15)

# without observed proportions
allspline_pred_plot <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = index_labels), 
            size = linesize, data = allspline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
              alpha = .2, data = allspline_pred) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH), legend_title) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of cancer",
       y = y_label,
       x = x_label)
ggsave(plot = allspline_pred_plot, filename = paste0(output_dir, "allspline_pred_plot.png"), 
       unit = "cm", width = 20, height = 15)

rm(list = ls())
gc()

