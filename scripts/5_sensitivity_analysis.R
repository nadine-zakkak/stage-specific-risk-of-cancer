# Purpose: sensitivity analysis - repeat stage (only) analysis for 
# (1) advanced stage (IV) vs non-advanced stage (I-III)
# (2) exclude screen detected cancers
# Author: Nadine Zakkak

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(stringr, tidyr, dplyr, marginaleffects, ggplot2, purrr, splines, broom)

source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_stage_dataprep.R") #import data from MySQL

if(!dir.exists("./results/sensitivity_analysis/")) {dir.create("./results/sensitivity_analysis/")}

# Predict function ----
model_predict <- function(model, pred_data){
  predictions(model,
              newdata = pred_data,
              vcov = "HC3" # robust standard error
  ) |>
    mutate(sx = factor(case_when(
      sx_eventtype_4  == 1 ~ "CIBH",
      sx_eventtype_14 == 1 ~ "Rectal Bleeding",
      sx_eventtype_14 == 0 & sx_eventtype_4 == 0 ~ "Neither"
    ), 
    levels = c("Rectal Bleeding", "CIBH",  "Neither")))
}

# 1) Advanced Stage (IV) vs Early Stage (I-III) ----
print("Advanced stage (IV)")
## Output directory ----
output_dir <- "./results/sensitivity_analysis/stage/"
if(!dir.exists(output_dir)) { dir.create(output_dir)}

## define outcomes ----
# stage binary cut-off (IV - advanced, I-III - early)
df <-
  df |>
  mutate(stage_bin = ifelse(stage_best == 4, "advanced",
                            ifelse(stage_best == -99, "missing",
                                   "early")))
df_cc <- df |> 
  filter(stage_bin != "missing") |> # complete case - excl. records with missing stage info  
  mutate(advanced = stage_bin == "advanced") |>
  mutate(sx_eventtype_4  = as.factor(sx_eventtype_4),
         sx_eventtype_14 = as.factor(sx_eventtype_14))


# observed proportions ----
observed_prop <-  
  df_cc |>
  mutate(age_decade = floor(age_decade)) |> # groups into 10y age category (30-40, 50-60,...)
  group_by(cr_eventtype, gender, sx, age_decade) |>
  summarise(n = sum(advanced),
            N = n()) |>
  mutate(prop = n / N ,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper) |>
  filter(sx != "Both") |>
  rowwise() |>
  mutate(cr_eventtype = names(cancer_codes[which(cancer_codes == cr_eventtype)]),
         gender       = ifelse(gender == 1, "Men", "Women"),
         age          = 60 + age_decade*10) |>
  ungroup() 

# Logistic regression with age as natural cubic spline with knot at 60yo and symptoms ------
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
  rowwise() |>
  mutate(cancer_site = names(cancer_codes[which(cancer_codes == cr_eventtype)]),
         gender = ifelse(gender == 1, "Men", "Women")) |>
  ungroup()

allmodels_spline_tidy <-
  allmodels_spline |>
  mutate(estimates = map(glms, tidy, conf.int = T, exponentiate = TRUE),
         stat = map(glms, glance)) |>
  select(cancer_site, gender, estimates, stat) |>
  unnest(cols = c(estimates, stat)) |>
  select(cancer_site, gender, term = term, OR = estimate, lower = conf.low, upper = conf.high, p = p.value, AIC)

write.csv(allmodels_spline_tidy, file = paste0(output_dir, "allmodels_spline.csv"), row.names = F)

#predictions
allmodels_spline_pred <- allmodels_spline |>
  mutate(
    predictions = map(.x = glms, 
                      ~ model_predict(., stage_predict_data))
  ) |>
  select(cancer_site, gender, predictions) |>
  unnest(predictions) |>
  mutate(age = 60 + age_decade*10)

write.csv(allmodels_spline_pred, file = paste0(output_dir, "allmodels_spline_pred.csv"), row.names = F)

# with observed
allmodels_spline_pred_plot_observ <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = sx), 
            size = .5, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = sx), 
              alpha = .2, data = allmodels_spline_pred) +
  geom_point(aes(x = age, y = 100*prop, color = sx),
             data = observed_prop) +
  geom_step(aes(x = age, y = 100*prop, color = sx),
            size = .5, data = observed_prop) +
  scale_fill_manual(values  = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  legend +
  facet_grid(gender ~ cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of stage IV cancer in cancer cases",
       y = y_label,
       x = x_label)
ggsave(plot = allmodels_spline_pred_plot_observ, filename = paste0(output_dir, "allmodels_spline_pred_plot_observ.png"), 
       unit = "cm", width = 20, height = 15)

allmodels_spline_pred_plot <-
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = sx), 
            size = .5, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = sx), 
              alpha = .2, data = allmodels_spline_pred) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme +
  labs(title = "Risk of stage IV cancer in cancer cases",
       y = y_label,
       x = x_label)
ggsave(plot = allmodels_spline_pred_plot, filename = paste0(output_dir, "allmodels_spline_pred_plot.png"), 
       unit = "cm", width = 20, height = 15)

rm(advanced_df, advanced_plot, early_df, early_plot, df, 
   list = ls(pattern = "^allmodels"),
   cancer_pred, stage_pred, df_cc, observed_prop, output_dir)


# 2) Exclude screen-detected cancers from stage cohort -----
print("Screen detected cancers")
 
# Load libraries 
if(!require('pacman'))install.packages('pacman')
pacman::p_load(stringr, tidyr, dplyr, marginaleffects, ggplot2, purrr, splines, broom)

source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_stage_dataprep.R")

## Output directory ----
output_dir <- "./results/sensitivity_analysis/excl_screen/"
if(!dir.exists(output_dir)) { dir.create(output_dir)}

## define outcomes ----
df <- df |> mutate(screen = final_route == "Screening")
df_cc <- df_cc |> mutate(screen = final_route == "Screening")

df_noscreen <- df_cc |> filter(!screen) # complete case - excl. records with missing stage info  

#What is the stage distribution by screening status? 
t <- df |>
  mutate(screen = ifelse(final_route == "Screening", 'Yes', 'No')) |>
  group_by(screen, stage_bin) |>
  summarise(n=n()) |>
  mutate(prop = 100*n/sum(n)) |>
  bind_rows(df |> group_by(stage_bin) |> summarise(n=n()) |> mutate(prop = 100*n/sum(n)) |> mutate(screen="Total")) |>
  mutate(screen    = factor(screen, levels = rev(c('Yes','No', 'Total'))),
         stage_bin = factor(stage_bin, levels = rev(c('early', 'advanced', 'missing')))) |>
  ggplot() +
  geom_bar(aes(y = prop, x = screen, fill = stage_bin), stat = "identity") +
  labs(y="Proportion (%)", x="Detected via screening") +
  guides(fill=guide_legend(title = "Stage")) +
  coord_flip()
ggsave(plot = t, paste0(output_dir, "screen_stage.png"),
       width = 30, height = 15, unit = "cm")

# observed proportions ----
observed_prop <-  
  df_noscreen |>
  mutate(age_decade = floor(age_decade)) |> # groups into 10y age category (30-40, 50-60,...)
  group_by(cr_eventtype, gender, sx, age_decade) |>
  summarise(n = sum(advanced),
            N = n()) |>
  mutate(prop = n / N ,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper
  ) |>
  filter(sx != "Both") |>
  rowwise() |>
  mutate(cancer_site = names(cancer_codes[which(cancer_codes == cr_eventtype)]),
         gender = ifelse(gender == 1, "Men", "Women")) |>
  ungroup() |>
  mutate(age = 60+age_decade*10) 

# Logistic regression with age as natural cubic spline with knot at 60yo and symptoms ------
allmodels_spline <-
  df_noscreen |>
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

# get model coefficients
allmodels_spline_tidy <-
  allmodels_spline |>
  mutate(estimates = map(glms, tidy, conf.int = T, exponentiate = TRUE),
         stat = map(glms, glance)) |>
  select(cancer_site, gender, estimates, stat) |>
  unnest(cols = c(estimates, stat)) |>
  select(cancer_site, gender, term = term, OR = estimate, lower = conf.low, upper = conf.high, p = p.value, AIC)

write.csv(allmodels_spline_tidy, file = paste0(output_dir, "allmodels_spline.csv"), row.names = F)


#predictions
allmodels_spline_pred <- allmodels_spline |>
  mutate(
    predictions = map(
      .x = glms, 
      ~ model_predict(.,stage_predict_data)
    )
  ) |>
  select(cancer_site, gender, predictions) |>
  unnest(predictions) |>
  mutate(age = 60 + age_decade*10)
write.csv(allmodels_spline_pred, file = paste0(output_dir, "allmodels_spline_pred.csv"), row.names = F)


# with observed
allmodels_spline_pred_plot_observ <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = sx), 
            size = linesize, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = sx), 
              alpha = .2, data = allmodels_spline_pred) +
  geom_point(aes(x = age, y = 100*prop, color = sx),
             data = observed_prop) +
  geom_step(aes(x = age, y = 100*prop, color = sx),
            size = .5, data = observed_prop) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  legend +
  facet_grid(gender~cancer_site)+
  theme_light() +
  theme + 
  labs(title = "Risk of advanced cancer in cancer cases",
       caption = "Excluding screen detected cancer",
       y = y_label,
       x = x_label)
ggsave(plot = allmodels_spline_pred_plot_observ, filename = paste0(output_dir, "allmodels_spline_pred_plot_observ.png"), 
       unit = "cm", width = 20, height = 15)

# without observed
allmodels_spline_pred_plot <- 
  ggplot() +
  geom_line(aes(x = age, y = 100*predicted, color = sx), 
            size = linesize, data = allmodels_spline_pred) +
  geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = sx), 
              alpha = .2, data = allmodels_spline_pred) +
  scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
  legend +
  facet_grid(gender ~ cancer_site)+
  theme_light() +
  theme + 
  labs(title = "Risk of advanced cancer in cancer cases",
       caption = "Excluding screen detected cancer",
       y = y_label,
       x = x_label)
ggsave(plot = allmodels_spline_pred_plot, filename = paste0(output_dir, "allmodels_spline_pred_plot.png"), 
       unit = "cm", width = 20, height = 15)
