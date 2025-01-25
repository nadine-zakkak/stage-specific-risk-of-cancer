# Purpose: Load and reshapee cancer registry data for analysis
# Author: Nadine Zakkak and Matthew Barclay

library(broom)
library(tidyr)
library(purrr)
library(splines)
library(marginaleffects)
library(boot)

source("./scripts/9a_load_cr_data.R")

# Modelling using Matt's method
nested_crc_data_mod <- 
  crc_data_mod |>
  group_by(cancer, male) |>
  nest() 

predictions_modelvars_only <- crc_data_mod |> 
  select(age_c) |> 
  unique() |> 
  arrange(age_c) |>
  select(age_c)

# Advanced 3 & 4
model_fit <- nested_crc_data_mod |>
  mutate(
    glms = map(
      data,  ~glm(cbind(stage34, stage12) ~ ns(age_c, knots = c(0)), family = binomial(link = "logit"), data = .x)
    )
  )

model_preds <- model_fit |>
  mutate(
    preds = map(
      .x = glms, 
      ~ predictions(., predictions_modelvars_only, vcov = TRUE, type = "link")
    )
  ) |>
  select(male, cancer, preds) |>
  unnest(cols = c(preds))

model_preds <-
  model_preds |>
  mutate(
    pr_lb = inv.logit(conf.low),
    pr_ub = inv.logit(conf.high),
    pr = inv.logit(predicted)) |>
  mutate(age = age_c*10+60,
         pr = 100*pr,
         pr_lb = 100*pr_lb,
         pr_ub = 100*pr_ub) |>
  select(rowid, male, cancer, age, pr, pr_lb, pr_ub, p.value)

model_preds

saveRDS(model_preds, "cr_model_preds.rds")
