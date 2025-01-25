## Purpose: descriptives of rectal bleeding & change in bowel habit cohorts
## Author:  Nadine Zakkak 

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(dplyr, tidyr, purrr, stringr, forcats, lubridate, scales, gt, flextable, ggplot2)

# Output directory ----
output_dir <- "./results/cancer_incidence_cohort/"
if(!dir.exists(output_dir)) { dir.create(output_dir) }

# Data prep ----
source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_incidence_dataprep.R") #import data from MySQL
df_crc <- df |> filter(crc == 1)

# % Stage Missingness by Year -----
missingness <- 
  as.data.frame(table(year(df_crc$first_crc_date), 
                      df_crc$stage_bin, useNA = "ifany")) |> 
  pivot_wider(id_cols = "Var1", names_from = "Var2", values_from = "Freq") |>
  rename(Year = "Var1") |> 
  mutate(total = rowSums(across(where(is.numeric))),
         `Missingness %` = round(100*missing/total, 1)) |> 
  select(-total) |> 
  select(Year, "non-advanced", advanced, missing, "Missingness %")

missingness |>
  mutate(`Missingness %` = paste0(`Missingness %`, "%")) |>
  flextable() |>
  add_header_row(values = c("", "Colorectal Cancer Stage", ""), colwidths = c(1,3,1)) |>
  add_table_theme() |>
  save_as_docx(path = paste0(output_dir, "missingness_year.docx"))

missingness_plot <- missingness |>
  ggplot() +
  geom_bar(aes(x = Year, y = `Missingness %`), stat = "identity")
ggsave(plot = missingness_plot, file = paste0(output_dir, "missingness_year.png"),
       width = 25, height = 13, unit = "cm")

completeness_plot <- missingness |>
  mutate(`Complete Stage %` = 100 - `Missingness %`) |>
  ggplot() +
  geom_bar(aes(x = Year, y = `Complete Stage %`), stat = "identity") +
  scale_y_continuous(breaks = seq(0, 100, by = 20))
ggsave(plot = completeness_plot, file = paste0(output_dir, "completestage_year.png"),
       width = 25, height = 13, unit = "cm")

# % stage missingness by year per cancer site -----
missingness_site_plot <- df_crc |>
  group_by(cancer_site, year(first_crc_date)) |>
  count(stage_bin) |>
  rename(Year = `year(first_crc_date)`) |>
  pivot_wider(id_cols = c("cancer_site", Year), names_from = stage_bin, values_from = n) |>
  mutate(total = rowSums(across(where(is.numeric))),
         `Missingness %` = round(100*missing/total, 1))  |>
  rowwise() |>
  mutate(site = names(cancer_codes)[cancer_codes == cancer_site]) |>
  ggplot() +
  geom_bar(aes(x = Year, y = `Missingness %`), stat = "identity") +
  facet_wrap(~site)

ggsave(plot = missingness_site_plot, file = paste0(output_dir, "missingness_year_site.png"),
       width = 25, height = 13, unit = "cm")


# Descriptive summaries ----
# variables to stratify by 
categ <- c(quo(gender), quo(age10_cat), quo(imd2015_5), quo(index_eventtype))
## % cancer site by covariate(s) ----
desc_summary <-
  lapply(categ, FUN = function(x){ # summary per covariate
    df |>
      mutate(!!x := as.character(!!x)) |>
      group_by(!!x, cancer_site) |>
      summarise(n = n()) |>
      mutate(N = sum(n),
             prop = n/N,
             lower = calculate_ci(n, N)$lower,
             upper = calculate_ci(n, N)$upper,
             summary = str_glue("{n} ({format_numb(100*prop)}%)")) |>
      mutate(!!x := as.factor(!!x)) |>
      ungroup()
  })

### Table ----
desc_summary_df <- bind_rows(desc_summary) |>
  bind_rows( # add info of total patients
    df |> group_by(cancer_site) |>
      summarise(n = n()) |>
      mutate(N = sum(n),
             prop = n/N,
             lower = calculate_ci(n, N)$lower,
             upper = calculate_ci(n, N)$upper,
             summary = str_glue("{n} ({format_numb(100*prop)}%)"),
             Total = "Total") |>
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
write.csv(desc_summary_df_wide, paste0(output_dir, "desc_all.csv"), row.names = F)

desc_summary_df_wide |>
  select(covariate, value, N, summary_Colon, summary_Rectal) |>
  rename_with(.fn = ~gsub("summary_", "", .x), cols=starts_with("summary_")) |>
  flextable() |>
  merge_v(j=1) |>
  add_table_theme() |>
  save_as_docx(path = paste0(output_dir, "desc_all.docx"))


### Dot Plot  ----
bind_rows(desc_summary) |>
  recode_gender() |>
  mutate(imd2015_5 = ifelse(imd2015_5 == 1, "1 - Least", ifelse(imd2015_5 == 5, "5 - Most", imd2015_5)),
         prop = 100*prop,
         lower = 100*lower, 
         upper = 100*upper) |>
  rowwise() |>
  mutate(index_eventtype = ifelse(!is.na(index_eventtype), names(symptom_codes)[symptom_codes == index_eventtype], NA), 
         cancer_site = ifelse(!is.na(cancer_site), names(cancer_codes)[cancer_codes == cancer_site], NA)) |>
  rename(Gender = gender,
         Age = age10_cat,
         `IMD (2015)` = imd2015_5,
         `Index Symptom` = index_eventtype) |>
  pivot_longer(cols = c(Age, `IMD (2015)`, `Index Symptom`, Gender), names_to = "covariate")|> 
  filter(!is.na(value), !is.na(cancer_site)) |>
  select(covariate, value, cancer_site, prop, lower, upper)|>
  ggplot() +
  geom_point(aes(value, prop)) +
  geom_errorbar(aes(x = value, y = prop, ymin = lower, ymax = upper), width = .1) +
  facet_grid(cancer_site~covariate, scales = "free_x") +
  labs(x = "",
       y = "Proportion of Patients (%)") +
  theme(axis.text.x = element_text(size = 8))
ggsave(paste0(output_dir, "desc_all.png"),
       width = 25, height = 13, unit = "cm", dpi = 150)

rm(list = ls())
gc()

