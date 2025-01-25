## Purpose: descriptives of colon & rectal cancer cohorts
## Author:  Nadine Zakkak 

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(dplyr, tidyr, stringr,  lubridate, gt, flextable, ggplot2)

source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_stage_dataprep.R") #import data from MySQL


# Output directory ----
output_dir <- "./results/cancer_stage_cohort/"
if(!dir.exists(output_dir)) { dir.create(output_dir) }

# % Stage Missingness by Year -----
missingness <- 
  as.data.frame(table(year(df$cr_eventdate), 
                      df$stage_bin, useNA = "ifany")) |> 
  pivot_wider(id_cols = "Var1", names_from = "Var2", values_from = "Freq") |>
  rename(Year = "Var1") |> 
  mutate(total = rowSums(across(where(is.numeric))),
         `Missingness %` = round(100*missing/total, 1)) |> 
  select(-total) |> 
  select(Year, "early", advanced, missing, "Missingness %")

missingness |>
  mutate(`Missingness %` = paste0(`Missingness %`, "%")) |>
  flextable() |>
  add_header_row(values = c("", "Colorectal Cancer Stage", ""), colwidths = c(1,3,1)) |>
  add_table_theme()|>
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
completeness_site_plot <- df |>
  group_by(cr_eventtype, year(cr_eventdate)) |>
  count(stage_bin) |>
  rename(Year = `year(cr_eventdate)`) |>
  pivot_wider(id_cols = c(cr_eventtype, Year), names_from = stage_bin, values_from = n) |>
  mutate(total = rowSums(across(where(is.numeric))),
         `Missingness %` = round(100*missing/total, 1),
         `Stage completeness %` = round(100-(100*missing/total), 1))  |>
  rowwise() |>
  mutate(site = names(cancer_codes)[cancer_codes == cr_eventtype]) |>
  ggplot() +
  geom_bar(aes(x = Year, y = `Stage completeness %`), stat = "identity") +
  facet_wrap(~site)

ggsave(plot = completeness_site_plot, file = paste0(output_dir, "completestage_year_site.png"),
       width = 25, height = 13, unit = "cm")

# Stage of cancer in each cancer site by covariate ----
## "Both" symptoms counted as such
df_long <-
  df |>
  recode_gender() |>
  mutate(
    age10_cat = ifelse(age10_cat %in% c("30-39", "40-49"), "<50", 
                       ifelse(age10_cat %in% c("80-89", "90+"), "80+", 
                              as.character(age10_cat))),
    imd2015_5 = ifelse(imd2015_5 == 1, "1 - Least", ifelse(imd2015_5 == 5, "5 - Most", imd2015_5))) |>
  rename(Gender = gender,
         Age = age10_cat,
         `IMD (2015)` = imd2015_5) |>
  pivot_longer(cols = c(Gender, Age, `IMD (2015)`), names_to = "covariate")

symptom_summary <-
  df |>
  select(gender, cr_eventtype, sx_eventtype_14, sx_eventtype_4, stage_best) |>
  mutate(Neither = sx_eventtype_14 == 0 & sx_eventtype_4 == 0,
         `Rectal Bleeding` = as.logical(as.numeric(sx_eventtype_14) - 1),
         CIBH = as.logical(as.numeric(sx_eventtype_4) -1)) |>
  pivot_longer(cols = c(`Rectal Bleeding`, CIBH, Neither), names_to = "symptom") |>
  group_by(cr_eventtype, symptom,stage_best) |>
  summarise(n = sum(value)) |>
  mutate(N = sum(n),
         prop = n/N,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper) |>
  rename(value = symptom) |>
  mutate(covariate = "Symptom") |>
  ungroup()

summary <-
  df_long |>
  group_by(cr_eventtype, covariate, value, stage_best) |>
  summarise(n = n()) |>
  mutate(N = sum(n),
         prop = n/N,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper) |>
  bind_rows(symptom_summary) |>
  bind_rows(df |> select(cr_eventtype, stage_best) |> 
              group_by(cr_eventtype, stage_best) |> # add total info
              summarise(n = n()) |>
              mutate(N = sum(n),
                     prop = n/N,
                     lower = calculate_ci(n, N)$lower,
                     upper = calculate_ci(n, N)$upper,
                     covariate = "Total",
                     value = "Total")) |>
  ungroup() |>
  mutate(prop = 100*prop,
         lower = 100*lower,
         upper= 100*upper) |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup()
write.csv(summary, paste0(output_dir, "raw_desc_stage_all.csv"), row.names = F)


summary_wide <-
  summary |> 
  mutate(prop = format_numb(prop),
         lower = format_numb(lower),
         upper = format_numb(upper),
         summary = str_glue("{n} ({prop}%)")) |>
  pivot_wider(id_cols =c(cr_eventtype, covariate, value, N),
              names_from = stage_best,
              values_from = summary, names_prefix = "summary_") |>
  mutate(covariate = factor(covariate, c("Gender", "Age", "IMD (2015)", "Symptom", "Total"))) |>
  arrange(covariate)

summary_wide |>
  rename_with(.fn = ~gsub("summary_", "", .x), cols=starts_with("summary_")) |>
  group_by(cancer_site, covariate) |>
  arrange(cancer_site) |>
  select(`Cancer site`=cancer_site, covariate, value, N, `1`, `2`, `3`, `4`, Missing=`-99`) |>
  flextable() |>
  add_header_row(values = c("", "Stage"), colwidths = c(4,5)) |>
  merge_v(j=c(1,2)) |>
  add_table_theme() |>
  save_as_docx(path=paste0(output_dir, "desc_stage_all.docx"))

total_label <-
  summary |>
  group_by(cancer_site, covariate, value) |>
  summarise(N = sum(n)) |>
  mutate(label = str_glue("N={N}")) |>
  ungroup()

p <- summary |>
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
  select(covariate, value, cancer_site, prop, lower, upper, stage_best, N) |>
  ggplot() +
  geom_bar(aes(x = value, y = prop, fill = stage_best), stat = "identity") +
  ylim(-12, 101)+
  geom_text(size = 3, hjust = 0, colour = "#5A5A5A", aes(x = value, y = -12, label = label), data = total_label) +
  scale_fill_manual(values = rev(c("#5340ff", "lightblue", "#C47B7B", "darkred", "grey"))) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  facet_grid(covariate~cancer_site, scales = "free_y") +
  labs(x = "",
       y = "Proportion of Patients (%)") +
  guides(fill=guide_legend(reverse=T, title="Stage")) +
  theme_minimal()+
  theme(strip.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 9, colour =  "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(plot = p, paste0(output_dir, "desc_all.png"),
       width = 30, height = 15, unit = "cm")

# Proportion of patients with both symptoms -----
prop_symptom <- prop.table(table(df$sx))*100
write.csv(prop_symptom, paste0(output_dir, "prop_symptom.csv"), row.names = F)

df |> group_by(random_sample, cohort, sx) |> summarise(n=n()) |> mutate(N=sum(n))

n_stage <- df |> group_by(stage_best) |> summarise(n=n()) |> mutate(N=sum(n), prop = format_numb(100*n/N))
write.csv(n_stage, paste0(output_dir, "n_stage.csv"), row.names = F)

n_stage_age <- 
  df |>
  mutate(age70_bin = case_when(
    age10_cat %in% c("30-39", "40-49", "50-59", "60-69") ~ "<70",
    age10_cat %in% c("70-79", "80-89", "90+") ~ ">=70")) |>
  group_by(age70_bin, stage_best) |> summarise(n=n()) |> mutate(N = sum(n), prop = format_numb(100*n/N))
write.csv(n_stage_age, paste0(output_dir, "n_stage_by_age.csv"), row.names = F)


df_cc <- df |> filter(stage_best != -99)

prop.table(table(df_cc$sx))*100

df_cc |> group_by(random_sample, cohort, sx) |> summarise(n=n()) |> mutate(N=sum(n))



# Descriptive Without missing stages (Complete case analysis) ----
df_cc <- df |> filter(stage_best != -99)

df_cc_long <-
  df_cc |>
  recode_gender() |>
  mutate(
    age10_cat = ifelse(age10_cat %in% c("30-39", "40-49"), "<50", 
                       ifelse(age10_cat %in% c("80-89", "90+"), "80+", 
                              as.character(age10_cat))),
    imd2015_5 = ifelse(imd2015_5 == 1, "1 - Least", ifelse(imd2015_5 == 5, "5 - Most", imd2015_5))) |>
  rename(Gender = gender,
         Age = age10_cat,
         `IMD (2015)` = imd2015_5) |>
  pivot_longer(cols = c(Gender, Age, `IMD (2015)`), names_to = "covariate")

symptom_summary_cc <-
  df_cc |>
  mutate(Neither = sx_eventtype_14 == 0 & sx_eventtype_4 == 0,
         `Rectal Bleeding` = as.logical(as.numeric(sx_eventtype_14) - 1),
         CIBH = as.logical(as.numeric(sx_eventtype_4) -1)) |>
  pivot_longer(cols = c(`Rectal Bleeding`, CIBH, Neither), names_to = "symptom") |>
  group_by(cr_eventtype, symptom, stage_best) |>
  summarise(n = sum(value)) |>
  mutate(N = sum(n),
         prop = n/N,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper) |>
  rename(value = symptom) |>
  mutate(covariate = "Symptom") |>
  ungroup()

summary_cc <-
  df_cc_long |>
  group_by(cr_eventtype, covariate, value, stage_best) |>
  summarise(n = n()) |>
  mutate(N = sum(n),
         prop = n/N,
         lower = calculate_ci(n, N)$lower,
         upper = calculate_ci(n, N)$upper) |>
  bind_rows(symptom_summary_cc) |>
  bind_rows(df_cc |>select(cr_eventtype, stage_best) |>
              group_by(cr_eventtype, stage_best) |> # add total info
              summarise(n = n()) |>
              mutate(N = sum(n),
                     prop = n/N,
                     lower = calculate_ci(n, N)$lower,
                     upper = calculate_ci(n, N)$upper,
                     covariate = "Total",
                     value = "Total")) |>
  ungroup() |>
  mutate(prop = 100*prop,
         lower = 100*lower,
         upper = 100*upper) |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup()
write.csv(summary_cc, paste0(output_dir, "raw_desc_stage_cc.csv"), row.names = F)

summary_cc_wide <-
  summary_cc |> 
  mutate(prop = format_numb(prop),
         lower = format_numb(lower),
         upper = format_numb(upper),
         summary = str_glue("{n} ({prop}%)")) |>
  pivot_wider(id_cols =c(cancer_site, covariate, value, N),
              names_from = stage_best,
              values_from = summary, names_prefix = "summary_") |>
  mutate(covariate = factor(covariate, c("Gender", "Age", "IMD (2015)", "Symptom", "Total"))) |>
  arrange(covariate)

summary_cc_wide |>
  rename_with(.fn = ~gsub("summary_", "", .x), cols=starts_with("summary_")) |>
  group_by(cancer_site, covariate) |>
  arrange(cancer_site) |>
  select(`Cancer site`=cancer_site, covariate, value, N, `1`, `2`, `3`, `4`) |>
  flextable() |>
  add_header_row(values = c("", "Stage"), colwidths = c(4,4)) |>
  merge_v(j=c(1,2)) |>
  add_table_theme() |>
  save_as_docx(path=paste0(output_dir, "desc_stage_cc.docx"))


total_label_cc <-
  summary_cc |>
  group_by(cancer_site, covariate, value) |>
  summarise(N = sum(n)) |>
  mutate(label = str_glue("N={N}")) |>
  ungroup()

p_cc <- summary_cc |>
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
  select(covariate, value, cancer_site, prop, lower, upper, stage_best, N) |>
  ggplot() +
  geom_bar(aes(x = value, y = prop, fill = stage_best), stat = "identity") +
  ylim(-12, 101)+
  geom_text(size = 3, hjust = 0, colour = "#5A5A5A", aes(x = value, y = -12, label = label), data = total_label_cc) +
  scale_fill_manual(values = rev(c("#5340ff", "lightblue", "#C47B7B", "darkred"))) +
  scale_x_discrete(limits = rev) +
  coord_flip() +
  facet_grid(covariate~cancer_site, scales = "free_y") +
  labs(x = "",
       y = "Proportion of Patients (%)") +
  guides(fill=guide_legend(reverse=T, title="Stage")) +
  theme_minimal()+
  theme(strip.text = element_text(size = 9, colour = "black"),
        axis.text.x = element_text(size = 8, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 9, colour = "black"),
        axis.title.x = element_text(size = 9, colour =  "black"),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())
ggsave(plot = p_cc, paste0(output_dir, "desc_cc.png"),
       width = 30, height = 15, unit = "cm")




# Desc not by stage ----
total_all <- df_cc |> group_by(cr_eventtype) |> # add total info
  summarise(n_all = n()) |>
  mutate(N_all = sum(n_all),
         prop_all = format_numb(100*n_all/N_all),
         `N (%)` = str_glue("{n_all} ({prop_all}%)"),
         covariate = "Total",
         value = "Total")

summary_cc_all <- 
  df_cc_long |> 
  left_join(total_all |> select(cr_eventtype, N_all=n_all)) |>
  group_by(cr_eventtype, N_all, covariate, value) |>
  summarise(n_all = n()) |>
  mutate(prop_all = format_numb(100*n_all/N_all),
         `N (%)` = str_glue("{n_all} ({prop_all}%)")) |>
  ungroup() |>
  # get symptom summary
  bind_rows(df_cc |>
              left_join(total_all |> select(cr_eventtype, N_all=n_all)) |>
              mutate(Neither = sx_eventtype_14 == 0 & sx_eventtype_4 == 0,
                     `Rectal Bleeding` = as.logical(as.numeric(sx_eventtype_14) - 1),
                     CIBH = as.logical(as.numeric(sx_eventtype_4) -1)) |>
              pivot_longer(cols = c(`Rectal Bleeding`, CIBH, Neither), names_to = "symptom") |>
              group_by(cr_eventtype, N_all, symptom) |>
              summarise(n_all = sum(value)) |>
              mutate(prop_all = format_numb(100*n_all/N_all),
                     `N (%)` = str_glue("{n_all} ({prop_all}%)")) |>
              rename(value = symptom) |>
              mutate(covariate = "Symptom") |>
              ungroup()) |>
  bind_rows(total_all) |>
  mutate(covariate = factor(covariate, c("Gender", "Age", "IMD (2015)", "Symptom", "Total"))) |>
  arrange(covariate) |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  pivot_wider(id_cols = c(covariate, value),
              names_from = cancer_site,
              values_from = `N (%)`)

summary_cc_all |>
  flextable() |>
  merge_v(j=c(1)) |>
  add_table_theme() |>
  save_as_docx(path=paste0(output_dir, "desc_overall_cc.docx"))

rm(list = ls())
gc()


# Patient Characteristics by origin of cohort (symptom cohort vs random sample) ----
# if in both --> treat as from random sample
# Load libraries 
if(!require('pacman'))install.packages('pacman')
pacman::p_load(dplyr, tidyr, stringr,  lubridate, gt, ggplot2, purrr, broom, flextable)

source("./scripts/00_functions.R") #load general helper functions
source("./scripts/00_global_var.R") #load global variables
source("./scripts/00_cancer_stage_dataprep.R") #import data from MySQL

## Output directory
output_dir <- "./results/cancer_stage_cohort/"
if(!dir.exists(output_dir)) { dir.create(output_dir) }
colnames(df)

df_cc <- df |> filter(stage_best != -99) 

df_cc <- df_cc |> mutate(cohort_origin = ifelse(random_sample == 1,
                                                "random sample",
                                                ifelse(random_sample == 0 & cohort == 1,
                                                       "symptom cohort", "")))
ftable(xtabs(~random_sample+cohort+cohort_origin, data = df_cc)) #all correct locations

df_cc_long <-
  df_cc |>
  rowwise() |>
  mutate(cancer_site = names(which(cancer_codes == cr_eventtype))) |>
  ungroup() |> #undo 'rowwise'
  recode_gender() |>
  mutate(
    EP = final_route == "Emergency presentation",
    advanced = stage_bin == "advanced",
    imd2015_5 = ifelse(imd2015_5 == 1, "1 - Least", ifelse(imd2015_5 == 5, "5 - Most", imd2015_5))) |>
  rename(Gender = gender,
         `Emergency presentation` = EP,
         `Cancer site` = cancer_site,
         `Advanced stage` = advanced,
         `IMD (2015)` = imd2015_5) |>
  mutate(across(everything(), as.character)) |>
  pivot_longer(cols = c(Gender, `Emergency presentation`, `Cancer site`, `Advanced stage`,`IMD (2015)`), 
               names_to = "covariate")

summary_compare <- 
  df_cc_long |> 
  group_by(cohort_origin, covariate, value) |> #categorical var
  summarise(n_all = n()) |>
  mutate(N_all = sum(n_all),
         prop_all = format_numb(100*n_all/N_all),
         `N (%)` = str_glue("{n_all} ({prop_all}%)")) |>
  ungroup() |>
  bind_rows(df_cc |> #age
              group_by(cohort_origin) |>
              summarise(median_age = format_numb(median(age_cr)),
                        iqr_age = format_numb(IQR(age_cr)),
                        `N (%)` = str_glue("{median_age} ({iqr_age})")) |>
              select(cohort_origin, `N (%)`) |>
              mutate(covariate = "Age",
                     value = "Median (IQR)") |>
              ungroup()) |>
  bind_rows(df_cc |> group_by(cohort_origin) |> # add total info
              summarise(n_all = n()) |>
              mutate(N_all = sum(n_all),
                     prop_all = format_numb(100*n_all/N_all),
                     `N (%)` = str_glue("{n_all} ({prop_all}%)"),
                     covariate = "Total",
                     value = "Total")) |>
  mutate(covariate = factor(covariate, 
                            c('Age', 'Gender', 'IMD (2015)', 
                              'Cancer site', 'Advanced stage', 'Emergency presentation',
                              'Total'))) |>
  arrange(covariate) |>
  pivot_wider(id_cols = c(covariate, value),
              names_from = cohort_origin,
              values_from = `N (%)`)

chisq <- df_cc_long |> 
  group_by(covariate) |>
  nest() |>
  mutate(chisq = map(data, ~chisq.test(table(.x$value, .x$cohort_origin))),
         stat = map(chisq, glance)) |>
  unnest(stat) |>
  select(-data, -chisq) |>
  arrange(p.value)

wilcox_age <- wilcox.test(age_cr ~ cohort_origin, data = df_cc)

summary_compare <-
  summary_compare |>
  left_join(
    bind_rows(chisq |> select(covariate, p.value, method),
              data.frame(covariate = "Age", p.value = wilcox_age$p.value, method="wilcoxon-rank sum test")),
    by = 'covariate'
  )
write.csv(summary_compare, paste0(output_dir, "compare_cohort_desc.csv"), row.names = F)

summary_compare |>
  mutate(p.value = round(p.value, 3)) |>
  select(-method) |>
  flextable() |>
  merge_v(c(1,5)) |>
  hline(i=c(1,3,8,10,12,14)) |>
  fontsize(part = "all", size = 9) |>
  bold(part = "header") |>
  font(part = "all", fontname =  "Calibri") |>
  save_as_docx(path = paste0(output_dir, "compare_cohort_desc.docx"))




