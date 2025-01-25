# Purpose: Load and reshapee cancer registry data for analysis
# Author: Matthew Barclay

library(readxl)
library(dplyr)
library(tidyr)

orig_data <- 
  read_excel("Counts by single year age, sex, deprivation, stage for 2018-19 combined.xlsx", sheet = 2)

colnames(orig_data) <- tolower(colnames(orig_data))

summary(orig_data)
table(orig_data$cancer)
summary(orig_data$age) #30-89
table(orig_data$stage) #1-2, 3, 4, X

orig_data <- orig_data |>
  mutate(stage34 = ifelse(stage == "1-2", 0, ifelse(stage %in% c("3", "4"), 1, NA)))


#only colon and rectal cancers
crc_data <- orig_data |> filter(cancer %in% c("Colon", "Rectal"))

# Reshape for main analysis -----------------------------------------------
# MATT  
crc_data_mod <- crc_data |> 
  filter(!is.na(stage34)) |> #exclude missing stages
  select(cancer, male, imd19_quintile_lsoas, age, stage34, tumour_count) |>
  group_by(cancer, male, imd19_quintile_lsoas, age, stage34) |>
  summarise(n = sum(tumour_count)) |>
  # group_by() |> 
  pivot_wider(
    names_from = stage34,
    names_prefix = "n",
    values_from = n,
    id_cols = c(cancer, male, imd19_quintile_lsoas, age)
  ) |>
  rename(stage12 = n0, stage34 = n1)

crc_data_mod <- crc_data_mod |>
  mutate(age_c = (age - 60)/10) |> 
  select(cancer, male, imd19_quintile_lsoas, age, age_c, stage12, stage34) |>
  ungroup()


data <- as.data.frame(lapply(crc_data, rep, crc_data$tumour_count))
