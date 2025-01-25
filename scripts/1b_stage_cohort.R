## Purpose: Prepare table for cancer stage analysis, using colon and rectal cancer cohorts
## Author:  Nadine Zakkak 

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(RMySQL, dplyr, lubridate, fastDummies, readr)

# Load global variables
source("./scripts/00_global_var.R")

# helper function -----
## Fill in the patient number data frame
get_pat_numb <- function(desc, brief, df, eventtype_col){
  return(data.frame(description = desc,
             brief_desc = brief,
             patients = length(unique(df$e_patid)),
             colon_patients = length(unique(df$e_patid[which(df[[eventtype_col]]==11)])),
             rectal_patients = length(unique(df$e_patid[which(df[[eventtype_col]]==12)])),
             tumours = nrow(df)))
}

# set up MySQL connection ---- 
# db <- dbConnect(MySQL(), host = "", 
#                 user = "", password = "",
#                 dbname = "", port = 3306)
# dbListTables(db)

# Import data from MySQL ---- 
rs <- dbSendQuery(db, "select * from e2_allcrc")
allcrc <- fetch(rs, -1)
rs1 <- dbSendQuery(db, "select * from e2_allcrc_initial_sx_new")
allsx <- fetch(rs1, -1)

## Empirical checks that data loaded correctly
length(unique(allcrc$e_patid)) 
nrow(allcrc)
allcrc |> group_by(cancer_site_desc) |> summarise(row = n(), pat = length(unique(e_patid))) 
min(allcrc$eventdate); max(allcrc$eventdate) 

length(unique(allsx$e_patid)) 
nrow(allsx) 
allsx |> group_by(eventtype) |> summarise(row = n(), pat = length(unique(e_patid))) 
min(allsx$eventdate); max(allsx$eventdate)

allcrc <- allcrc |> rename(old_elig = elig)

## date columns
date_cols <- c('eventdate',
               'dob', 'uts', 'crd', 'crd_oneyear', 'tod', 'lcd', 'deathdate', 
               'age30_date', 'age100_date')
allcrc <- allcrc |> mutate(across(all_of(date_cols), as.Date))
allsx$eventdate <- as.Date(allsx$eventdate)

patients_numb <- data.frame()

patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer (1/1/2012-31/12/2018)",
               brief = "CRC in 1/1/2012-31/12/2018",
               allcrc, 'eventtype'))

#///////////////////////////////////////////////////////
# Step 3: Apply inclusion/exclusion criteria for CRC ----
#///////////////////////////////////////////////////////

## a) Age 31-99 years old ---- 
allcrc <- allcrc |> filter(eventdate >= age30_date %m+% years(1), eventdate < age100_date)

patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer after or at 31 years old and before 100 years old",
               brief = "CRC after >=31 years, <100yo",
               allcrc, 'eventtype'))

## b) 1 year after up to standard date ----
allcrc <- allcrc |> filter(eventdate >= uts %m+% years(1))

patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer after 1 year of uts",
               brief = "CRC after uts+1y",
               allcrc, 'eventtype'))

## c) 2 years after current registration date ----
allcrc <- allcrc |> filter(eventdate >= crd_oneyear %m+% years(1))
patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer after 2 years after crd",
               brief = "CRC after crd+2y",
               allcrc, 'eventtype'))

## d) Before last collection date ----
allcrc <- allcrc |> filter(eventdate <= lcd)
patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer before or on last collection date",
               brief = "CRC at <=lcd",
               allcrc, 'eventtype'))

## e) Before transfer out date, if any ----
allcrc <- allcrc |> filter(is.na(tod) |(eventdate <= tod))
patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer before or on transfer out date",
               brief = "CRC at <=tod",
               allcrc, 'eventtype'))

## f) Before death date, if died ----
allcrc <- allcrc |> filter(is.na(deathdate) |(eventdate <= deathdate))
patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc  = "colon and/or rectal cancer before or on death date",
               brief = "CRC at <=dod",
               allcrc, 'eventtype'))

#///////////////////////////////////////////////////////
# Step 4: Keep symptoms for patients with eligible CRC only ----
#///////////////////////////////////////////////////////

allsx <- allsx |> filter(e_patid %in% allcrc$e_patid)

# Empirical checks based on older code
length(unique(allsx$e_patid))
nrow(allsx)
allsx |> group_by(eventtype) |> summarise(row = n(), pat = length(unique(e_patid))) 

min(allsx$eventdate); max(allsx$eventdate)

patient_info <- 
  allcrc |>
  distinct(e_patid, dob, uts, crd, crd_oneyear, lcd, tod, deathdate, 
           age30_date, age100_date, imd2015_10, gender, random_sample, cohort)

#///////////////////////////////////////////////////////
# Step 5: Flag symptomatic CRC ----
#///////////////////////////////////////////////////////

## a) Join all CRC and symptoms into 1 dataset ----
df <- allcrc |> 
  select(-c(dob, uts, crd, crd_oneyear, lcd, tod, deathdate, 
            age30_date, age100_date, imd2015_10, gender, random_sample, cohort)) |>
  rename(cr_origin_id = origin_id,
         cr_eventdate = eventdate,
         cr_eventtype = eventtype) |>
  left_join(allsx |> rename(sx_eventdate = eventdate,
                            sx_origin_id = origin_id,
                            sx_eventtype = eventtype), by = c("e_patid"="e_patid"),
            multiple = "all") # join cancer records with sx records 
head(df)
str(df)
df_orig <- df
# Data preprocessing -----

## b) Flag symptomatic CRC -----
# if rectal bleeding and/or change in bowel habit occurred in the year prior to CRC diagnosis
df <- df |>
  mutate(prev_sx_1y = !is.na(cr_eventdate - sx_eventdate) &
           between(as.double(cr_eventdate - sx_eventdate), 0, 365)) # sx occurs within 1 yr prior to cancer?
nrow(df)
length(unique(df$e_patid))
length(unique(df$e_patid[which(df$prev_sx_1y)])) 

# remove the symptom date & event id and keep only the unique symptom events
df <- df |> 
  distinct(e_patid, cr_origin_id, cr_eventdate, final_route, stage_best, cr_eventtype, sx_eventtype, prev_sx_1y)
nrow(df)
length(unique(df$e_patid)) 
length(unique(df$e_patid[which(df$prev_sx_1y)]))

## c) Convert to wide dataset with 1 row per CRC ----
df_wide <- df |> 
  mutate(sx_eventtype = ifelse(is.na(sx_eventtype) | !prev_sx_1y, 0, sx_eventtype)) |> # set as 0 if no sx or sx not within 1y
  distinct(e_patid, cr_origin_id, cr_eventdate, stage_best, cr_eventtype, sx_eventtype, final_route) |>
  dummy_cols(select_columns = "sx_eventtype", remove_selected_columns = T) |>
  group_by(e_patid, cr_origin_id, cr_eventdate, cr_eventtype, stage_best, final_route) |> # for each CRC
  summarise(sx_eventtype_4 = max(sx_eventtype_4),
            sx_eventtype_14 = max(sx_eventtype_14)) |> # merge same CRC together
  mutate(symptomatic = sx_eventtype_4 == 1 | sx_eventtype_14 == 1) |> #flag if symptomatic CRC
  ungroup() |>
  distinct() |>
  arrange(e_patid)

df_wide |>
  count(cr_origin_id) |>
  filter(n > 1) 

df_wide |>
  count(e_patid) |>
  filter(n > 1) #  121 patients with > CRC tumour

length(unique(df_wide$e_patid))

#///////////////////////////////////////////////////////
# Step 6: Select only 1 CRC per patient -----
# Symptomatic CRC or randomly if all symptomatic or neither symptomatic
#///////////////////////////////////////////////////////
set.seed(142153)

df_wide_single <-
  df_wide |>
  group_by(e_patid) |>
  slice(which.max(rank(symptomatic, ties.method = "random"))) |> #choose symptomatic CRC or randomly 
  ungroup()

length(unique(df_wide_single$e_patid)) 
length(unique(df_wide_single$e_patid)) == nrow(df_wide_single)

patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc = "Single CRC per patient",
             brief  = "Single CRC per patient",
             df_wide_single, 'cr_eventtype'))


#///////////////////////////////////////////////////////
# Step 7: Create final comprehensive dataset ----
# with 1 CRC/patient and all the patient info
# 1) exclude patients with missing IMD
# 2) create some derived variables that will be used in analysis
#///////////////////////////////////////////////////////
final <- 
  df_wide_single |>
  left_join(patient_info, by = c("e_patid"="e_patid"))

## Exclude patients with missing IMD ----
table(final$imd2015_10, useNA = "always")
length(which(is.na(final$imd2015_10)))
final <-
  final |>
  filter(!is.na(imd2015_10))
length(which(is.na(final$imd2015_10)))

nrow(final)
nrow(final) == length(unique(final$e_patid)) 

patients_numb <- patients_numb |> bind_rows(
  get_pat_numb(desc = "Patients with missing IMD",
             brief  = "Patients with missing IMD",
             final, "cr_eventtype"))

## Derived variables ----
### Stage recode
table(final$stage_best, useNA = "always")
final$stage_best <- parse_number(final$stage_best, na=c("?", "U", "")) #only extract number
table(final$stage_best, useNA = "always")
final$stage_best[final$cr_eventtype %in% cancer_codes & is.na(final$stage_best)] <- -99 #assign with missing stage to -99
table(final$stage_best, useNA = "always")

final <-
  final |>
  mutate(stage_bin = factor(case_when(
    final$stage_best %in% c('1', '2') ~ "early",
    final$stage_best %in% c('3', '4') ~ "advanced",
    final$stage_best == "-99" ~ "missing"),
    levels = c("early", "advanced", "missing")))

final <- final |> 
  mutate(imd2015_5 = ceiling(imd2015_10/2),
         age_cr = as.double(as.Date(cr_eventdate) - as.Date(dob))/365)
final <- final |>
  mutate(age10_cat = factor(case_when(
    final$age_cr < 40 ~ "30-39",
    final$age_cr < 50 ~ "40-49",
    final$age_cr < 60 ~ "50-59",
    final$age_cr < 70 ~ "60-69",
    final$age_cr < 80 ~ "70-79",
    final$age_cr < 90 ~ "80-89",
    final$age_cr >= 90 ~ "90+"
  ),
  levels = c("30-39", "40-49", "50-59", "60-69", "70-79", "80-89", "90+")),
  sx = factor(case_when(
    final$sx_eventtype_4 == 1 & final$sx_eventtype_14 == 1 ~ "Both",
    final$sx_eventtype_4 == 1 ~ "CIBH",
    final$sx_eventtype_14 == 1 ~ "Rectal Bleeding",
    final$sx_eventtype_14 == 0 ~ "Neither"
  ), 
  levels = c("Both", "CIBH", "Rectal Bleeding", "Neither")))


# Patient flowchart file ----
patients_numb <- patients_numb |>
  mutate(patients_diff = patients - lag(patients),
         colon_diff    = colon_patients - lag(colon_patients),
         rectal_diff   = rectal_patients - lag(rectal_patients),
         tumours_diff  = tumours - lag(tumours)) 

write.csv(patients_numb, "./results/stage_cohort_flowchart.csv", row.names = F)

# Extract back to MySQL ----
dbWriteTable(conn = db,
             name = "e2_allcrc_sx_final_dataset_new",
             value = final |> mutate(across(where(is.logical), as.integer)), overwrite = TRUE, row.names = FALSE)

dbDisconnect(db)

