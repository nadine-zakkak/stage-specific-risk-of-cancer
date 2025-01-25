## Purpose: Prepare table for cancer incidence analysis, using rectal bleeding & change in bowel habit cohorts
## Author:  Nadine Zakkak 
## Adopted from Becky White

# Load libraries ----
if(!require('pacman'))install.packages('pacman')
pacman::p_load(RMySQL, dplyr, ggplot2, stringr, data.table, tidyr, readr, forcats, lubridate)

# Load global variables
source("./scripts/00_global_var.R")

# set up MySQL connection ---- 
# db <- dbConnect(MySQL(), host = "", 
#                 user = "", password = "",
#                 dbname = "", port = 3306)
dbListTables(db)

# Import data from MySQL ---- 
rs <- dbSendQuery(db, "select * from e2_cohort_allsx")
allsx <- fetch(rs, -1)
rs1 <- dbSendQuery(db, "select * from e2_cohort_allsx_crc")
allcrc <- fetch(rs1, -1)
rs2 <- dbSendQuery(db, "select * from e2_cohort_allsx_othercancer")
othercancer <- fetch(rs2, -1)

## Empirical checks that data loaded correctly
length(unique(allsx$e_patid)) 
nrow(allsx) 
allsx |> group_by(eventtype) |> summarise(row = n(), pat = length(unique(e_patid))) 

min(allsx$eventdate); max(allsx$eventdate)

length(unique(allcrc$e_patid))
nrow(allcrc) 
min(allcrc$eventdate); max(allcrc$eventdate) 

length(unique(othercancer$e_patid)) 
nrow(othercancer) 
min(othercancer$eventdate); max(othercancer$eventdate) 

# Clean up data for R ----
allsx$eventdate <- as.Date(allsx$eventdate)
allcrc$eventdate <- as.Date(allcrc$eventdate)
othercancer$eventdate <- as.Date(othercancer$eventdate)

patient_info <- 
  allsx |>
  distinct(e_patid, dob, uts, crd, crd_oneyear, lcd, tod, deathdate, 
           age30_date, age100_date, imd2015_10, gender)
df <- allsx |> 
  select(-c(dob, uts, crd, crd_oneyear, lcd, tod, deathdate, 
            age30_date, age100_date, imd2015_10, gender)) |>
  bind_rows(allcrc) # join symptoms records with crc records 
df |> group_by(sx) |> summarise(row = n(), pat = length(unique(e_patid)))

head(df)
str(df)
df_orig <- df

patients_numb <- data.frame()

patients_numb <- patients_numb |> bind_rows(
  data.frame(description = "rectal bleed and/or CIBH event (1/1/2007-31/12/2017)",
             brief_desc = "RB/CIBH in 1/1/2007-31/12/2017",
             patients = length(unique(allsx$e_patid)),
             rb_patients = length(unique(allsx$e_patid[which(allsx$eventtype == 14)])),
             cibh_patients = length(unique(allsx$e_patid[which(allsx$eventtype == 4)]))))

patients_numb <- patients_numb |> bind_rows(
  data.frame(description   = "rectal bleed and/or CIBH event in max(uts, crd+1y, >=30yo) --> min(lcd, tod, dod, <100yo)",
             brief_desc    = "potential eligible sx",
             patients      = length(unique(allsx$e_patid[which(allsx$potent_elig == 1)])),
             rb_patients   = length(unique(allsx$e_patid[which(allsx$potent_elig == 1 & allsx$eventtype == 14)])),
             cibh_patients = length(unique(allsx$e_patid[which(allsx$potent_elig == 1 & allsx$eventtype == 4)]))))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# F) All sx + cancer records were sorted chronologically and eligible 'rectal bleeding'/CIBH events were flagged. ----
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# ‘Potentially eligible’ rectal bleeding/CIBH records were flagged as ‘eligible’ records if there was no previous colon&rectal cancer diagnosis within a year
# sort chronologically and carry forward colon&rectal cancer event dates by patient
df <-
  df |> 
  group_by(e_patid) |> 
  arrange(e_patid, eventdate, origin_id) |> 
  mutate(eventtype_lag = lag(eventtype),
         sx_lag = lag(sx)) |> 
  mutate(prev_crc_date = ifelse(sx_lag == 0 &
                                  eventtype_lag %in% cancer_codes, 
                                lag(as.character(eventdate)), NA)) |> 
  fill(prev_crc_date, .direction = "down") |> # carry forward previous event dates for CRC
  mutate(prev_crc_1yr = 
           !is.na(as.double(eventdate - as.Date(prev_crc_date))) 
         & as.double(eventdate - as.Date(prev_crc_date)) < 365) # CRC less than a year before sx?

# flag eligible sx records
df <-
  df |> 
  mutate(sx_eligible = 
           sx == 1 
         & eventtype %in% symptom_codes 
         & !is.na(potent_elig) 
         & potent_elig == 1 
         & prev_crc_1yr == 0)

length(which(df$sx_eligible == 1))
length(unique(df$e_patid[which(df$sx_eligible == 1)])) 

table(df$sx_eligible, df$eventtype)

patients_numb <- patients_numb |> bind_rows(
  data.frame(description   = "rectal bleed and/or CIBH event w/o CRC in prev year",
             brief_desc    = "eligible sx",
             patients      = length(unique(df$e_patid[which(df$sx_eligible == 1)])),
             rb_patients   = length(unique(df$e_patid[which(df$eventtype == 14 & df$sx_eligible == 1)])),
             cibh_patients = length(unique(df$e_patid[which(df$eventtype == 4 & df$sx_eligible == 1)]))))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# G) ‘Eligible’ sx records were sorted again chronologically, and a random event was flagged as the ‘index sx event’ ----
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#
sx_df <-
  df |> filter(sx == 1 & eventtype %in% symptom_codes) 
nrow(sx_df) # expected 286998

rm(df)

# randomly choose an eligible symptom to be index event
set.seed(171420)

sx_eligible_df <-
  sx_df |>
  filter(sx_eligible == 1) |>
  group_by(e_patid) |>
  mutate(random_numb = runif(row_number(), 1, row_number())) |> # generate random numbers for eligible sx per patient
  mutate(index_sx = random_numb == sample(random_numb, size = 1)) # randomly choose a number and set as index sx event

# number of 'index' events should now be equal to the number of patients who had >=1 'eligible' symtpom
length(which(sx_eligible_df$index_sx == 1)) == length(unique(df$e_patid[which(df$sx_eligible == 1)]))
length(which(sx_eligible_df$index_sx == 1)) 
sx_eligible_df |> filter(index_sx) |> count(e_patid) |> filter(n > 1) # each patient had exactly 1 index symptom


patients_numb <- patients_numb |> bind_rows(
  data.frame(description   = "Single random 'eligible' sx for each patient",
             brief_desc    = "1 'index' event for each patient",
             patients      = length(unique(sx_eligible_df$e_patid[which(sx_eligible_df$index_sx==1)])),
             rb_patients   = length(unique(sx_eligible_df$e_patid[which(sx_eligible_df$index_sx==1 & sx_eligible_df$eventtype==14)])),
             cibh_patients = length(unique(sx_eligible_df$e_patid[which(sx_eligible_df$index_sx==1 & sx_eligible_df$eventtype==4)]))))

# join with all sx events
final_df <-
  sx_df |> 
  filter(sx_eligible == 0) |> 
  bind_rows(sx_eligible_df) |> 
  arrange(e_patid, eventdate, origin_id)

length(unique(final_df$e_patid)) == length(unique(sx_df$e_patid))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# H) Flag the latest symptoms that occurred before/on index date, per symptom ----
# -- Exclusions:
# --- Same symptom as index symptom recorded on index date (i.e index sx recorded >1 on index date)
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --#
# Fill in index date/eventtype for all sx per patient
final_df <-
  final_df |> 
  group_by(e_patid) |> 
  arrange(e_patid, desc(index_sx), eventdate) |> # index event first (desc(index)) then sort by eventdate
  mutate(index_date = ifelse(index_sx == 1, as.character(eventdate), NA), 
         index_eventtype = ifelse(index_sx == 1, eventtype, NA)) |>
  fill(index_date, index_eventtype) # fill in index info 

length(unique(final_df$e_patid)) == length(unique(sx_df$e_patid))
min(final_df$index_date, na.rm = T); max(final_df$index_date,  na.rm = T) # some NA bcz not all pat had eligible index event

final_df$index_date <- as.Date(final_df$index_date)
min(final_df$index_date, na.rm = T); max(final_df$index_date,  na.rm = T) # expected 2007-01-01, 2017-12-30

final_df_backup <- final_df

sx_before <-
  final_df |>
  group_by(e_patid) |>
  mutate(sx_before_index =
           !(index_eventtype == eventtype & index_date == eventdate) 
         & eventdate <= index_date) |> # sx occurred before or on index date ?
  ungroup() |>
  group_by(e_patid, eventtype, sx_before_index) |>
  mutate(sx_latest = 
           sx_before_index == 1 
         & eventdate == max(eventdate)) |> # latest event before index date 
  filter(sx_latest == 1) |>
  ungroup() |>
  select(e_patid, eventtype, eventdate) |>
  distinct()

final_df <-
  final_df |> 
  filter(index_sx == 1) |> 
  select(e_patid, prev_crc_date, index_date, index_eventtype) |>
  left_join(sx_before, by = c("e_patid" = "e_patid")) |> 
  ungroup() |> 
  as_tibble() |>
  pivot_wider(id_cols= c(e_patid, prev_crc_date, index_date, index_eventtype), 
              names_from = eventtype, 
              values_from = eventdate,
              names_prefix = "sx_") |>
  select(-'sx_NA')

nrow(final_df) == length(unique(final_df$e_patid))
nrow(final_df)
nrow(final_df) == length(unique(sx_eligible_df$e_patid[sx_eligible_df$index_sx == 1]))
min(final_df$index_date); max(final_df$index_date)

# clean data frame
final_df <-
  final_df |> 
  select(e_patid, prev_crc_date, index_date, index_eventtype, matches("sx")) |> 
  as_tibble()

rm(df_orig, sx_df, sx_eligible_df, sx_before)

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# I) Fill in index dates for CRC records by patient ---- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
allcrc <- allcrc |> 
  select(-sx) |>
  rename(cancer_origin_id = origin_id,
         cancer_eventtype = eventtype,
         cancer_eventdate = eventdate) ## rename columns to make it specific for cancer

final_df_backup <- final_df

# join index sx events with all CRC events
final_df <-
  final_df |>
  left_join(allcrc, by = c("e_patid" = "e_patid")) |>
  arrange(e_patid, desc(index_date), cancer_eventdate) |>
  group_by(e_patid) |>
  fill(index_date, index_eventtype) |>
  as_tibble()

length(unique(final_df$e_patid))
min(final_df$index_date); max(final_df$index_date)  

final_df$index_date <- as.Date(final_df$index_date)
min(final_df$index_date); max(final_df$index_date)  

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# J) First 'eligible' CRC outcome(s) was(were) flagged and the number of 'eligible' CRC events were counted ----
# 'eligible' CRC: within a year following index date
# # treat combined (find earliest of colon & rectal cancers)
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
final_df_backup <- final_df

final_df <-
  final_df |> 
  mutate(outcome_elig = between(as.double(cancer_eventdate - index_date), 1, 365)) |> # eligible CRC outcome (event within a year after index date)
  group_by(e_patid, outcome_elig) |>
  arrange(e_patid, cancer_eventdate) |> # order by eventdate
  mutate(crc_outcome = 
           outcome_elig == 1 
         & cancer_eventdate == min(cancer_eventdate)) |> # choose earliest eligible CRC outcome !!can have > 1 per patient!!
  ungroup() |> 
  group_by(e_patid) |> 
  mutate(count_crc_1y = sum(outcome_elig, na.rm = T), # all eligible crc events within 1y of index
         count_crc_sync = sum(crc_outcome, na.rm = T)) |> # all earliest crc events (on same day)
  ungroup() |>
  as_tibble()

# Split final_df into 2 dataframes: index sx and outcome
index_df <- final_df |> distinct(e_patid, prev_crc_date, index_date, index_eventtype, sx_14, sx_4, count_crc_1y, count_crc_sync)
length(unique(index_df$e_patid)) == nrow(index_df)

outcome_df <- final_df |> filter(crc_outcome == 1) |>
  select(e_patid, cancer_origin_id, cancer_eventdate, stage_best, site_icd10_o2, final_route, cancer_eventtype, outcome_elig, crc_outcome)# 3970 combined colon & rectal cancers
length(unique(outcome_df$e_patid))
nrow(outcome_df)

# explore multiple outcomes:
## 1) > 1 colon cancers on same day
## 2) > 1 rectal cancers on same day
## 3) colon & rectal cancers on same day 
mult_outcome <- 
  outcome_df |> 
  filter(e_patid %in% unique(outcome_df$e_patid[duplicated(outcome_df$e_patid)])) 
length(unique(mult_outcome$e_patid)) 
table(mult_outcome$e_patid, mult_outcome$cancer_eventtype)

mult_outcome <- 
  mult_outcome |>
  group_by(e_patid, cancer_eventtype) |>
  summarise(n = n()) |>
  mutate(N = sum(n)) |> 
  mutate(mult_colon = ifelse(cancer_eventtype == 11 & n > 1, 1, 0),
         mult_rec = ifelse(cancer_eventtype == 12 & n > 1, 1, 0), 
         col_rec = ifelse(n < N, 1, 0))
length(unique(mult_outcome$e_patid[which(mult_outcome$col_rec == 1)])) 
length(unique(mult_outcome$e_patid[which(mult_outcome$mult_colon == 1)])) 
length(unique(mult_outcome$e_patid[which(mult_outcome$mult_rec == 1)]))


# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# K) Choose 1 'eligible' CRC outcome if >1 CRC on same day ----
# 'eligible' CRC: within a year following index date
# conditions: 
# 1) Choose CRC with the highest stage
# 2) If no stage recorded or same stage --> choose randomly
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Stage recode
table(outcome_df$stage_best, useNA = "always")
outcome_df$stage_best <- parse_number(outcome_df$stage_best, na=c("?", "U", "")) #only extract number
table(outcome_df$stage_best, useNA = "always")
outcome_df$stage_best[outcome_df$cancer_eventtype %in% cancer_codes & is.na(outcome_df$stage_best)] <- -99 #assign with missing stage to -99
table(outcome_df$stage_best, useNA = "always")

set.seed(160644)

outcome_df <- 
  outcome_df |> 
  group_by(e_patid) |>
  slice(which.max(rank(stage_best, ties.method = "random", na.last = FALSE))) |> # (NA will always rank 1)
  as_tibble()

nrow(outcome_df)
nrow(outcome_df) == length(unique(outcome_df$e_patid)) 

# join index sx and single CRC outcome
final_df <- 
  index_df |>
  left_join(outcome_df, by = c("e_patid"="e_patid"))

nrow(final_df) == length(unique(final_df$e_patid)) 

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# L) Fill in index dates for all other cancer events ---- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
othercancer_elig <-
  othercancer |> 
  filter(e_patid %in% unique(final_df$e_patid))|> # only for patients with eligible sx
  select(-sx, -stage_best) |>
  rename(cancer_origin_id = origin_id,
         cancer_eventtype = eventtype,
         cancer_eventdate = eventdate) ## rename columns to make it specific for cancer

nrow(othercancer_elig)

# join other cancers to index sx and outcome rc
final_df_cancer <-
  final_df |>
  bind_rows(othercancer_elig)

head(final_df_cancer)
str(final_df_cancer)
nrow(final_df_cancer)
length(unique(final_df_cancer$e_patid))

final_df_cancer <-
  final_df_cancer |> 
  arrange(e_patid, desc(index_date), cancer_eventdate) |> # index event first then sort by date
  group_by(e_patid) |> 
  fill(index_date, index_eventtype, matches("sx_.*", perl = T), prev_crc_date) |> 
  arrange(e_patid, desc(index_date), cancer_eventdate) |> 
  as_tibble()

min(final_df_cancer$index_date); max(final_df_cancer$index_date) 

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- ---
# M) count number of `other` cancer within +-1year of index event ---- 
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
final_df_cancer_back <- final_df_cancer
final_df_cancer <-
  final_df_cancer |> 
  mutate(cancer_after_sx = 
           !(cancer_eventtype %in% cancer_codes) 
         & cancer_eventdate >= index_date
         & as.double(cancer_eventdate - index_date) <= 365, # cancer event (except CRC) within a year after index event
         cancer_before_sx = 
           !(cancer_eventtype %in% cancer_codes)
         & cancer_eventdate < index_date
         & as.double(index_date - cancer_eventdate) <= 365) |> # cancer event (except CRC) within a year before index event
  group_by(e_patid) |> 
  mutate(count_cancer_after_1y = sum(cancer_after_sx, na.rm = T), # count all cancers except CRC within a year following index RB
         count_cancer_before_1y = sum(cancer_before_sx, na.rm = T)) |>  # count all cancers except CRC a year before index RB
  fill(count_crc_1y, count_crc_sync, count_cancer_after_1y, count_cancer_before_1y, .direction = "downup") |>
  ungroup() |>
  as_tibble()

final_df_cancer |> filter(cancer_before_sx | cancer_after_sx | crc_outcome) |> nrow() 
final_df_cancer |> filter(cancer_before_sx | cancer_after_sx | crc_outcome) |> distinct(e_patid) |> nrow() 
final_df_cancer |> filter(cancer_before_sx | cancer_after_sx | crc_outcome) |> distinct(cancer_origin_id) |> nrow() 
final_df_cancer |> filter(cancer_before_sx | cancer_after_sx | crc_outcome) |> count(cancer_origin_id) |> filter(n > 1)


# Remove any patid that have diff patid but same tid (i.e. same patient but switched prac)
final_df_cancer <- 
  final_df_cancer |>
  filter(!e_patid %in% 
           (final_df_cancer |>
              filter(cancer_origin_id %in% 
                       (final_df_cancer |> 
                          filter(cancer_before_sx | cancer_after_sx | crc_outcome) |> 
                          count(cancer_origin_id) |> 
                          filter(n > 1) |>
                          pull(cancer_origin_id))) |> # get tid that occurs >1
              distinct(e_patid) |> 
              pull())) # get patid related to the tid

nrow(final_df_cancer) 
length(unique(final_df_cancer$e_patid))


patients_numb <- patients_numb |> bind_rows(
  data.frame(description   = "Patients that did not share the same tumour id",
             brief_desc    = "Patients that did not share the same tid",
             patients      = length(unique(final_df_cancer$e_patid)),
             rb_patients   = length(unique(final_df_cancer$e_patid[which(final_df_cancer$index_eventtype == 14)])),
             cibh_patients = length(unique(final_df_cancer$e_patid[which(final_df_cancer$index_eventtype == 4)]))))

# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# N) Create final comprehensive dataset ----
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# Final comprehensive datasest
full <- 
  final_df_cancer |>
  distinct(e_patid,
           index_date, index_eventtype,
           prev_crc_date, sx_4, sx_14,
           count_crc_1y, count_crc_sync, 
           count_cancer_after_1y, count_cancer_before_1y) |> # extract index sx with info
  inner_join(patient_info) |> # add patient info
  left_join(outcome_df |>  select(e_patid, cancer_eventdate, stage_best, site_icd10_o2, 
                                  final_route, cancer_eventtype), by = c("e_patid")) |> # add CRC outcome info
  rename(indexdate = index_date, 
         first_crc_date = cancer_eventdate,
         stage = stage_best,
         first_crc_desc = site_icd10_o2,
         first_crc_route = final_route,
         cancer_site = cancer_eventtype)

nrow(full)
nrow(full) == length(unique(full$e_patid))
min(full$indexdate); max(full$indexdate) 
min(full$first_crc_date, na.rm = T); max(full$first_crc_date, na.rm = T)

# Exclude patients with missing IMD ----
table(full$imd2015_10, useNA = "always")
length(which(is.na(full$imd2015_10)))
full <-
  full |>
  filter(!is.na(imd2015_10))
length(which(is.na(full$imd2015_10)))

nrow(full) 
nrow(full) == length(unique(full$e_patid)) 

patients_numb <- patients_numb |> bind_rows(
  data.frame(description   = "Patients with missing IMD",
             brief_desc    = "Patients with missing iMD",
             patients      = length(unique(full$e_patid)),
             rb_patients   = length(unique(full$e_patid[which(full$index_eventtype == 14)])),
             cibh_patients = length(unique(full$e_patid[which(full$index_eventtype == 4)]))))

# Derive new variables ----
full <- full |> mutate(crc = cancer_site %in% cancer_codes)
full$age_sx <- as.integer(as.Date(full$indexdate) - as.Date(full$dob))/365
full$days_to_crc <- as.double(as.Date(full$first_crc_date) - as.Date(full$indexdate))

full <- full |> mutate(imd2015_5 = ceiling(imd2015_10/2))

full <- full |> 
  mutate(age10_cat = case_when(
    full$age_sx < 40 ~ "30-39",
    full$age_sx < 50 ~ "40-49",
    full$age_sx < 60 ~ "50-59",
    full$age_sx < 70 ~ "60-69",
    full$age_sx < 80 ~ "70-79",
    full$age_sx < 90 ~ "80-89",
    full$age_sx >= 90 ~ "90+"
  ),
  age10_cat = factor(age10_cat))

full <- 
  full |> 
  mutate(stage_bin = as.factor(case_when(
    full$stage %in% c('1', '2') ~ "non-advanced",
    full$stage %in% c('3', '4') ~ "advanced",
    full$stage == "-99" ~ "missing")),
    stage_bin = fct_relevel(stage_bin, "non-advanced", "advanced", "missing"))

str(full)

# Patient flowchart file ----
patients_numb <- patients_numb |>
  mutate(patients_diff = patients - lag(patients),
         rb_diff = rb_patients - lag(rb_patients),
         cibh_diff = cibh_patients - lag(cibh_patients)) 

write.csv(patients_numb, "./results/cohort_flowchart.csv", row.names = F)

# Extract back to MySQL ----
dbWriteTable(conn = db,
             name = "e2_cohort_allsx_crc_final_dataset",
             value = full |> mutate(across(where(is.logical), as.integer)), overwrite = TRUE, row.names = FALSE)

dbDisconnect(db)

