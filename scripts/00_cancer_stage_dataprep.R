## Purpose: Import and prepare cancer incidence data for stage analyses ##
## Author:  Nadine Zakkak ##

# Add libraries here
if(!require('pacman'))install.packages('pacman')
pacman::p_load(RMySQL, forcats, dplyr, lubridate)

############
## Stage of colon and rectal cancers e2 data preparation
##
## Created by Nadine Zakkak 1/Nov/2022
############
# set up MySQL connection ---- 
# db <- dbConnect(MySQL(), host = "", 
#                 user = "", password="", 
#                 dbname = "", port = 3306)
rs <- dbSendQuery(db, "select * from e2_allcrc_sx_final_dataset_new")
df <- fetch(rs, -1)


# Clean up df for R ----
## date columns
date_cols <- c('cr_eventdate', 'dob', 'uts', 'crd', 'lcd', 'tod', 'deathdate', 'age30_date', 'age100_date', 'crd_oneyear')
df <- df |> mutate(across(all_of(date_cols), as.Date))
## Categorical columns
fact_cols <- c("cr_eventtype", "stage_best", "stage_bin", "gender", 'imd2015_10', "imd2015_5", "age10_cat", "sx", "sx_eventtype_4", "sx_eventtype_14")
df <- df|> mutate(across(all_of(fact_cols), as.factor))
df$stage_best <- fct_relevel(df$stage_best, "1", "2", "3", "4", "-99")
df$stage_bin <- fct_relevel(df$stage_bin, "early", "advanced", "missing")
df$imd2015_5 <- fct_relevel(df$imd2015_5, '1', '2', '3', '4', '5')
df$sx <- fct_relevel(df$sx, "Both", "Rectal Bleeding", "CIBH", "Neither")

rm(date_cols, fact_cols)

## prepare data for modelling ----
df_backup <- df
df <- 
  df |>
  mutate(age_decade = (age_cr - 60)/10) #centre age at 60yo and present as decades

# Complete case analysis ----
df_cc <- df |> 
  filter(stage_best != -99) |> # complete case - excl. records with missing stage info  
  mutate(advanced = stage_bin == "advanced")

## prediction dataframe ----
stage_predict_data <-
  expand.grid(age_decade=seq(-3,3,length.out=100),
              sx_eventtype_4 = factor(c(0,1)),
              sx_eventtype_14 = factor(c(0,1))) |>
  filter(!(sx_eventtype_4 == 1 & sx_eventtype_14 == 1))

dbDisconnect(db)





