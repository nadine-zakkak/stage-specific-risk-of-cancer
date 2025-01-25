## Purpose: Import and prepare cancer incidence data for incidence analyses ##
## Author:  Nadine Zakkak ##

if(!require('pacman'))install.packages('pacman')
pacman::p_load(RMySQL, forcats, dplyr, lubridate)

# set up MySQL connection ---- 
print("Fetching data from MySQL")
# db <- dbConnect(MySQL(), host = "", 
#                 user = "", password = "", 
#                 dbname = "", port = 3306)
rs <- dbSendQuery(db, "select * from e2_cohort_allsx_crc_final_dataset")
df <- fetch(rs, -1)

# Clean up df for R ----
print("Cleaning up data")
# str(df)
attach(df)
## date columns
date_cols <- c('deathdate',
               'uts', 'lcd', 'crd', 'crd_oneyear', 'tod', 'dob', 'age30_date', 'age100_date',
               'indexdate', 'first_crc_date',
               'sx_14', 'sx_4')
df <- df|>mutate(across(date_cols, as.Date))
## Categorical columns
fact_cols <- c("gender", 'imd2015_10', "imd2015_5", "stage", "index_eventtype")
df <- df|> mutate(across(fact_cols, as.factor))
df$stage <- fct_relevel(df$stage, "1", "2", "3", "4", "-99")
df$imd2015_10 <- fct_relevel(df$imd2015_10, '1', '2', '3', '4', '5', '6', '7', '8', '9', '10')
df$imd2015_5 <- fct_relevel(df$imd2015_5, '1', '2', '3', '4', '5')
df$index_eventtype <- fct_relevel(df$index_eventtype, '14', '4') #Rectal bleeding to be the reference
detach(df)
## Add new column + recategorise age
attach(df)
df$died <- deathdate <= indexdate + years(1) 
df <- df |> 
  mutate(age10_cat = ifelse(age10_cat %in% c("30-39", "40-49"), "<50", 
                            ifelse(age10_cat %in% c("80-89", "90+"), "80+", 
                                   age10_cat)))
detach(df)

# str(df)

rm(date_cols, fact_cols)

## prediction data -----
incidence_predict_data <-  expand.grid(age_decade = seq(-3,3,length.out=100), 
                                       index_eventtype = factor(symptom_codes))

dbDisconnect(db)

