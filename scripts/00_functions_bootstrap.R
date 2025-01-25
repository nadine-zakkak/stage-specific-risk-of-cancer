## Purpose: Functions to run bootstrap ##
## Author:  Nadine Zakkak ##

# Predict functions ----
model_predict <- function(model, pred_data){
  predictions(model,
              newdata = pred_data,
              vcov = "HC3") |> # robust standard error
    rowwise() |>
    mutate(index_labels = factor(names(symptom_codes)[symptom_codes == index_eventtype])) |>
    mutate(index_labels = fct_relevel(index_labels, 'Rectal Bleeding'))
}

model_stage_predict <- function(model, pred_data){
  predictions(model,
              newdata = pred_data,  
              vcov = "HC3" # robust standard error
  ) |>
    mutate(sx = factor(case_when(
      sx_eventtype_4 == 1 ~ "CIBH",
      sx_eventtype_14 == 1 ~ "Rectal Bleeding",
      sx_eventtype_14 == 0 ~ "Neither"
    ), 
    levels = c("Rectal Bleeding", "CIBH",  "Neither")))
}

# bootstrap functions -----
bootstrap <- function(df, reg_formula, groups, stage, pred_data) {
  bootstrap_model(df |> 
                    group_by(!!!groups) |> 
                    nest() |>
                    mutate(datarand = map(data, ~sample_n(.x, nrow(.x), replace = T))), 
                  reg_formula, stage, pred_data)
}

bootstrap_model <- function(df, reg_formula, stage, pred_data){
  bootstrap_predict(df |> 
                      mutate(glms = map(
                        datarand, 
                        ~glm(reg_formula,
                             family = binomial(link = "logit"),
                             data = .x)
                      )), stage, pred_data)
}

bootstrap_predict <- function(df, stage, pred_data){
  if(!stage){
    df |> 
      mutate(predictions = map(
        glms,
        ~model_predict(., pred_data)
      )) |>
      select(gender, predictions) |>
      unnest(predictions) |>
      select(gender, index_labels, age_decade, predicted) |>
      mutate(age = 60+age_decade*10) |>
      rowwise() |>
      mutate(gender = ifelse(gender==1, "Men", "Women")) |>
      ungroup()
  }else{
    df |> 
      mutate(predictions = map(
        glms,
        ~model_stage_predict(., pred_data)
      )) |>
      select(cr_eventtype, gender, predictions) |>
      unnest(predictions) |>
      select(cr_eventtype, gender, sx, age_decade, predicted) |>
      mutate(age = 60+age_decade*10) |>
      rowwise() |>
      mutate(cancer_site = names(cancer_codes[which(cancer_codes == cr_eventtype)]),
             gender = ifelse(gender == 1, "Men", "Women")) |>
      ungroup()
  }
}
