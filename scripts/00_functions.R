## Purpose: Common helper functions ##
## Author:  Nadine Zakkak ##

# Helper function ----

#' format_numb
#'
#' @param numb - number to format 
#'
#' @return formatted number to include 3 digits. if <10: 2 decimal places, (10,100): 1 decimal place and >100 0 decimal places
#'
#' @examples: format_numb(0.128) --> 0.13, format_numb(90.4736) --> 90.5, format_numb(123.2737) --> 123
format_numb <- function(numb) {
  
  numb = ifelse(numb<10,
                sprintf(numb, fmt='%.2f'),
                ifelse(numb<100,
                       sprintf(numb, fmt = '%.1f'),
                       sprintf(numb, fmt = '%1.0f')))
}

#' calculate_ci (95% CI) - from Matt
#'
#' @param n - number of samples (numerator)
#' @param N - total number of samples (denominator)
#'
#' @return: a list with 2 values: "lower" = lower CI, "upper" = upper CI
calculate_ci <- function(n, N){
  
  lower = (1/(1+(qnorm(0.975)^2)/N))*((n/N)+(qnorm(0.975)^2)/(2*N)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/N)))*sqrt((n/N)*(1-n/N)/N + (qnorm(0.975)^2)/(4*(N^2)))
  upper = (1/(1+(qnorm(0.975)^2)/N))*((n/N)+(qnorm(0.975)^2)/(2*N)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/N)))*sqrt((n/N)*(1-n/N)/N + (qnorm(0.975)^2)/(4*(N^2)))
  ls <- list(lower, upper)
  names(ls) <- c("lower", "upper")
  return(ls)
}


# Recode gender
recode_gender <- function(data) {
  data |> mutate(gender = case_when(gender == 1 ~ "Men", gender == 2 ~ "Women"))
}

# Add table theme (using flextable - need to make sure it is loaded)
add_table_theme <- function(table) {
  table |> 
    fontsize(part = "all", size = 9) |>
    bold(part = "header") |>
    font(part = "all", fontname =  "Calibri")
}

# Plot prediction data
plot_preds_data <- function(data, ylab){
  ggplot() +
    geom_line(aes(x = age, y = 100*predicted, color = index_labels), 
              linewidth = linesize, data = data) +
    geom_ribbon(aes(x = age, ymin = 100*conf.low, ymax = 100*conf.high, fill = index_labels), 
                alpha = .2, data = data) +
    scale_fill_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
    scale_color_manual(values = c(colours$`Rectal Bleeding`, colours$CIBH, colours$Neither), legend_title) +
    legend +
    facet_grid(gender~cancer_site)+
    theme_light() +
    theme +
    theme(text=element_text(size=7)) +
    labs(y = ylab,
         x = x_label)
}
