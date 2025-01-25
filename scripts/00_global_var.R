## Purpose: Common global variables ##
## Author:  Nadine Zakkak ##

# Symptom and Cancer sites codes & descriptions ----
symptom_codes <- c("CIBH" = 4, "Rectal Bleeding" = 14)
cancer_codes <- c("Colon" = 11, "Rectal" = 12)

# Figures ----
library(ggplot2)
colours <- list("#54A0C9", "#62A02C", "black")
names(colours) <- c("Rectal Bleeding", "CIBH", "Neither")
linesize <- .5
theme <-  theme(strip.background = element_rect(fill = "#F5F5F5"),
                strip.text = element_text(colour = "black"),
                text = element_text(size = 12))
legend <-   guides(color = guide_legend(override.aes=list(size = 2)),
                   fill = "none")
legend_title <- "Symptom"
x_label <- "Age (years)"
y_label <- "Probability (%)"
