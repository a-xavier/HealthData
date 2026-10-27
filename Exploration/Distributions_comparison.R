library(tidyverse)
library(data.table)
gs_df = fread("/Users/alex/Documents/ICD10_Project_Official/Generation_Scotland/Final_Data_GS.csv")
ukb_df = fread("/Users/alex/Documents/ICD10_Project_Official/UKB/Final_Data_UKB.csv")

# Age 


# Continuous is PC1 - 10 - dpw_moderate_activity - dpw_vigorous_activity - age_at_recruitment
# drinking_events_per_month - BMI

continuous_covariates = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                          "dpw_moderate_activity", "dpw_vigorous_activity",
                          "age_at_recruitment", "drinking_events_per_month", "BMI")

quantile(ukb_df$age_at_recruitment, na.rm = T)

for (covar in continuous_covariates) {
  gs_covar = gs_df[, c("eid", covar), with = FALSE]
  gs_covar$cohort = "GS"
  ukb_covar = ukb_df[, c("eid", covar), with = FALSE]
  ukb_covar$cohort = "UKB"
  
  merged_covar = rbind(gs_covar, ukb_covar)
  
  # Density plot with fill is cohort 
  p <- ggplot(data = merged_covar, aes_string(x = covar, fill = "cohort")) +
    geom_density(alpha = 0.5) + 
    theme_minimal() + 
    ggtitle(paste("Density plot of", covar)) +
    xlab(covar) +
    ylab("Density")
  print(p)
}



# Cat is is_male and smokever
cat_covars = c("is_male", "smokever")
for (covar in cat_covars) {
  gs_covar = gs_df[, c("eid", covar), with = FALSE]
  gs_covar$cohort = "GS"
  ukb_covar = ukb_df[, c("eid", covar), with = FALSE]
  ukb_covar$cohort = "UKB"
  
  merged_covar = rbind(gs_covar, ukb_covar)
  
  # Bar plot with fill is cohort 
  p <- ggplot(data = merged_covar, aes_string(x = covar, fill = "cohort")) +
    geom_bar(position = "dodge") + 
    theme_minimal() + 
    ggtitle(paste("Bar plot of", covar)) +
    xlab(covar) +
    ylab("Count")
  print(p)
}
