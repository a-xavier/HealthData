#!/usr/bin/env Rscript
library(rcompanion) # to get pseudoR2 calculations

args <- commandArgs(trailingOnly = TRUE)
code <- args[1]
sex_flag <- args[2]
cohort_code <- args[3]
final_df_file <- args[4]

print(paste("Processing code:",code, "with sex flag:", sex_flag, "for cohort:",
cohort_code, "with final df file:", final_df_file))

df = read.csv(final_df_file)

if (sex_flag == "NA") {
    "Do nothing for NA sex flag"
} else if (sex_flag == "M") {
    df = df[df$is_male == 1, ]
} else if (sex_flag == "F") {
    df = df[df$is_male == 0, ]
}

# create phenotype column grepping code from column all_levels and 1 when found 0 otherwise 
df$phenotype = ifelse(grepl(code, df$all_levels), 1, 0)

# Full model

if (cohort_code == "GS") {
    full_formula = as.formula(
        "phenotype ~ is_male + age_at_recruitment + BMI + smokever + dpw_moderate_activity + dpw_vigorous_activity + drinking_events_per_month + PC1 + PC2 + PC3 + PC4 + PC5 + famid"
    )
} else{
     full_formula = as.formula(
        "phenotype ~ is_male + age_at_recruitment + BMI + smokever + dpw_moderate_activity + dpw_vigorous_activity + drinking_events_per_month + PC1 + PC2 + PC3 + PC4 + PC5"
    )
}

# All covariates
all_covariates_to_test = c("is_male", "age_at_recruitment", "BMI", "smokever", "dpw_moderate_activity", "dpw_vigorous_activity", "drinking_events_per_month", "PC1", "PC2", "PC3", "PC4", "PC5")
all_covariates_model = glm(full_formula, data = df, family = binomial)
all_covariates_coefficients = summary(all_covariates_model)$coefficients

# print(all_covariates_coefficients)

result_holder = list()

for (covariate in all_covariates_to_test) {

    # formula is full_formula minus the current covariate
    reduced_formula = as.formula(
        paste("phenotype ~", paste(setdiff(all_covariates_to_test, covariate), collapse = " + "))
    )
    reduced_model = glm(reduced_formula, data = df, family = binomial)

    # Get the R2 (all 3 measured) 
    r2_results = nagelkerke(all_covariates_model, null = reduced_model)$Pseudo.R.squared.for.model.vs.null

    # get results for the current covariate
    # print(paste("R2 for covariate", covariate, ":"))
    # print(r2_results)
    # print(r2_results["McFadden", "Pseudo.R.squared"])

    # Need to do for each covariate McFadden Cox and Snell (ML) and Nagelkerke (Cragg and Uhler) Pseudo.R.squared
    # plus also the full model Estimate   Std. Error      z value      Pr(>|z|)
    covar_row_result = data.frame(
        McFadden_R2 = r2_results["McFadden", "Pseudo.R.squared"],
        Cox_Snell_R2 = r2_results["Cox and Snell (ML)", "Pseudo.R.squared"],
        Nagelkerke_R2 = r2_results["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"],
        Estimate = all_covariates_coefficients[covariate, "Estimate"],
        Std_Error = all_covariates_coefficients[covariate, "Std. Error"],
        z_value = all_covariates_coefficients[covariate, "z value"],
        Pr_z = all_covariates_coefficients[covariate, "Pr(>|z|)"]
    )

    # add covar_prefix
    result_holder[[covariate]] = covar_row_result
}

# add into one row (cbind)
final_result = do.call(cbind, result_holder)

# add the code
final_result$icd10_code = code



# Output
final_output_name = paste0("processed_", code, "_", cohort_code, ".csv")

write.csv(final_result, final_output_name, row.names = FALSE)