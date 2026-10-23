#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(fastglm)
})


# ==================================== #
#          ARGUMENTS                    #
# ==================================== #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {
    stop(
        paste(
            "Expected 4 arguments:",
            "<pair_batch_file>",
            "<cohort_code>",
            "<final_df>",
            "<output_file>"
        ),
        call. = FALSE
    )
}

pair_batch_file <- args[[1]]
cohort_code <- args[[2]]
final_df <- args[[3]]
output_file <- args[[4]]

if (!file.exists(pair_batch_file)) {
    stop(
        paste("Pair batch file does not exist:", pair_batch_file),
        call. = FALSE
    )
}

if (!file.exists(final_df)) {
    stop(
        paste("Input dataframe does not exist:", final_df),
        call. = FALSE
    )
}

message(
    "Processing pair batch: ", pair_batch_file,
    " for cohort ", cohort_code,
    " using ", final_df
)


# ==================================== #
#          READ PAIR BATCH              #
# ==================================== #

pairs <- fread(
    pair_batch_file,
    header = FALSE,
    sep = "\t",
    col.names = c(
        "code_1",
        "code_2",
        "observed_both_count"
    ),
    colClasses = c(
        code_1 = "character",
        code_2 = "character",
        observed_both_count = "numeric"
    ),
    data.table = TRUE,
    showProgress = FALSE,
    nThread = 1L
)

if (nrow(pairs) == 0L) {
    stop(
        paste("The pair batch is empty:", pair_batch_file),
        call. = FALSE
    )
}

if (anyNA(pairs$code_1) || anyNA(pairs$code_2)) {
    stop(
        paste("The pair batch contains missing code values:", pair_batch_file),
        call. = FALSE
    )
}

message("Pairs in batch: ", nrow(pairs))


# ==================================== #
#          COVARIATES                   #
# ==================================== #

base_covariates <- c(
    "is_male",
    "age_at_recruitment",
    "BMI",
    "smokever",
    "dpw_moderate_activity",
    "dpw_vigorous_activity",
    "drinking_events_per_month",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5"
)

covariates <- if (cohort_code == "GS") {
    c(base_covariates, "famid")
} else {
    base_covariates
}

required_columns <- unique(c(
    "all_levels",
    covariates
))


# ==================================== #
#          LOAD COHORT DATA             #
# ==================================== #

load_start <- proc.time()[["elapsed"]]

df <- fread(
    final_df,
    select = required_columns,
    data.table = FALSE,
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

load_elapsed <- proc.time()[["elapsed"]] - load_start

message(
    "Loaded ", nrow(df),
    " rows in ",
    round(load_elapsed, 3),
    " seconds"
)

missing_columns <- setdiff(required_columns, names(df))

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "Required columns are missing:",
            paste(missing_columns, collapse = ", ")
        ),
        call. = FALSE
    )
}


# ==================================== #
#           CODE STRINGS                #
# ==================================== #

# A missing all_levels value is treated as having no recorded codes.
all_levels <- df$all_levels
all_levels[is.na(all_levels)] <- ""

participant_count <- length(all_levels)


# ==================================== #
#       COMPLETE-CASE POPULATION        #
# ==================================== #

complete_case_start <- proc.time()[["elapsed"]]

complete_rows <- complete.cases(
    df[, covariates, drop = FALSE]
)

glm_n <- sum(complete_rows)
glm_excluded_n <- length(complete_rows) - glm_n

complete_case_elapsed <-
    proc.time()[["elapsed"]] - complete_case_start

message(
    "GLM complete cases: ", glm_n,
    "; excluded due to missing covariates: ", glm_excluded_n
)

if (glm_n == 0L) {
    stop(
        "No complete observations remain for linear modelling",
        call. = FALSE
    )
}


# ==================================== #
#       COVARIATE MODEL MATRIX          #
# ==================================== #

matrix_start <- proc.time()[["elapsed"]]

covariate_formula <- reformulate(covariates)

covariate_matrix <- model.matrix(
    covariate_formula,
    data = df[complete_rows, , drop = FALSE]
)

matrix_elapsed <- proc.time()[["elapsed"]] - matrix_start

message(
    "Covariate matrix dimensions: ",
    nrow(covariate_matrix),
    " x ",
    ncol(covariate_matrix)
)

message(
    "Covariate matrix size: ",
    round(
        as.numeric(object.size(covariate_matrix)) / 1024^2,
        2
    ),
    " MB"
)

message(
    "Covariate matrix constructed in ",
    round(matrix_elapsed, 3),
    " seconds"
)


# ==================================== #
#       CODE MEMBERSHIP CACHE           #
# ==================================== #

# Membership vectors are cached within this batch.
#
# If the same code occurs in several pairs, grepl() is performed only
# once for that code during this R invocation.
code_membership_cache <- new.env(
    hash = TRUE,
    parent = emptyenv()
)

get_code_membership <- function(code) {

    cache_key <- paste0("code_", code)

    if (!exists(
        cache_key,
        envir = code_membership_cache,
        inherits = FALSE
    )) {
        membership <- grepl(
            pattern = code,
            x = all_levels,
            fixed = TRUE
        )

        assign(
            cache_key,
            membership,
            envir = code_membership_cache
        )
    }

    get(
        cache_key,
        envir = code_membership_cache,
        inherits = FALSE
    )
}


# ==================================== #
#       GENERAL HELPERS                 #
# ==================================== #

safe_divide <- function(numerator, denominator) {

    if (
        length(numerator) == 0L ||
        length(denominator) == 0L ||
        is.na(numerator) ||
        is.na(denominator) ||
        denominator == 0
    ) {
        return(NA_real_)
    }

    numerator / denominator
}


p_value_floor <- 5e-308
maximum_finite_value <- 1e308


normalise_p_value <- function(value) {

    if (
        length(value) == 0L ||
        is.na(value) ||
        is.nan(value)
    ) {
        return(NA_real_)
    }

    if (value == 0) {
        return(p_value_floor)
    }

    as.numeric(value)
}


normalise_ratio <- function(value) {

    if (
        length(value) == 0L ||
        is.na(value) ||
        is.nan(value)
    ) {
        return(NA_real_)
    }

    if (is.infinite(value)) {
        if (value > 0) {
            return(maximum_finite_value)
        }

        return(-maximum_finite_value)
    }

    if (value == 0) {
        return(p_value_floor)
    }

    as.numeric(value)
}


empty_glm_statistics <- function() {

    c(
        Estimate = NA_real_,
        `Std. Error` = NA_real_,
        `z value` = NA_real_,
        `Pr(>|z|)` = NA_real_
    )
}


combine_warning_messages <- function(...) {

    messages <- unlist(list(...), use.names = FALSE)

    messages <- messages[
        !is.na(messages) &
        nzchar(messages)
    ]

    if (length(messages) == 0L) {
        return(NA_character_)
    }

    paste(unique(messages), collapse = " | ")
}


# ==================================== #
#          FAST GLM HELPER              #
# ==================================== #

fit_target_glm <- function(
    design_matrix,
    response,
    target_column,
    model_label
) {
    result <- list(
        statistics = empty_glm_statistics(),
        converged = FALSE,
        iterations = NA_integer_,
        warning = NA_character_
    )

    if (length(unique(response)) < 2L) {
        result$warning <- paste(
            model_label,
            "was skipped because its response has no variation"
        )

        return(result)
    }

    target_values <- design_matrix[, target_column]

    if (length(unique(target_values)) < 2L) {
        result$warning <- paste(
            model_label,
            "was skipped because its target predictor has no variation"
        )

        return(result)
    }

    captured_warnings <- character()

    fit <- tryCatch(
        withCallingHandlers(
            fastglm::fastglm(
                x = design_matrix,
                y = response,
                family = binomial(),
                method = 0L,
                control = fastglm::fastglm_control(
                    tol = 1e-8,
                    maxit = 25L
                )
            ),
            warning = function(w) {
                captured_warnings <<- c(
                    captured_warnings,
                    conditionMessage(w)
                )

                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) {
            result$warning <<- paste(
                model_label,
                "failed:",
                conditionMessage(e)
            )

            NULL
        }
    )

    if (is.null(fit)) {
        return(result)
    }

    if (!is.null(fit$converged)) {
        result$converged <- isTRUE(fit$converged)
    }

    if (!is.null(fit$iter)) {
        result$iterations <- as.integer(fit$iter)
    } else if (!is.null(fit$iterations)) {
        result$iterations <- as.integer(fit$iterations)
    }

    if (length(captured_warnings) > 0L) {
        result$warning <- combine_warning_messages(
            result$warning,
            captured_warnings
        )
    }

    coefficient_table <- tryCatch(
        summary(fit)$coefficients,
        error = function(e) {
            result$warning <<- combine_warning_messages(
                result$warning,
                paste(
                    "Could not calculate coefficient statistics:",
                    conditionMessage(e)
                )
            )

            NULL
        }
    )

    if (is.null(coefficient_table)) {
        return(result)
    }

    if (!target_column %in% rownames(coefficient_table)) {
        result$warning <- combine_warning_messages(
            result$warning,
            paste(
                "Coefficient",
                target_column,
                "was not estimable in",
                model_label
            )
        )

        return(result)
    }

    result$statistics <- coefficient_table[
        target_column,
        c(
            "Estimate",
            "Std. Error",
            "z value",
            "Pr(>|z|)"
        )
    ]

    result
}


# ==================================== #
#       PROCESS ONE CODE PAIR           #
# ==================================== #

process_pair <- function(
    code_1,
    code_2,
    observed_both_count
) {
    pair_start <- proc.time()[["elapsed"]]

    # -------------------------------- #
    # Code membership                   #
    # -------------------------------- #

    matching_start <- proc.time()[["elapsed"]]

    code_1_present <- get_code_membership(code_1)
    code_2_present <- get_code_membership(code_2)

    code_1_count <- sum(code_1_present)
    code_2_count <- sum(code_2_present)
    both_count <- sum(code_1_present & code_2_present)

    cell_c <- as.numeric(both_count)
    cell_a <- as.numeric(code_1_count - both_count)
    cell_b <- as.numeric(code_2_count - both_count)

    cell_d <- as.numeric(
        participant_count -
        cell_a -
        cell_b -
        cell_c
    )

    matching_elapsed <-
        proc.time()[["elapsed"]] - matching_start


    # -------------------------------- #
    # Fisher exact test                 #
    # -------------------------------- #

    fisher_matrix <- matrix(
        c(
            cell_c,
            cell_a,
            cell_b,
            cell_d
        ),
        nrow = 2L
    )

    fisher_warning <- NA_character_

    fisher_result <- tryCatch(
        withCallingHandlers(
            fisher.test(fisher_matrix),
            warning = function(w) {
                fisher_warning <<- combine_warning_messages(
                    fisher_warning,
                    conditionMessage(w)
                )

                invokeRestart("muffleWarning")
            }
        ),
        error = function(e) {
            fisher_warning <<- combine_warning_messages(
                fisher_warning,
                paste(
                    "Fisher test failed:",
                    conditionMessage(e)
                )
            )

            NULL
        }
    )

    if (is.null(fisher_result)) {
        odds_ratio <- NA_real_
        fisher_p_value <- NA_real_
        fisher_conf_int <- c(
            NA_real_,
            NA_real_
        )
    } else {
        odds_ratio <- unname(fisher_result$estimate)
        fisher_p_value <- fisher_result$p.value
        fisher_conf_int <- unname(fisher_result$conf.int)
    }


    # -------------------------------- #
    # Relative risks                    #
    # -------------------------------- #

    risk_code_1_with_code_2 <- safe_divide(
        cell_c,
        cell_c + cell_b
    )

    risk_code_1_without_code_2 <- safe_divide(
        cell_a,
        cell_a + cell_d
    )

    risk_code_2_with_code_1 <- safe_divide(
        cell_c,
        cell_c + cell_a
    )

    risk_code_2_without_code_1 <- safe_divide(
        cell_b,
        cell_b + cell_d
    )

    rr_code1_given_code2 <- safe_divide(
        risk_code_1_with_code_2,
        risk_code_1_without_code_2
    )

    rr_code2_given_code1 <- safe_divide(
        risk_code_2_with_code_1,
        risk_code_2_without_code_1
    )


    # -------------------------------- #
    # Phi coefficient                   #
    # -------------------------------- #

    phi_denominator <- (
        (cell_a + cell_b) *
        (cell_c + cell_d) *
        (cell_a + cell_c) *
        (cell_b + cell_d)
    )

    if (
        !is.finite(phi_denominator) ||
        phi_denominator <= 0
    ) {
        phi_coefficient <- NA_real_
    } else {
        phi_coefficient <- (
            (cell_c * cell_d) -
            (cell_a * cell_b)
        ) / sqrt(phi_denominator)
    }


    # -------------------------------- #
    # Jaccard index                     #
    # -------------------------------- #

    jaccard_denominator <- (
        cell_a +
        cell_b +
        cell_c
    )

    if (jaccard_denominator == 0) {
        jaccard_index <- NA_real_
    } else {
        jaccard_index <- (
            cell_c /
            jaccard_denominator
        )
    }


    # -------------------------------- #
    # Logistic regression               #
    # -------------------------------- #

    glm_start <- proc.time()[["elapsed"]]

    glm_code_1_response <- as.integer(
        code_1_present[complete_rows]
    )

    glm_code_2_response <- as.integer(
        code_2_present[complete_rows]
    )

    # Construct one design matrix and reuse its final column.
    design_matrix <- cbind(
        covariate_matrix,
        code_2 = glm_code_2_response
    )

    glm_code_1_result <- fit_target_glm(
        design_matrix = design_matrix,
        response = glm_code_1_response,
        target_column = "code_2",
        model_label = "code_1 ~ code_2 + covariates"
    )

    # Reuse the same allocated design matrix for the reverse model.
    design_matrix[, ncol(design_matrix)] <-
        glm_code_1_response

    colnames(design_matrix)[ncol(design_matrix)] <-
        "code_1"

    glm_code_2_result <- fit_target_glm(
        design_matrix = design_matrix,
        response = glm_code_2_response,
        target_column = "code_1",
        model_label = "code_2 ~ code_1 + covariates"
    )

    rm(design_matrix)

    glm_elapsed <- proc.time()[["elapsed"]] - glm_start

    code_1_vs_code_2 <- glm_code_1_result$statistics
    code_2_vs_code_1 <- glm_code_2_result$statistics


    # -------------------------------- #
    # Numerical normalisation           #
    # -------------------------------- #

    fisher_p_value <- normalise_p_value(
        fisher_p_value
    )

    odds_ratio <- normalise_ratio(
        odds_ratio
    )

    rr_code1_given_code2 <- normalise_ratio(
        rr_code1_given_code2
    )

    rr_code2_given_code1 <- normalise_ratio(
        rr_code2_given_code1
    )

    fisher_conf_int[[1]] <- normalise_ratio(
        fisher_conf_int[[1]]
    )

    fisher_conf_int[[2]] <- normalise_ratio(
        fisher_conf_int[[2]]
    )

    code_1_vs_code_2["Pr(>|z|)"] <-
        normalise_p_value(
            code_1_vs_code_2["Pr(>|z|)"]
        )

    code_2_vs_code_1["Pr(>|z|)"] <-
        normalise_p_value(
            code_2_vs_code_1["Pr(>|z|)"]
        )


    # -------------------------------- #
    # Result row                        #
    # -------------------------------- #

    pair_elapsed <- proc.time()[["elapsed"]] - pair_start

    data.table(
        code_1 = code_1,
        code_2 = code_2,

        input_observed_both_count =
            observed_both_count,

        code_1_only = cell_a,
        code_2_only = cell_b,
        both_codes = cell_c,
        neither_code = cell_d,

        fisher_odds_ratio = odds_ratio,
        fisher_p_value = fisher_p_value,
        fisher_conf_int_low =
            fisher_conf_int[[1]],
        fisher_conf_int_high =
            fisher_conf_int[[2]],
        fisher_warning = fisher_warning,

        rr_code1_given_code2 =
            rr_code1_given_code2,
        rr_code2_given_code1 =
            rr_code2_given_code1,

        phi_coefficient = phi_coefficient,
        jaccard_index = jaccard_index,

        glm_sample_size = glm_n,
        glm_excluded_missing_covariates =
            glm_excluded_n,

        code_1_vs_code_2_estimate =
            unname(
                code_1_vs_code_2["Estimate"]
            ),

        code_1_vs_code_2_std_error =
            unname(
                code_1_vs_code_2["Std. Error"]
            ),

        code_1_vs_code_2_z_value =
            unname(
                code_1_vs_code_2["z value"]
            ),

        code_1_vs_code_2_p_value =
            unname(
                code_1_vs_code_2["Pr(>|z|)"]
            ),

        code_1_vs_code_2_converged =
            glm_code_1_result$converged,

        code_1_vs_code_2_iterations =
            glm_code_1_result$iterations,

        code_1_vs_code_2_warning =
            glm_code_1_result$warning,

        code_2_vs_code_1_estimate =
            unname(
                code_2_vs_code_1["Estimate"]
            ),

        code_2_vs_code_1_std_error =
            unname(
                code_2_vs_code_1["Std. Error"]
            ),

        code_2_vs_code_1_z_value =
            unname(
                code_2_vs_code_1["z value"]
            ),

        code_2_vs_code_1_p_value =
            unname(
                code_2_vs_code_1["Pr(>|z|)"]
            ),

        code_2_vs_code_1_converged =
            glm_code_2_result$converged,

        code_2_vs_code_1_iterations =
            glm_code_2_result$iterations,

        code_2_vs_code_1_warning =
            glm_code_2_result$warning,

        code_matching_seconds =
            matching_elapsed,

        glm_seconds = glm_elapsed,

        pair_total_seconds = pair_elapsed
    )
}


# ==================================== #
#       PROCESS COMPLETE BATCH          #
# ==================================== #

batch_start <- proc.time()[["elapsed"]]

results <- vector(
    mode = "list",
    length = nrow(pairs)
)

for (pair_index in seq_len(nrow(pairs))) {

    code_1 <- pairs$code_1[[pair_index]]
    code_2 <- pairs$code_2[[pair_index]]

    observed_both_count <-
        pairs$observed_both_count[[pair_index]]

    message(
        "Processing pair ",
        pair_index,
        " of ",
        nrow(pairs),
        ": ",
        code_1,
        " / ",
        code_2
    )

    results[[pair_index]] <- process_pair(
        code_1 = code_1,
        code_2 = code_2,
        observed_both_count = observed_both_count
    )
}

report_df <- rbindlist(
    results,
    use.names = TRUE,
    fill = TRUE
)

batch_elapsed <- proc.time()[["elapsed"]] - batch_start

report_df[, data_load_seconds := load_elapsed]
report_df[, complete_case_seconds := complete_case_elapsed]
report_df[, model_matrix_seconds := matrix_elapsed]
report_df[, batch_total_seconds := batch_elapsed]
report_df[, batch_pair_count := nrow(pairs)]
report_df[, cohort_code := cohort_code]


# ==================================== #
#             OUTPUT                    #
# ==================================== #

fwrite(
    report_df,
    file = output_file,
    sep = ",",
    quote = TRUE,
    na = "NA",
    nThread = 1L
)

message(
    "Processed ", nrow(report_df),
    " pairs in ",
    round(batch_elapsed, 3),
    " seconds"
)

message(
    "Cached membership vectors: ",
    length(ls(
        envir = code_membership_cache,
        all.names = TRUE
    ))
)

message("Wrote output: ", output_file)