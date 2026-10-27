#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})


# ==================================== #
#              ARGUMENTS                #
# ==================================== #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2L) {
    stop(
        paste(
            "Expected 2 arguments:",
            "<cohort_code>",
            "<output_file>"
        ),
        call. = FALSE
    )
}

cohort_code <- trimws(args[[1]])
output_file <- args[[2]]

if (
    is.na(cohort_code) ||
    cohort_code == ""
) {
    stop(
        "cohort_code cannot be empty",
        call. = FALSE
    )
}

message("Merging temporal statistics for cohort: ", cohort_code)


# ==================================== #
#          HELPER FUNCTIONS             #
# ==================================== #

format_count <- function(value) {
    format(
        value,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
    )
}


minimum_or_na <- function(values) {

    values <- values[
        !is.na(values)
    ]

    if (length(values) == 0L) {
        return(NA_real_)
    }

    as.numeric(
        min(values)
    )
}


maximum_or_na <- function(values) {

    values <- values[
        !is.na(values)
    ]

    if (length(values) == 0L) {
        return(NA_real_)
    }

    as.numeric(
        max(values)
    )
}


quantile_or_na <- function(
    values,
    probability
) {
    values <- values[
        !is.na(values)
    ]

    if (length(values) == 0L) {
        return(NA_real_)
    }

    as.numeric(
        quantile(
            values,
            probs = probability,
            names = FALSE,
            type = 7L
        )
    )
}


mean_from_summary <- function(
    count,
    sum_value
) {
    if (
        is.na(count) ||
        count <= 0L
    ) {
        return(NA_real_)
    }

    sum_value / count
}


sd_from_summary <- function(
    count,
    sum_value,
    sum_squared
) {
    if (
        is.na(count) ||
        count <= 1L
    ) {
        return(NA_real_)
    }

    variance <- (
        sum_squared -
        (
            sum_value^2 /
            count
        )
    ) / (
        count - 1L
    )

    # Protect against tiny negative values introduced by
    # floating-point rounding.
    variance <- max(
        variance,
        0
    )

    sqrt(variance)
}


exact_binomial_interval <- function(
    successes,
    total,
    confidence_level = 0.95
) {
    if (
        is.na(successes) ||
        is.na(total) ||
        total <= 0L ||
        successes < 0L ||
        successes > total
    ) {
        return(
            list(
                estimate = NA_real_,
                lower = NA_real_,
                upper = NA_real_
            )
        )
    }

    alpha <- 1 - confidence_level
    failures <- total - successes

    estimate <- successes / total

    lower <- if (successes == 0L) {
        0
    } else {
        qbeta(
            alpha / 2,
            successes,
            failures + 1L
        )
    }

    upper <- if (failures == 0L) {
        1
    } else {
        qbeta(
            1 - alpha / 2,
            successes + 1L,
            failures
        )
    }

    list(
        estimate = estimate,
        lower = lower,
        upper = upper
    )
}


exact_direction_p_value <- function(
    code_1_first,
    code_2_first
) {
    total <- (
        code_1_first +
        code_2_first
    )

    if (
        is.na(total) ||
        total <= 0L
    ) {
        return(NA_real_)
    }

    smaller_count <- min(
        code_1_first,
        code_2_first
    )

    min(
        1,
        2 * pbinom(
            smaller_count,
            size = total,
            prob = 0.5
        )
    )
}


# ==================================== #
#          DISCOVER INPUT FILES         #
# ==================================== #

partial_stats_files <- sort(
    Sys.glob(
        "partial_stats_*.csv"
    )
)

partial_delay_files <- sort(
    Sys.glob(
        "partial_delays_*.csv"
    )
)

if (length(partial_stats_files) == 0L) {
    stop(
        "No partial_stats_*.csv files were found",
        call. = FALSE
    )
}

if (length(partial_delay_files) == 0L) {
    stop(
        "No partial_delays_*.csv files were found",
        call. = FALSE
    )
}

message(
    "Partial statistics files: ",
    format_count(
        length(partial_stats_files)
    )
)

message(
    "Partial delay files: ",
    format_count(
        length(partial_delay_files)
    )
)


# ==================================== #
#       LOAD PARTIAL STATISTICS         #
# ==================================== #

partial_stats_list <- lapply(
    partial_stats_files,
    function(input_file) {
        fread(
            input_file,
            sep = ",",
            na.strings = c(
                "NA",
                ""
            ),
            showProgress = FALSE,
            nThread = 1L,
            check.names = FALSE
        )
    }
)

partial_stats <- rbindlist(
    partial_stats_list,
    use.names = TRUE,
    fill = FALSE
)

rm(partial_stats_list)

message(
    "Partial statistics rows loaded: ",
    format_count(
        nrow(partial_stats)
    )
)

if (nrow(partial_stats) == 0L) {
    stop(
        "The combined partial statistics table is empty",
        call. = FALSE
    )
}


# ==================================== #
#       VALIDATE PARTIAL STATISTICS     #
# ==================================== #

required_stats_columns <- c(
    "cohort_code",
    "code_1",
    "code_2",
    "participants_with_both",
    "code_1_before_code_2_count",
    "code_2_before_code_1_count",
    "same_date_count",
    "temporally_ordered_count",

    "code_1_to_code_2_sum_days",
    "code_1_to_code_2_sum_squared_days",
    "code_1_to_code_2_min_days",
    "code_1_to_code_2_max_days",

    "code_2_to_code_1_sum_days",
    "code_2_to_code_1_sum_squared_days",
    "code_2_to_code_1_min_days",
    "code_2_to_code_1_max_days"
)

missing_stats_columns <- setdiff(
    required_stats_columns,
    names(partial_stats)
)

if (length(missing_stats_columns) > 0L) {
    stop(
        paste(
            "Partial statistics are missing required columns:",
            paste(
                missing_stats_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

partial_stats[
    ,
    cohort_code := trimws(cohort_code)
]

partial_stats[
    ,
    code_1 := toupper(
        trimws(code_1)
    )
]

partial_stats[
    ,
    code_2 := toupper(
        trimws(code_2)
    )
]

input_cohorts <- unique(
    partial_stats$cohort_code
)

if (
    length(input_cohorts) != 1L ||
    is.na(input_cohorts[[1]]) ||
    input_cohorts[[1]] != cohort_code
) {
    stop(
        paste(
            "Partial statistics contain a cohort that does not match",
            cohort_code
        ),
        call. = FALSE
    )
}


# ==================================== #
#      IDENTIFY MERGEABLE COLUMNS       #
# ==================================== #

identifier_columns <- c(
    "cohort_code",
    "code_1",
    "code_2"
)

minimum_columns <- names(partial_stats)[
    grepl(
        "(_min_days|_age_min)$",
        names(partial_stats)
    )
]

maximum_columns <- names(partial_stats)[
    grepl(
        "(_max_days|_age_max)$",
        names(partial_stats)
    )
]

numeric_columns <- names(partial_stats)[
    vapply(
        partial_stats,
        is.numeric,
        logical(1)
    )
]

sum_columns <- setdiff(
    numeric_columns,
    c(
        minimum_columns,
        maximum_columns
    )
)

message(
    "Summed columns: ",
    length(sum_columns)
)

message(
    "Minimum columns: ",
    length(minimum_columns)
)

message(
    "Maximum columns: ",
    length(maximum_columns)
)


# ==================================== #
#      MERGE SUMMABLE STATISTICS        #
# ==================================== #

merged_sums <- partial_stats[
    ,
    lapply(
        .SD,
        sum,
        na.rm = TRUE
    ),
    by = .(
        cohort_code,
        code_1,
        code_2
    ),
    .SDcols = sum_columns
]


# ==================================== #
#           MERGE MINIMUMS              #
# ==================================== #

if (length(minimum_columns) > 0L) {

    merged_minimums <- partial_stats[
        ,
        lapply(
            .SD,
            minimum_or_na
        ),
        by = .(
            cohort_code,
            code_1,
            code_2
        ),
        .SDcols = minimum_columns
    ]

    merged_stats <- merge(
        merged_sums,
        merged_minimums,
        by = c(
            "cohort_code",
            "code_1",
            "code_2"
        ),
        all = TRUE,
        sort = FALSE
    )

    rm(merged_minimums)

} else {
    merged_stats <- merged_sums
}

rm(merged_sums)


# ==================================== #
#           MERGE MAXIMUMS              #
# ==================================== #

if (length(maximum_columns) > 0L) {

    merged_maximums <- partial_stats[
        ,
        lapply(
            .SD,
            maximum_or_na
        ),
        by = .(
            cohort_code,
            code_1,
            code_2
        ),
        .SDcols = maximum_columns
    ]

    merged_stats <- merge(
        merged_stats,
        merged_maximums,
        by = c(
            "cohort_code",
            "code_1",
            "code_2"
        ),
        all = TRUE,
        sort = FALSE
    )

    rm(merged_maximums)
}

rm(partial_stats)


# ==================================== #
#         VALIDATE MERGED COUNTS        #
# ==================================== #

merged_stats[
    ,
    count_validation := (
        code_1_before_code_2_count +
        code_2_before_code_1_count +
        same_date_count
    )
]

invalid_count_rows <- merged_stats[
    count_validation !=
        participants_with_both
]

if (nrow(invalid_count_rows) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_count_rows)
            ),
            "merged pairs that failed temporal count validation"
        ),
        call. = FALSE
    )
}

merged_stats[
    ,
    count_validation := NULL
]

if (any(
    merged_stats$temporally_ordered_count !=
        (
            merged_stats$code_1_before_code_2_count +
            merged_stats$code_2_before_code_1_count
        )
)) {
    stop(
        paste(
            "temporally_ordered_count failed validation",
            "after merging"
        ),
        call. = FALSE
    )
}


# ==================================== #
#         LOAD PARTIAL DELAYS           #
# ==================================== #

partial_delay_list <- lapply(
    partial_delay_files,
    function(input_file) {
        fread(
            input_file,
            sep = ",",
            na.strings = c(
                "NA",
                ""
            ),
            colClasses = list(
                character = c(
                    "cohort_code",
                    "code_1",
                    "code_2",
                    "temporal_class"
                ),
                integer = c(
                    "days_between",
                    "is_male"
                ),
                numeric = c(
                    "age_from",
                    "age_to"
                )
            ),
            showProgress = FALSE,
            nThread = 1L,
            check.names = FALSE
        )
    }
)

partial_delays <- rbindlist(
    partial_delay_list,
    use.names = TRUE,
    fill = FALSE
)

rm(partial_delay_list)

message(
    "Directed delay rows loaded: ",
    format_count(
        nrow(partial_delays)
    )
)

if (nrow(partial_delays) == 0L) {
    stop(
        "The combined directed-delay table is empty",
        call. = FALSE
    )
}


# ==================================== #
#        VALIDATE DELAY TABLE           #
# ==================================== #

required_delay_columns <- c(
    "cohort_code",
    "code_1",
    "code_2",
    "temporal_class",
    "days_between",
    "age_from",
    "age_to",
    "is_male"
)

missing_delay_columns <- setdiff(
    required_delay_columns,
    names(partial_delays)
)

if (length(missing_delay_columns) > 0L) {
    stop(
        paste(
            "Partial delay files are missing required columns:",
            paste(
                missing_delay_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

partial_delays[
    ,
    cohort_code := trimws(cohort_code)
]

partial_delays[
    ,
    code_1 := toupper(
        trimws(code_1)
    )
]

partial_delays[
    ,
    code_2 := toupper(
        trimws(code_2)
    )
]

valid_temporal_classes <- c(
    "code_1_to_code_2",
    "code_2_to_code_1"
)

if (any(
    !partial_delays$temporal_class %in%
        valid_temporal_classes
)) {
    stop(
        "Unknown temporal_class values were found",
        call. = FALSE
    )
}

if (any(
    partial_delays$days_between <= 0L
)) {
    stop(
        "Directed delay rows contain non-positive delays",
        call. = FALSE
    )
}


# ==================================== #
#       EXACT DELAY QUANTILES           #
# ==================================== #

message("Calculating exact delay and age quantiles...")


# ------------------------------------ #
#      Code 1 to code 2                 #
# ------------------------------------ #

quantiles_1_to_2 <- partial_delays[
    temporal_class == "code_1_to_code_2",
    .(
        code_1_to_code_2_q1_days =
            quantile_or_na(
                days_between,
                0.25
            ),

        code_1_to_code_2_median_days =
            quantile_or_na(
                days_between,
                0.50
            ),

        code_1_to_code_2_q3_days =
            quantile_or_na(
                days_between,
                0.75
            ),

        code_1_to_code_2_source_age_q1 =
            quantile_or_na(
                age_from,
                0.25
            ),

        code_1_to_code_2_source_age_median =
            quantile_or_na(
                age_from,
                0.50
            ),

        code_1_to_code_2_source_age_q3 =
            quantile_or_na(
                age_from,
                0.75
            ),

        code_1_to_code_2_target_age_q1 =
            quantile_or_na(
                age_to,
                0.25
            ),

        code_1_to_code_2_target_age_median =
            quantile_or_na(
                age_to,
                0.50
            ),

        code_1_to_code_2_target_age_q3 =
            quantile_or_na(
                age_to,
                0.75
            )
    ),
    by = .(
        cohort_code,
        code_1,
        code_2
    )
]


# ------------------------------------ #
#      Code 2 to code 1                 #
# ------------------------------------ #

quantiles_2_to_1 <- partial_delays[
    temporal_class == "code_2_to_code_1",
    .(
        code_2_to_code_1_q1_days =
            quantile_or_na(
                days_between,
                0.25
            ),

        code_2_to_code_1_median_days =
            quantile_or_na(
                days_between,
                0.50
            ),

        code_2_to_code_1_q3_days =
            quantile_or_na(
                days_between,
                0.75
            ),

        code_2_to_code_1_source_age_q1 =
            quantile_or_na(
                age_from,
                0.25
            ),

        code_2_to_code_1_source_age_median =
            quantile_or_na(
                age_from,
                0.50
            ),

        code_2_to_code_1_source_age_q3 =
            quantile_or_na(
                age_from,
                0.75
            ),

        code_2_to_code_1_target_age_q1 =
            quantile_or_na(
                age_to,
                0.25
            ),

        code_2_to_code_1_target_age_median =
            quantile_or_na(
                age_to,
                0.50
            ),

        code_2_to_code_1_target_age_q3 =
            quantile_or_na(
                age_to,
                0.75
            )
    ),
    by = .(
        cohort_code,
        code_1,
        code_2
    )
]

rm(partial_delays)


# ==================================== #
#          ATTACH QUANTILES             #
# ==================================== #

merged_stats <- merge(
    merged_stats,
    quantiles_1_to_2,
    by = c(
        "cohort_code",
        "code_1",
        "code_2"
    ),
    all.x = TRUE,
    sort = FALSE
)

merged_stats <- merge(
    merged_stats,
    quantiles_2_to_1,
    by = c(
        "cohort_code",
        "code_1",
        "code_2"
    ),
    all.x = TRUE,
    sort = FALSE
)

rm(
    quantiles_1_to_2,
    quantiles_2_to_1
)


# ==================================== #
#       DERIVE TIME STATISTICS          #
# ==================================== #

merged_stats[
    ,
    code_1_to_code_2_mean_days := mapply(
        mean_from_summary,
        code_1_before_code_2_count,
        code_1_to_code_2_sum_days
    )
]

merged_stats[
    ,
    code_1_to_code_2_sd_days := mapply(
        sd_from_summary,
        code_1_before_code_2_count,
        code_1_to_code_2_sum_days,
        code_1_to_code_2_sum_squared_days
    )
]

merged_stats[
    ,
    code_2_to_code_1_mean_days := mapply(
        mean_from_summary,
        code_2_before_code_1_count,
        code_2_to_code_1_sum_days
    )
]

merged_stats[
    ,
    code_2_to_code_1_sd_days := mapply(
        sd_from_summary,
        code_2_before_code_1_count,
        code_2_to_code_1_sum_days,
        code_2_to_code_1_sum_squared_days
    )
]


# ==================================== #
#        DERIVE AGE STATISTICS          #
# ==================================== #

merged_stats[
    ,
    code_1_to_code_2_source_age_mean :=
        mapply(
            mean_from_summary,
            code_1_to_code_2_source_age_count,
            code_1_to_code_2_source_age_sum
        )
]

merged_stats[
    ,
    code_1_to_code_2_source_age_sd :=
        mapply(
            sd_from_summary,
            code_1_to_code_2_source_age_count,
            code_1_to_code_2_source_age_sum,
            code_1_to_code_2_source_age_sum_squared
        )
]

merged_stats[
    ,
    code_1_to_code_2_target_age_mean :=
        mapply(
            mean_from_summary,
            code_1_to_code_2_target_age_count,
            code_1_to_code_2_target_age_sum
        )
]

merged_stats[
    ,
    code_1_to_code_2_target_age_sd :=
        mapply(
            sd_from_summary,
            code_1_to_code_2_target_age_count,
            code_1_to_code_2_target_age_sum,
            code_1_to_code_2_target_age_sum_squared
        )
]

merged_stats[
    ,
    code_2_to_code_1_source_age_mean :=
        mapply(
            mean_from_summary,
            code_2_to_code_1_source_age_count,
            code_2_to_code_1_source_age_sum
        )
]

merged_stats[
    ,
    code_2_to_code_1_source_age_sd :=
        mapply(
            sd_from_summary,
            code_2_to_code_1_source_age_count,
            code_2_to_code_1_source_age_sum,
            code_2_to_code_1_source_age_sum_squared
        )
]

merged_stats[
    ,
    code_2_to_code_1_target_age_mean :=
        mapply(
            mean_from_summary,
            code_2_to_code_1_target_age_count,
            code_2_to_code_1_target_age_sum
        )
]

merged_stats[
    ,
    code_2_to_code_1_target_age_sd :=
        mapply(
            sd_from_summary,
            code_2_to_code_1_target_age_count,
            code_2_to_code_1_target_age_sum,
            code_2_to_code_1_target_age_sum_squared
        )
]


# ==================================== #
#       DIRECTION PROPORTIONS           #
# ==================================== #

merged_stats[
    ,
    code_1_first_among_all_with_both :=
        code_1_before_code_2_count /
        participants_with_both
]

merged_stats[
    ,
    code_2_first_among_all_with_both :=
        code_2_before_code_1_count /
        participants_with_both
]

merged_stats[
    ,
    same_date_proportion :=
        same_date_count /
        participants_with_both
]

merged_stats[
    ,
    code_1_first_among_ordered := fifelse(
        temporally_ordered_count > 0L,

        code_1_before_code_2_count /
            temporally_ordered_count,

        NA_real_
    )
]

merged_stats[
    ,
    code_2_first_among_ordered := fifelse(
        temporally_ordered_count > 0L,

        code_2_before_code_1_count /
            temporally_ordered_count,

        NA_real_
    )
]


# ==================================== #
#       DIRECTIONAL EFFECT SIZE         #
# ==================================== #

merged_stats[
    ,
    direction_balance := fifelse(
        temporally_ordered_count > 0L,

        (
            code_1_before_code_2_count -
            code_2_before_code_1_count
        ) /
            temporally_ordered_count,

        NA_real_
    )
]

merged_stats[
    ,
    direction_odds_corrected := (
        code_1_before_code_2_count +
        0.5
    ) / (
        code_2_before_code_1_count +
        0.5
    )
]

merged_stats[
    ,
    log_direction_odds := log(
        direction_odds_corrected
    )
]

merged_stats[
    ,
    dominant_recorded_direction := fcase(
        temporally_ordered_count == 0L,
        "same_date_only",

        code_1_before_code_2_count >
            code_2_before_code_1_count,
        paste0(
            code_1,
            "_before_",
            code_2
        ),

        code_2_before_code_1_count >
            code_1_before_code_2_count,
        paste0(
            code_2,
            "_before_",
            code_1
        ),

        default = "equal"
    )
]


# ==================================== #
#        EXACT DIRECTION TEST           #
# ==================================== #

direction_results <- rbindlist(
    lapply(
        seq_len(nrow(merged_stats)),
        function(row_index) {

            code_1_first <- (
                merged_stats$
                    code_1_before_code_2_count[
                        row_index
                    ]
            )

            code_2_first <- (
                merged_stats$
                    code_2_before_code_1_count[
                        row_index
                    ]
            )

            ordered_count <- (
                code_1_first +
                code_2_first
            )

            interval <- exact_binomial_interval(
                successes = code_1_first,
                total = ordered_count
            )

            data.table(
                code_1_first_ci_low =
                    interval$lower,

                code_1_first_ci_high =
                    interval$upper,

                direction_p_value =
                    exact_direction_p_value(
                        code_1_first,
                        code_2_first
                    )
            )
        }
    )
)

merged_stats <- cbind(
    merged_stats,
    direction_results
)

rm(direction_results)

merged_stats[
    ,
    direction_fdr := p.adjust(
        direction_p_value,
        method = "BH"
    )
]


# ==================================== #
#       FIXED-WINDOW PROPORTIONS        #
# ==================================== #

window_names <- c(
    "7_days",
    "30_days",
    "90_days",
    "1_year",
    "2_years",
    "5_years",
    "10_years"
)

for (window_name in window_names) {

    count_column_1_to_2 <- paste0(
        "code_1_to_code_2_within_",
        window_name
    )

    proportion_column_1_to_2 <- paste0(
        count_column_1_to_2,
        "_proportion"
    )

    merged_stats[
        ,
        (proportion_column_1_to_2) :=
            fifelse(
                code_1_before_code_2_count >
                    0L,

                get(
                    count_column_1_to_2
                ) /
                    code_1_before_code_2_count,

                NA_real_
            )
    ]

    count_column_2_to_1 <- paste0(
        "code_2_to_code_1_within_",
        window_name
    )

    proportion_column_2_to_1 <- paste0(
        count_column_2_to_1,
        "_proportion"
    )

    merged_stats[
        ,
        (proportion_column_2_to_1) :=
            fifelse(
                code_2_before_code_1_count >
                    0L,

                get(
                    count_column_2_to_1
                ) /
                    code_2_before_code_1_count,

                NA_real_
            )
    ]
}


# ==================================== #
#        SEX-SPECIFIC PROPORTIONS       #
# ==================================== #

merged_stats[
    ,
    male_ordered_count := (
        male_code_1_before_code_2 +
        male_code_2_before_code_1
    )
]

merged_stats[
    ,
    female_ordered_count := (
        female_code_1_before_code_2 +
        female_code_2_before_code_1
    )
]

merged_stats[
    ,
    male_code_1_first_among_ordered :=
        fifelse(
            male_ordered_count > 0L,

            male_code_1_before_code_2 /
                male_ordered_count,

            NA_real_
        )
]

merged_stats[
    ,
    female_code_1_first_among_ordered :=
        fifelse(
            female_ordered_count > 0L,

            female_code_1_before_code_2 /
                female_ordered_count,

            NA_real_
        )
]


# ==================================== #
#      REMOVE INTERNAL ACCUMULATORS     #
# ==================================== #

internal_columns <- c(
    "code_1_to_code_2_sum_days",
    "code_1_to_code_2_sum_squared_days",
    "code_2_to_code_1_sum_days",
    "code_2_to_code_1_sum_squared_days",

    "code_1_to_code_2_source_age_sum",
    "code_1_to_code_2_source_age_sum_squared",
    "code_1_to_code_2_target_age_sum",
    "code_1_to_code_2_target_age_sum_squared",

    "code_2_to_code_1_source_age_sum",
    "code_2_to_code_1_source_age_sum_squared",
    "code_2_to_code_1_target_age_sum",
    "code_2_to_code_1_target_age_sum_squared"
)

existing_internal_columns <- intersect(
    internal_columns,
    names(merged_stats)
)

if (length(existing_internal_columns) > 0L) {
    merged_stats[
        ,
        (existing_internal_columns) := NULL
    ]
}


# ==================================== #
#              SORT OUTPUT              #
# ==================================== #

setorder(
    merged_stats,
    -participants_with_both,
    -temporally_ordered_count,
    code_1,
    code_2
)


# ==================================== #
#          FINAL VALIDATION             #
# ==================================== #

if (any(
    merged_stats$code_1_first_among_ordered <
        0,
    na.rm = TRUE
)) {
    stop(
        "Negative direction proportions were produced",
        call. = FALSE
    )
}

if (any(
    merged_stats$code_1_first_among_ordered >
        1,
    na.rm = TRUE
)) {
    stop(
        "Direction proportions above one were produced",
        call. = FALSE
    )
}

if (any(
    merged_stats$same_date_proportion <
        0 |
    merged_stats$same_date_proportion >
        1,
    na.rm = TRUE
)) {
    stop(
        "Invalid same-date proportions were produced",
        call. = FALSE
    )
}


# ==================================== #
#              WRITE OUTPUT             #
# ==================================== #

fwrite(
    merged_stats,
    file = output_file,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)


# ==================================== #
#               SUMMARY                 #
# ==================================== #

message(
    "Final unique temporal code pairs: ",
    format_count(
        nrow(merged_stats)
    )
)

message(
    "Total participants-with-both observations: ",
    format_count(
        sum(
            merged_stats$participants_with_both
        )
    )
)

message(
    "Total directed observations: ",
    format_count(
        sum(
            merged_stats$temporally_ordered_count
        )
    )
)

message(
    "Total same-date observations: ",
    format_count(
        sum(
            merged_stats$same_date_count
        )
    )
)

message(
    "Wrote merged temporal summary: ",
    output_file
)