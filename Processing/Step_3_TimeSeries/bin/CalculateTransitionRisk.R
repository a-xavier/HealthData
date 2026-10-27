#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})


# ==================================== #
#              ARGUMENTS                #
# ==================================== #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7L) {
    stop(
        paste(
            "Expected 7 arguments:",
            "<first_events>",
            "<code_summary>",
            "<cohort_summary>",
            "<temporal_summary>",
            "<participant_metadata>",
            "<cohort_code>",
            "<output_file>"
        ),
        call. = FALSE
    )
}

first_events_file <- args[[1]]
code_summary_file <- args[[2]]
cohort_summary_file <- args[[3]]
temporal_summary_file <- args[[4]]
participant_metadata_file <- args[[5]]
cohort_code <- trimws(args[[6]])
output_file <- args[[7]]

input_files <- c(
    first_events_file,
    code_summary_file,
    cohort_summary_file,
    temporal_summary_file,
    participant_metadata_file
)

missing_files <- input_files[
    !file.exists(input_files)
]

if (length(missing_files) > 0L) {
    stop(
        paste(
            "Input files do not exist:",
            paste(
                missing_files,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

if (
    is.na(cohort_code) ||
    cohort_code == ""
) {
    stop(
        "cohort_code cannot be empty",
        call. = FALSE
    )
}


# ==================================== #
#              HELPERS                  #
# ==================================== #

format_count <- function(value) {
    format(
        value,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
    )
}


read_required_columns <- function(
    input_file,
    required_columns,
    ...
) {
    result <- fread(
        input_file,
        showProgress = FALSE,
        nThread = 1L,
        check.names = FALSE,
        ...
    )

    missing_columns <- setdiff(
        required_columns,
        names(result)
    )

    if (length(missing_columns) > 0L) {
        stop(
            paste(
                basename(input_file),
                "is missing required columns:",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }

    result
}


add_binomial_intervals <- function(
    data,
    events_column,
    eligible_column,
    prefix
) {
    events <- as.numeric(
        data[[events_column]]
    )

    eligible <- as.numeric(
        data[[eligible_column]]
    )

    valid <- (
        !is.na(events) &
        !is.na(eligible) &
        eligible > 0 &
        events >= 0 &
        events <= eligible
    )

    risk <- rep(
        NA_real_,
        nrow(data)
    )

    lower <- rep(
        NA_real_,
        nrow(data)
    )

    upper <- rep(
        NA_real_,
        nrow(data)
    )

    risk[valid] <- (
        events[valid] /
        eligible[valid]
    )

    valid_zero <- (
        valid &
        events == 0
    )

    valid_nonzero <- (
        valid &
        events > 0
    )

    valid_complete <- (
        valid &
        events == eligible
    )

    valid_incomplete <- (
        valid &
        events < eligible
    )

    lower[valid_zero] <- 0

    lower[valid_nonzero] <- qbeta(
        0.025,
        events[valid_nonzero],
        eligible[valid_nonzero] -
            events[valid_nonzero] +
            1
    )

    upper[valid_complete] <- 1

    upper[valid_incomplete] <- qbeta(
        0.975,
        events[valid_incomplete] + 1,
        eligible[valid_incomplete] -
            events[valid_incomplete]
    )

    data[
        ,
        (paste0(prefix, "_risk")) :=
            risk
    ]

    data[
        ,
        (paste0(prefix, "_risk_ci_low")) :=
            lower
    ]

    data[
        ,
        (paste0(prefix, "_risk_ci_high")) :=
            upper
    ]
}


add_poisson_rate_intervals <- function(
    data,
    events_column,
    person_years_column,
    prefix,
    multiplier = 1000
) {
    events <- as.numeric(
        data[[events_column]]
    )

    person_years <- as.numeric(
        data[[person_years_column]]
    )

    valid <- (
        !is.na(events) &
        !is.na(person_years) &
        events >= 0 &
        person_years > 0
    )

    rate <- rep(
        NA_real_,
        nrow(data)
    )

    lower <- rep(
        NA_real_,
        nrow(data)
    )

    upper <- rep(
        NA_real_,
        nrow(data)
    )

    rate[valid] <- (
        events[valid] /
        person_years[valid]
    ) * multiplier

    lower_event_limit <- rep(
        NA_real_,
        nrow(data)
    )

    lower_event_limit[
        valid &
        events == 0
    ] <- 0

    nonzero <- (
        valid &
        events > 0
    )

    lower_event_limit[nonzero] <- (
        0.5 *
        qchisq(
            0.025,
            2 * events[nonzero]
        )
    )

    upper_event_limit <- rep(
        NA_real_,
        nrow(data)
    )

    upper_event_limit[valid] <- (
        0.5 *
        qchisq(
            0.975,
            2 * (
                events[valid] + 1
            )
        )
    )

    lower[valid] <- (
        lower_event_limit[valid] /
        person_years[valid]
    ) * multiplier

    upper[valid] <- (
        upper_event_limit[valid] /
        person_years[valid]
    ) * multiplier

    data[
        ,
        (paste0(
            prefix,
            "_incidence_per_1000_person_years"
        )) := rate
    ]

    data[
        ,
        (paste0(
            prefix,
            "_incidence_ci_low"
        )) := lower
    ]

    data[
        ,
        (paste0(
            prefix,
            "_incidence_ci_high"
        )) := upper
    ]
}


# ==================================== #
#          LOAD COHORT SUMMARY          #
# ==================================== #

cohort_summary <- read_required_columns(
    cohort_summary_file,
    required_columns = c(
        "cohort_code",
        "cohort_size",
        "observation_start",
        "observation_end"
    ),
    sep = ",",
    colClasses = c(
        cohort_code = "character",
        cohort_size = "integer",
        observation_start = "character",
        observation_end = "character"
    )
)

if (nrow(cohort_summary) != 1L) {
    stop(
        "cohort_summary must contain exactly one row",
        call. = FALSE
    )
}

cohort_summary[
    ,
    cohort_code := trimws(cohort_code)
]

if (
    cohort_summary$cohort_code[[1]] !=
    cohort_code
) {
    stop(
        "cohort_summary does not match cohort_code",
        call. = FALSE
    )
}

cohort_size <- cohort_summary$cohort_size[[1]]

observation_start <- as.IDate(
    cohort_summary$observation_start[[1]],
    format = "%Y-%m-%d"
)

observation_end <- as.IDate(
    cohort_summary$observation_end[[1]],
    format = "%Y-%m-%d"
)

if (
    is.na(observation_start) ||
    is.na(observation_end) ||
    observation_end <= observation_start
) {
    stop(
        "Invalid cohort observation period",
        call. = FALSE
    )
}

message("Cohort: ", cohort_code)
message("Cohort size: ", format_count(cohort_size))

message(
    "Observation period: ",
    observation_start,
    " to ",
    observation_end
)


# ==================================== #
#       VALIDATE PARTICIPANTS           #
# ==================================== #

participants <- read_required_columns(
    participant_metadata_file,
    required_columns = c(
        "eid",
        "is_male",
        "dob"
    ),
    sep = ",",
    colClasses = c(
        eid = "character",
        is_male = "integer",
        dob = "character"
    )
)

participants[
    ,
    eid := trimws(eid)
]

if (anyDuplicated(participants$eid)) {
    stop(
        "participant_metadata contains duplicate eid values",
        call. = FALSE
    )
}

if (
    uniqueN(participants$eid) !=
    cohort_size
) {
    stop(
        paste(
            "participant_metadata contains",
            uniqueN(participants$eid),
            "participants but cohort_summary reports",
            cohort_size
        ),
        call. = FALSE
    )
}


# ==================================== #
#          LOAD FIRST EVENTS            #
# ==================================== #

first_events <- read_required_columns(
    first_events_file,
    required_columns = c(
        "eid",
        "code",
        "first_date",
        "is_male",
        "age_at_first_diagnosis"
    ),
    sep = ",",
    colClasses = c(
        eid = "character",
        code = "character",
        first_date = "character",
        is_male = "integer",
        age_at_first_diagnosis = "numeric"
    )
)

first_events[
    ,
    eid := trimws(eid)
]

first_events[
    ,
    code := toupper(
        trimws(code)
    )
]

first_events[
    ,
    first_date := as.IDate(
        first_date,
        format = "%Y-%m-%d"
    )
]

if (anyNA(first_events$first_date)) {
    stop(
        "Invalid first_date values were found",
        call. = FALSE
    )
}

if (anyDuplicated(
    first_events,
    by = c(
        "eid",
        "code"
    )
)) {
    stop(
        "first_events contains duplicate eid-code combinations",
        call. = FALSE
    )
}


# ==================================== #
#      SOURCE-CODE FOLLOW-UP TOTALS     #
# ==================================== #

first_events[
    ,
    followup_days_after_code := pmax(
        0,
        as.numeric(
            observation_end -
            first_date
        )
    )
]

code_followup <- first_events[
    ,
    .(
        source_total_participants = .N,

        source_positive_followup_count =
            sum(
                followup_days_after_code >
                0
            ),

        source_total_followup_days =
            sum(
                followup_days_after_code
            ),

        source_mean_age_at_diagnosis =
            if (
                all(
                    is.na(
                        age_at_first_diagnosis
                    )
                )
            ) {
                NA_real_
            } else {
                mean(
                    age_at_first_diagnosis,
                    na.rm = TRUE
                )
            },

        source_median_age_at_diagnosis =
            if (
                all(
                    is.na(
                        age_at_first_diagnosis
                    )
                )
            ) {
                NA_real_
            } else {
                median(
                    age_at_first_diagnosis,
                    na.rm = TRUE
                )
            },

        source_male_count =
            sum(
                is_male == 1L,
                na.rm = TRUE
            ),

        source_female_count =
            sum(
                is_male == 0L,
                na.rm = TRUE
            ),

        source_missing_sex_count =
            sum(
                is.na(is_male)
            )
    ),
    by = code
]


# ==================================== #
#          LOAD CODE SUMMARY            #
# ==================================== #

code_summary <- read_required_columns(
    code_summary_file,
    required_columns = c(
        "cohort_code",
        "code",
        "participants_with_code",
        "code_prevalence"
    ),
    sep = ",",
    colClasses = c(
        cohort_code = "character",
        code = "character",
        participants_with_code = "integer",
        code_prevalence = "numeric"
    )
)

code_summary[
    ,
    code := toupper(
        trimws(code)
    )
]

code_summary[
    ,
    cohort_code := trimws(cohort_code)
]

if (any(
    code_summary$cohort_code !=
        cohort_code
)) {
    stop(
        "code_summary contains a different cohort_code",
        call. = FALSE
    )
}


# ==================================== #
#        LOAD TEMPORAL SUMMARY          #
# ==================================== #

temporal_summary <- fread(
    temporal_summary_file,
    sep = ",",
    na.strings = c(
        "NA",
        ""
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

required_temporal_summary_columns <- c(
    "cohort_code",
    "code_1",
    "code_2",
    "participants_with_both",
    "code_1_before_code_2_count",
    "code_2_before_code_1_count",
    "same_date_count",
    "temporally_ordered_count"
)

missing_columns <- setdiff(
    required_temporal_summary_columns,
    names(temporal_summary)
)

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "temporal_summary is missing columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

temporal_summary[
    ,
    code_1 := toupper(
        trimws(code_1)
    )
]

temporal_summary[
    ,
    code_2 := toupper(
        trimws(code_2)
    )
]


# ==================================== #
#    DISCOVER TEMPORAL-PAIR BATCHES     #
# ==================================== #

temporal_pair_files <- sort(
    Sys.glob(
        "temporal_pairs_*.csv"
    )
)

if (length(temporal_pair_files) == 0L) {
    stop(
        paste(
            "No temporal_pairs_*.csv files",
            "were found in the process directory"
        ),
        call. = FALSE
    )
}

message(
    "Temporal-pair batch files: ",
    format_count(
        length(temporal_pair_files)
    )
)


# ==================================== #
#    CALCULATE PAIR-SPECIFIC OFFSETS    #
# ==================================== #

pair_adjustment_list <- vector(
    mode = "list",
    length = length(temporal_pair_files)
)

for (
    file_index in
    seq_along(temporal_pair_files)
) {
    input_file <- temporal_pair_files[
        file_index
    ]

    message(
        "Processing temporal-pair file ",
        file_index,
        " of ",
        length(temporal_pair_files),
        ": ",
        basename(input_file)
    )

    pair_data <- fread(
        input_file,
        sep = ",",
        na.strings = c(
            "NA",
            ""
        ),
        colClasses = c(
            cohort_code = "character",
            eid = "character",
            code_1 = "character",
            code_2 = "character",
            code_1_date = "character",
            code_2_date = "character",
            is_male = "integer"
        ),
        showProgress = FALSE,
        nThread = 1L,
        check.names = FALSE
    )

    required_pair_columns <- c(
        "cohort_code",
        "eid",
        "code_1",
        "code_2",
        "code_1_date",
        "code_2_date",
        "is_male"
    )

    missing_columns <- setdiff(
        required_pair_columns,
        names(pair_data)
    )

    if (length(missing_columns) > 0L) {
        stop(
            paste(
                basename(input_file),
                "is missing columns:",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }

    pair_data[
        ,
        code_1 := toupper(
            trimws(code_1)
        )
    ]

    pair_data[
        ,
        code_2 := toupper(
            trimws(code_2)
        )
    ]

    pair_data[
        ,
        code_1_date := as.IDate(
            code_1_date,
            format = "%Y-%m-%d"
        )
    ]

    pair_data[
        ,
        code_2_date := as.IDate(
            code_2_date,
            format = "%Y-%m-%d"
        )
    ]

    if (
        anyNA(pair_data$code_1_date) ||
        anyNA(pair_data$code_2_date)
    ) {
        stop(
            paste(
                "Invalid diagnosis dates in",
                basename(input_file)
            ),
            call. = FALSE
        )
    }

    pair_data[
        ,
        code_1_before_code_2 := (
            code_1_date <
            code_2_date
        )
    ]

    pair_data[
        ,
        code_2_before_code_1 := (
            code_2_date <
            code_1_date
        )
    ]

    pair_data[
        ,
        same_date := (
            code_1_date ==
            code_2_date
        )
    ]

    pair_adjustment_list[[
        file_index
    ]] <- pair_data[
        ,
        .(
            # -------------------------------- #
            # Code 1 -> code 2                 #
            # -------------------------------- #

            code_1_to_code_2_event_count =
                sum(
                    code_1_before_code_2
                ),

            code_1_to_code_2_ineligible_count =
                sum(
                    code_2_before_code_1 |
                    same_date
                ),

            code_1_to_code_2_ineligible_followup_days =
                sum(
                    pmax(
                        0,
                        as.numeric(
                            observation_end -
                            code_1_date[
                                code_2_before_code_1 |
                                same_date
                            ]
                        )
                    )
                ),

            code_1_to_code_2_post_event_days =
                sum(
                    pmax(
                        0,
                        as.numeric(
                            observation_end -
                            code_2_date[
                                code_1_before_code_2
                            ]
                        )
                    )
                ),

            # -------------------------------- #
            # Code 2 -> code 1                 #
            # -------------------------------- #

            code_2_to_code_1_event_count =
                sum(
                    code_2_before_code_1
                ),

            code_2_to_code_1_ineligible_count =
                sum(
                    code_1_before_code_2 |
                    same_date
                ),

            code_2_to_code_1_ineligible_followup_days =
                sum(
                    pmax(
                        0,
                        as.numeric(
                            observation_end -
                            code_2_date[
                                code_1_before_code_2 |
                                same_date
                            ]
                        )
                    )
                ),

            code_2_to_code_1_post_event_days =
                sum(
                    pmax(
                        0,
                        as.numeric(
                            observation_end -
                            code_1_date[
                                code_2_before_code_1
                            ]
                        )
                    )
                )
        ),
        by = .(
            code_1,
            code_2
        )
    ]

    rm(pair_data)
    invisible(gc())
}

pair_adjustments <- rbindlist(
    pair_adjustment_list,
    use.names = TRUE,
    fill = FALSE
)

rm(pair_adjustment_list)
invisible(gc())


# ==================================== #
#       MERGE BATCH ADJUSTMENTS         #
# ==================================== #

pair_adjustments <- pair_adjustments[
    ,
    lapply(
        .SD,
        sum
    ),
    by = .(
        code_1,
        code_2
    ),
    .SDcols = setdiff(
        names(pair_adjustments),
        c(
            "code_1",
            "code_2"
        )
    )
]


# ==================================== #
#      VALIDATE EVENT COUNTS            #
# ==================================== #

count_validation <- merge(
    temporal_summary[
        ,
        .(
            code_1,
            code_2,
            code_1_before_code_2_count,
            code_2_before_code_1_count,
            same_date_count
        )
    ],
    pair_adjustments,
    by = c(
        "code_1",
        "code_2"
    ),
    all.x = TRUE,
    sort = FALSE
)

if (anyNA(
    count_validation$
        code_1_to_code_2_event_count
)) {
    stop(
        paste(
            "Some temporal-summary pairs were not found",
            "in the participant-level pair files"
        ),
        call. = FALSE
    )
}

if (any(
    count_validation$
        code_1_before_code_2_count !=
    count_validation$
        code_1_to_code_2_event_count
)) {
    stop(
        paste(
            "Code 1 to code 2 event counts differ",
            "between temporal_summary and pair batches"
        ),
        call. = FALSE
    )
}

if (any(
    count_validation$
        code_2_before_code_1_count !=
    count_validation$
        code_2_to_code_1_event_count
)) {
    stop(
        paste(
            "Code 2 to code 1 event counts differ",
            "between temporal_summary and pair batches"
        ),
        call. = FALSE
    )
}

rm(count_validation)


# ==================================== #
#       BUILD DIRECTED TRANSITIONS      #
# ==================================== #

# ------------------------------------ #
#      Code 1 -> code 2                 #
# ------------------------------------ #

directed_1_to_2 <- merge(
    temporal_summary,
    pair_adjustments,
    by = c(
        "code_1",
        "code_2"
    ),
    all.x = TRUE,
    sort = FALSE
)

directed_1_to_2[
    ,
    source_code := code_1
]

directed_1_to_2[
    ,
    target_code := code_2
]

directed_1_to_2[
    ,
    transition_event_count :=
        code_1_before_code_2_count
]

directed_1_to_2[
    ,
    reverse_order_count :=
        code_2_before_code_1_count
]

directed_1_to_2[
    ,
    transition_label := paste0(
        source_code,
        "_to_",
        target_code
    )
]


# ------------------------------------ #
#      Code 2 -> code 1                 #
# ------------------------------------ #

directed_2_to_1 <- copy(
    directed_1_to_2
)

directed_2_to_1[
    ,
    source_code := code_2
]

directed_2_to_1[
    ,
    target_code := code_1
]

directed_2_to_1[
    ,
    transition_event_count :=
        code_2_before_code_1_count
]

directed_2_to_1[
    ,
    reverse_order_count :=
        code_1_before_code_2_count
]

directed_2_to_1[
    ,
    transition_label := paste0(
        source_code,
        "_to_",
        target_code
    )
]


# ==================================== #
#       SELECT DIRECTION VALUES         #
# ==================================== #

directional_fields <- function(
    data,
    direction
) {
    if (direction == "1_to_2") {

        data[
            ,
            ineligible_count :=
                code_1_to_code_2_ineligible_count
        ]

        data[
            ,
            ineligible_followup_days :=
                code_1_to_code_2_ineligible_followup_days
        ]

        data[
            ,
            post_event_days :=
                code_1_to_code_2_post_event_days
        ]

        prefix <- "code_1_to_code_2"

    } else {

        data[
            ,
            ineligible_count :=
                code_2_to_code_1_ineligible_count
        ]

        data[
            ,
            ineligible_followup_days :=
                code_2_to_code_1_ineligible_followup_days
        ]

        data[
            ,
            post_event_days :=
                code_2_to_code_1_post_event_days
        ]

        prefix <- "code_2_to_code_1"
    }

    copied_suffixes <- c(
        "min_days",
        "q1_days",
        "median_days",
        "mean_days",
        "sd_days",
        "q3_days",
        "max_days",

        "within_7_days",
        "within_30_days",
        "within_90_days",
        "within_1_year",
        "within_2_years",
        "within_5_years",
        "within_10_years",

        "within_7_days_proportion",
        "within_30_days_proportion",
        "within_90_days_proportion",
        "within_1_year_proportion",
        "within_2_years_proportion",
        "within_5_years_proportion",
        "within_10_years_proportion",

        "source_age_count",
        "source_age_mean",
        "source_age_sd",
        "source_age_min",
        "source_age_q1",
        "source_age_median",
        "source_age_q3",
        "source_age_max",

        "target_age_count",
        "target_age_mean",
        "target_age_sd",
        "target_age_q1",
        "target_age_median",
        "target_age_q3"
    )

    for (suffix in copied_suffixes) {

        source_column <- paste0(
            prefix,
            "_",
            suffix
        )

        if (source_column %in% names(data)) {
            data[
                ,
                (suffix) := get(
                    source_column
                )
            ]
        }
    }

    data
}


directed_1_to_2 <- directional_fields(
    directed_1_to_2,
    "1_to_2"
)

directed_2_to_1 <- directional_fields(
    directed_2_to_1,
    "2_to_1"
)


# ==================================== #
#       COMBINE DIRECTED RESULTS        #
# ==================================== #

directed_results <- rbindlist(
    list(
        directed_1_to_2,
        directed_2_to_1
    ),
    use.names = TRUE,
    fill = TRUE
)

rm(
    directed_1_to_2,
    directed_2_to_1,
    pair_adjustments
)

invisible(gc())


# ==================================== #
#       ATTACH SOURCE FOLLOW-UP         #
# ==================================== #

source_followup <- copy(
    code_followup
)

setnames(
    source_followup,
    "code",
    "source_code"
)

directed_results <- merge(
    directed_results,
    source_followup,
    by = "source_code",
    all.x = TRUE,
    sort = FALSE
)


# ==================================== #
#       ATTACH TARGET TOTALS            #
# ==================================== #

target_summary <- code_summary[
    ,
    .(
        target_code = code,
        target_total_participants =
            participants_with_code,
        target_prevalence =
            code_prevalence
    )
]

directed_results <- merge(
    directed_results,
    target_summary,
    by = "target_code",
    all.x = TRUE,
    sort = FALSE
)


# ==================================== #
#          ELIGIBLE POPULATION          #
# ==================================== #

# Eligible for source -> target:
#
# 1. Participant has the source code.
# 2. Target was not recorded before source.
# 3. Target was not recorded on the source date.
#
# Under the current cohort-wide cutoff assumption, every source
# diagnosis has positive follow-up because observation_end is one
# month after the latest diagnosis date.

directed_results[
    ,
    eligible_source_count := (
        source_total_participants -
        ineligible_count
    )
]

directed_results[
    ,
    expected_ineligible_count := (
        reverse_order_count +
        same_date_count
    )
]

ineligible_mismatch <- directed_results[
    ineligible_count !=
        expected_ineligible_count
]

if (nrow(ineligible_mismatch) > 0L) {
    fwrite(
        ineligible_mismatch,
        file = "invalid_ineligible_counts.csv",
        sep = ",",
        quote = FALSE,
        na = "NA"
    )

    stop(
        paste(
            nrow(ineligible_mismatch),
            "directed transitions have inconsistent ineligible counts.",
            "Details written to invalid_ineligible_counts.csv"
        ),
        call. = FALSE
    )
}

invalid_eligibility <- directed_results[
    eligible_source_count < 0L |
    transition_event_count >
        eligible_source_count
]

if (nrow(invalid_eligibility) > 0L) {

    diagnostic_columns <- c(
        "source_code",
        "target_code",
        "source_total_participants",
        "source_positive_followup_count",
        "participants_with_both",
        "transition_event_count",
        "reverse_order_count",
        "same_date_count",
        "ineligible_count",
        "expected_ineligible_count",
        "eligible_source_count"
    )

    diagnostic_columns <- intersect(
        diagnostic_columns,
        names(invalid_eligibility)
    )

    fwrite(
        invalid_eligibility[
            ,
            ..diagnostic_columns
        ],
        file = "invalid_transition_eligibility.csv",
        sep = ",",
        quote = FALSE,
        na = "NA"
    )

    print(
        invalid_eligibility[
            ,
            ..diagnostic_columns
        ][
            1:min(.N, 20L)
        ]
    )

    stop(
        paste(
            nrow(invalid_eligibility),
            "directed transitions have invalid eligible populations.",
            "Details written to invalid_transition_eligibility.csv"
        ),
        call. = FALSE
    )
}

directed_results[
    ,
    expected_ineligible_count := NULL
]


# ==================================== #
#          PERSON-TIME AT RISK          #
# ==================================== #

# Start with all available follow-up after the source diagnosis.
#
# Remove:
# 1. all source follow-up from participants whose target diagnosis
#    was recorded before or on the source date;
# 2. follow-up after the target diagnosis among transition events.

directed_results[
    ,
    person_days_at_risk := pmax(
        0,
        source_total_followup_days -
        ineligible_followup_days -
        post_event_days
    )
]

directed_results[
    ,
    person_years_at_risk := (
        person_days_at_risk /
        365.25
    )
]


# ==================================== #
#       ABSOLUTE OBSERVED RISK          #
# ==================================== #

add_binomial_intervals(
    data = directed_results,
    events_column = "transition_event_count",
    eligible_column = "eligible_source_count",
    prefix = "observed_transition"
)


# ==================================== #
#            INCIDENCE RATE             #
# ==================================== #

add_poisson_rate_intervals(
    data = directed_results,
    events_column = "transition_event_count",
    person_years_column = "person_years_at_risk",
    prefix = "transition",
    multiplier = 1000
)


# ==================================== #
#     CRUDE RATE PER SOURCE PERSON      #
# ==================================== #

directed_results[
    ,
    transition_events_per_1000_source_participants :=
        fifelse(
            source_total_participants > 0L,

            transition_event_count /
                source_total_participants *
                1000,

            NA_real_
        )
]


# ==================================== #
#        COHORT PROVENANCE              #
# ==================================== #

directed_results[
    ,
    cohort_code := cohort_code
]

directed_results[
    ,
    cohort_size := cohort_size
]

directed_results[
    ,
    observation_start := observation_start
]

directed_results[
    ,
    observation_end := observation_end
]


# ==================================== #
#        SELECT FINAL COLUMNS           #
# ==================================== #

preferred_columns <- c(
    "cohort_code",
    "source_code",
    "target_code",
    "transition_label",

    "cohort_size",
    "observation_start",
    "observation_end",

    "source_total_participants",
    "target_total_participants",
    "target_prevalence",
    "participants_with_both",

    "source_positive_followup_count",
    "eligible_source_count",
    "ineligible_count",

    "transition_event_count",
    "reverse_order_count",
    "same_date_count",

    "observed_transition_risk",
    "observed_transition_risk_ci_low",
    "observed_transition_risk_ci_high",

    "person_days_at_risk",
    "person_years_at_risk",

    "transition_incidence_per_1000_person_years",
    "transition_incidence_ci_low",
    "transition_incidence_ci_high",

    "transition_events_per_1000_source_participants",

    "min_days",
    "q1_days",
    "median_days",
    "mean_days",
    "sd_days",
    "q3_days",
    "max_days",

    "within_7_days",
    "within_30_days",
    "within_90_days",
    "within_1_year",
    "within_2_years",
    "within_5_years",
    "within_10_years",

    "within_7_days_proportion",
    "within_30_days_proportion",
    "within_90_days_proportion",
    "within_1_year_proportion",
    "within_2_years_proportion",
    "within_5_years_proportion",
    "within_10_years_proportion",

    "source_age_count",
    "source_age_mean",
    "source_age_sd",
    "source_age_min",
    "source_age_q1",
    "source_age_median",
    "source_age_q3",
    "source_age_max",

    "target_age_count",
    "target_age_mean",
    "target_age_sd",
    "target_age_q1",
    "target_age_median",
    "target_age_q3",

    "source_male_count",
    "source_female_count",
    "source_missing_sex_count",

    "temporally_ordered_count",
    "code_1_first_among_ordered",
    "code_1_first_ci_low",
    "code_1_first_ci_high",
    "direction_balance",
    "direction_odds_corrected",
    "log_direction_odds",
    "direction_p_value",
    "direction_fdr",
    "dominant_recorded_direction"
)

final_columns <- intersect(
    preferred_columns,
    names(directed_results)
)

directed_results <- directed_results[
    ,
    ..final_columns
]


# ==================================== #
#              FINAL CHECKS             #
# ==================================== #

if (any(
    directed_results$
        observed_transition_risk < 0 |
    directed_results$
        observed_transition_risk > 1,
    na.rm = TRUE
)) {
    stop(
        "Observed transition risks outside 0 to 1 were produced",
        call. = FALSE
    )
}

if (any(
    directed_results$
        person_days_at_risk < 0,
    na.rm = TRUE
)) {
    stop(
        "Negative person-time was produced",
        call. = FALSE
    )
}

if (anyDuplicated(
    directed_results,
    by = c(
        "cohort_code",
        "source_code",
        "target_code"
    )
)) {
    stop(
        "Duplicate directed transition rows were produced",
        call. = FALSE
    )
}


# ==================================== #
#              SORT AND WRITE           #
# ==================================== #

setorder(
    directed_results,
    -transition_event_count,
    source_code,
    target_code
)

fwrite(
    directed_results,
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
    "Directed transition results: ",
    format_count(
        nrow(directed_results)
    )
)

message(
    "Transitions with at least one event: ",
    format_count(
        directed_results[
            transition_event_count > 0L,
            .N
        ]
    )
)

message(
    "Total directed transition events: ",
    format_count(
        sum(
            directed_results$
                transition_event_count
        )
    )
)

message(
    "Wrote transition-risk results: ",
    output_file
)

