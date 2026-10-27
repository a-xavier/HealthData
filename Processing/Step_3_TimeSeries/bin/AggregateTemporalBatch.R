#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})


# ==================================== #
#              ARGUMENTS                #
# ==================================== #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4L) {
    stop(
        paste(
            "Expected 4 arguments:",
            "<temporal_pair_batch>",
            "<cohort_code>",
            "<partial_stats_output>",
            "<partial_delays_output>"
        ),
        call. = FALSE
    )
}

temporal_pair_batch_file <- args[[1]]
cohort_code <- trimws(args[[2]])
partial_stats_output <- args[[3]]
partial_delays_output <- args[[4]]

if (!file.exists(temporal_pair_batch_file)) {
    stop(
        paste(
            "Temporal-pair batch does not exist:",
            temporal_pair_batch_file
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

message(
    "Aggregating temporal-pair batch: ",
    temporal_pair_batch_file
)

message(
    "Cohort: ",
    cohort_code
)


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


safe_minimum <- function(values) {

    values <- values[
        !is.na(values)
    ]

    if (length(values) == 0L) {
        return(NA_real_)
    }

    as.numeric(min(values))
}


safe_maximum <- function(values) {

    values <- values[
        !is.na(values)
    ]

    if (length(values) == 0L) {
        return(NA_real_)
    }

    as.numeric(max(values))
}


# ==================================== #
#              LOAD DATA                #
# ==================================== #

temporal_pairs <- fread(
    temporal_pair_batch_file,
    sep = ",",
    na.strings = c(
        "NA",
        ""
    ),
    colClasses = list(
        character = c(
            "cohort_code",
            "eid",
            "code_1",
            "code_2",
            "code_from",
            "code_to",
            "direction"
        ),
        integer = "is_male",
        character = c(
            "code_1_date",
            "code_2_date",
            "date_from",
            "date_to"
        ),
        numeric = c(
            "code_1_age",
            "code_2_age",
            "age_from",
            "age_to"
        ),
        integer = "days_between",
        logical = "same_date"
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

message(
    "Participant-level pair rows loaded: ",
    format_count(nrow(temporal_pairs))
)

if (nrow(temporal_pairs) == 0L) {
    stop(
        "The temporal-pair batch is empty",
        call. = FALSE
    )
}


# ==================================== #
#          VALIDATE COLUMNS             #
# ==================================== #

required_columns <- c(
    "cohort_code",
    "eid",
    "is_male",
    "code_1",
    "code_2",
    "code_1_date",
    "code_2_date",
    "code_1_age",
    "code_2_age",
    "code_from",
    "code_to",
    "date_from",
    "date_to",
    "age_from",
    "age_to",
    "days_between",
    "same_date",
    "direction"
)

missing_columns <- setdiff(
    required_columns,
    names(temporal_pairs)
)

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "Temporal-pair batch is missing required columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

temporal_pairs <- temporal_pairs[
    ,
    ..required_columns
]


# ==================================== #
#          STANDARDISE VALUES           #
# ==================================== #

temporal_pairs[
    ,
    cohort_code := trimws(cohort_code)
]

temporal_pairs[
    ,
    eid := trimws(eid)
]

temporal_pairs[
    ,
    code_1 := toupper(
        trimws(code_1)
    )
]

temporal_pairs[
    ,
    code_2 := toupper(
        trimws(code_2)
    )
]

temporal_pairs[
    !is.na(code_from),
    code_from := toupper(
        trimws(code_from)
    )
]

temporal_pairs[
    !is.na(code_to),
    code_to := toupper(
        trimws(code_to)
    )
]

temporal_pairs[
    ,
    code_1_date := as.IDate(
        code_1_date,
        format = "%Y-%m-%d"
    )
]

temporal_pairs[
    ,
    code_2_date := as.IDate(
        code_2_date,
        format = "%Y-%m-%d"
    )
]

temporal_pairs[
    ,
    date_from := as.IDate(
        date_from,
        format = "%Y-%m-%d"
    )
]

temporal_pairs[
    ,
    date_to := as.IDate(
        date_to,
        format = "%Y-%m-%d"
    )
]


# ==================================== #
#          VALIDATE COHORT              #
# ==================================== #

input_cohorts <- unique(
    temporal_pairs$cohort_code
)

if (
    length(input_cohorts) != 1L ||
    is.na(input_cohorts[[1]]) ||
    input_cohorts[[1]] != cohort_code
) {
    stop(
        paste(
            "Temporal-pair batch cohort does not match",
            "the supplied cohort_code"
        ),
        call. = FALSE
    )
}


# ==================================== #
#          VALIDATE VALUES              #
# ==================================== #

invalid_identifiers <- temporal_pairs[
    is.na(eid) |
    eid == "" |
    is.na(code_1) |
    code_1 == "" |
    is.na(code_2) |
    code_2 == ""
]

if (nrow(invalid_identifiers) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_identifiers)
            ),
            "rows with missing eid, code_1, or code_2"
        ),
        call. = FALSE
    )
}

if (any(
    temporal_pairs$code_1 ==
        temporal_pairs$code_2
)) {
    stop(
        "Self-pairs were found where code_1 equals code_2",
        call. = FALSE
    )
}

if (any(
    temporal_pairs$code_1 >
        temporal_pairs$code_2
)) {
    stop(
        paste(
            "Canonical pair ordering failed:",
            "code_1 must be alphabetically before code_2"
        ),
        call. = FALSE
    )
}

if (
    anyNA(temporal_pairs$code_1_date) ||
    anyNA(temporal_pairs$code_2_date)
) {
    stop(
        paste(
            "Missing code_1_date or code_2_date",
            "values were found"
        ),
        call. = FALSE
    )
}

if (anyNA(temporal_pairs$days_between)) {
    stop(
        "Missing days_between values were found",
        call. = FALSE
    )
}

if (any(
    temporal_pairs$days_between < 0L
)) {
    stop(
        "Negative days_between values were found",
        call. = FALSE
    )
}

if (anyNA(temporal_pairs$same_date)) {
    stop(
        "Missing same_date values were found",
        call. = FALSE
    )
}

if (any(
    !is.na(temporal_pairs$is_male) &
    !temporal_pairs$is_male %in% c(0L, 1L)
)) {
    stop(
        "is_male must contain only 0, 1, or NA",
        call. = FALSE
    )
}


# ==================================== #
#       VALIDATE UNIQUE PAIR ROWS       #
# ==================================== #

duplicate_rows <- temporal_pairs[
    ,
    .N,
    by = .(
        eid,
        code_1,
        code_2
    )
][
    N > 1L
]

if (nrow(duplicate_rows) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(duplicate_rows)
            ),
            "duplicate participant-level code pairs"
        ),
        call. = FALSE
    )
}


# ==================================== #
#        CLASSIFY TEMPORAL ORDER        #
# ==================================== #

temporal_pairs[
    ,
    code_1_before_code_2 := (
        code_1_date <
        code_2_date
    )
]

temporal_pairs[
    ,
    code_2_before_code_1 := (
        code_2_date <
        code_1_date
    )
]

temporal_pairs[
    ,
    same_date := (
        code_1_date ==
        code_2_date
    )
]

temporal_pairs[
    ,
    days_between := as.integer(
        abs(
            code_2_date -
            code_1_date
        )
    )
]

temporal_pairs[
    ,
    classification_count := (
        as.integer(code_1_before_code_2) +
        as.integer(code_2_before_code_1) +
        as.integer(same_date)
    )
]

invalid_classification <- temporal_pairs[
    classification_count != 1L
]

if (nrow(invalid_classification) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_classification)
            ),
            "rows that could not be classified into exactly",
            "one temporal category"
        ),
        call. = FALSE
    )
}

temporal_pairs[
    ,
    classification_count := NULL
]


# ==================================== #
#        VALIDATE DIRECTION DATA        #
# ==================================== #

invalid_same_date <- temporal_pairs[
    same_date &
    (
        days_between != 0L |
        !is.na(code_from) |
        !is.na(code_to) |
        !is.na(date_from) |
        !is.na(date_to) |
        !is.na(age_from) |
        !is.na(age_to)
    )
]

if (nrow(invalid_same_date) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_same_date)
            ),
            "invalid same-date rows"
        ),
        call. = FALSE
    )
}

invalid_ordered <- temporal_pairs[
    !same_date &
    (
        days_between <= 0L |
        is.na(code_from) |
        is.na(code_to) |
        is.na(date_from) |
        is.na(date_to) |
        date_from >= date_to
    )
]

if (nrow(invalid_ordered) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_ordered)
            ),
            "invalid temporally ordered rows"
        ),
        call. = FALSE
    )
}


# ==================================== #
#       CREATE DIRECTION LABELS         #
# ==================================== #

temporal_pairs[
    ,
    temporal_class := fcase(
        code_1_before_code_2,
        "code_1_to_code_2",

        code_2_before_code_1,
        "code_2_to_code_1",

        same_date,
        "same_date",

        default = NA_character_
    )
]


# ==================================== #
#       AGGREGATE PARTIAL STATS         #
# ==================================== #

message(
    "Aggregating mergeable pair statistics..."
)

partial_stats <- temporal_pairs[
    ,
    {
        direction_1_to_2 <-
            code_1_before_code_2

        direction_2_to_1 <-
            code_2_before_code_1

        same_day <-
            same_date

        days_1_to_2 <- days_between[
            direction_1_to_2
        ]

        days_2_to_1 <- days_between[
            direction_2_to_1
        ]

        source_age_1_to_2 <- code_1_age[
            direction_1_to_2
        ]

        target_age_1_to_2 <- code_2_age[
            direction_1_to_2
        ]

        source_age_2_to_1 <- code_2_age[
            direction_2_to_1
        ]

        target_age_2_to_1 <- code_1_age[
            direction_2_to_1
        ]

        .(
            participants_with_both =
                .N,

            code_1_before_code_2_count =
                sum(direction_1_to_2),

            code_2_before_code_1_count =
                sum(direction_2_to_1),

            same_date_count =
                sum(same_day),

            temporally_ordered_count =
                sum(
                    direction_1_to_2 |
                    direction_2_to_1
                ),

            # ------------------------- #
            # Fixed-window counts       #
            # ------------------------- #

            code_1_to_code_2_within_7_days =
                sum(
                    direction_1_to_2 &
                    days_between <= 7L
                ),

            code_1_to_code_2_within_30_days =
                sum(
                    direction_1_to_2 &
                    days_between <= 30L
                ),

            code_1_to_code_2_within_90_days =
                sum(
                    direction_1_to_2 &
                    days_between <= 90L
                ),

            code_1_to_code_2_within_1_year =
                sum(
                    direction_1_to_2 &
                    days_between <= 365L
                ),

            code_1_to_code_2_within_2_years =
                sum(
                    direction_1_to_2 &
                    days_between <= 731L
                ),

            code_1_to_code_2_within_5_years =
                sum(
                    direction_1_to_2 &
                    days_between <= 1826L
                ),

            code_1_to_code_2_within_10_years =
                sum(
                    direction_1_to_2 &
                    days_between <= 3653L
                ),

            code_2_to_code_1_within_7_days =
                sum(
                    direction_2_to_1 &
                    days_between <= 7L
                ),

            code_2_to_code_1_within_30_days =
                sum(
                    direction_2_to_1 &
                    days_between <= 30L
                ),

            code_2_to_code_1_within_90_days =
                sum(
                    direction_2_to_1 &
                    days_between <= 90L
                ),

            code_2_to_code_1_within_1_year =
                sum(
                    direction_2_to_1 &
                    days_between <= 365L
                ),

            code_2_to_code_1_within_2_years =
                sum(
                    direction_2_to_1 &
                    days_between <= 731L
                ),

            code_2_to_code_1_within_5_years =
                sum(
                    direction_2_to_1 &
                    days_between <= 1826L
                ),

            code_2_to_code_1_within_10_years =
                sum(
                    direction_2_to_1 &
                    days_between <= 3653L
                ),

            # ------------------------- #
            # Mergeable time summaries  #
            # ------------------------- #

            code_1_to_code_2_sum_days =
                sum(
                    as.numeric(
                        days_1_to_2
                    )
                ),

            code_1_to_code_2_sum_squared_days =
                sum(
                    as.numeric(
                        days_1_to_2
                    )^2
                ),

            code_1_to_code_2_min_days =
                safe_minimum(
                    days_1_to_2
                ),

            code_1_to_code_2_max_days =
                safe_maximum(
                    days_1_to_2
                ),

            code_2_to_code_1_sum_days =
                sum(
                    as.numeric(
                        days_2_to_1
                    )
                ),

            code_2_to_code_1_sum_squared_days =
                sum(
                    as.numeric(
                        days_2_to_1
                    )^2
                ),

            code_2_to_code_1_min_days =
                safe_minimum(
                    days_2_to_1
                ),

            code_2_to_code_1_max_days =
                safe_maximum(
                    days_2_to_1
                ),

            # ------------------------- #
            # Sex-specific counts       #
            # ------------------------- #

            male_code_1_before_code_2 =
                sum(
                    direction_1_to_2 &
                    is_male == 1L,
                    na.rm = TRUE
                ),

            female_code_1_before_code_2 =
                sum(
                    direction_1_to_2 &
                    is_male == 0L,
                    na.rm = TRUE
                ),

            missing_sex_code_1_before_code_2 =
                sum(
                    direction_1_to_2 &
                    is.na(is_male)
                ),

            male_code_2_before_code_1 =
                sum(
                    direction_2_to_1 &
                    is_male == 1L,
                    na.rm = TRUE
                ),

            female_code_2_before_code_1 =
                sum(
                    direction_2_to_1 &
                    is_male == 0L,
                    na.rm = TRUE
                ),

            missing_sex_code_2_before_code_1 =
                sum(
                    direction_2_to_1 &
                    is.na(is_male)
                ),

            male_same_date =
                sum(
                    same_day &
                    is_male == 1L,
                    na.rm = TRUE
                ),

            female_same_date =
                sum(
                    same_day &
                    is_male == 0L,
                    na.rm = TRUE
                ),

            missing_sex_same_date =
                sum(
                    same_day &
                    is.na(is_male)
                ),

            # ------------------------- #
            # Source-age accumulators    #
            # ------------------------- #

            code_1_to_code_2_source_age_count =
                sum(
                    !is.na(
                        source_age_1_to_2
                    )
                ),

            code_1_to_code_2_source_age_sum =
                sum(
                    source_age_1_to_2,
                    na.rm = TRUE
                ),

            code_1_to_code_2_source_age_sum_squared =
                sum(
                    source_age_1_to_2^2,
                    na.rm = TRUE
                ),

            code_1_to_code_2_source_age_min =
                safe_minimum(
                    source_age_1_to_2
                ),

            code_1_to_code_2_source_age_max =
                safe_maximum(
                    source_age_1_to_2
                ),

            code_1_to_code_2_target_age_count =
                sum(
                    !is.na(
                        target_age_1_to_2
                    )
                ),

            code_1_to_code_2_target_age_sum =
                sum(
                    target_age_1_to_2,
                    na.rm = TRUE
                ),

            code_1_to_code_2_target_age_sum_squared =
                sum(
                    target_age_1_to_2^2,
                    na.rm = TRUE
                ),

            code_2_to_code_1_source_age_count =
                sum(
                    !is.na(
                        source_age_2_to_1
                    )
                ),

            code_2_to_code_1_source_age_sum =
                sum(
                    source_age_2_to_1,
                    na.rm = TRUE
                ),

            code_2_to_code_1_source_age_sum_squared =
                sum(
                    source_age_2_to_1^2,
                    na.rm = TRUE
                ),

            code_2_to_code_1_source_age_min =
                safe_minimum(
                    source_age_2_to_1
                ),

            code_2_to_code_1_source_age_max =
                safe_maximum(
                    source_age_2_to_1
                ),

            code_2_to_code_1_target_age_count =
                sum(
                    !is.na(
                        target_age_2_to_1
                    )
                ),

            code_2_to_code_1_target_age_sum =
                sum(
                    target_age_2_to_1,
                    na.rm = TRUE
                ),

            code_2_to_code_1_target_age_sum_squared =
                sum(
                    target_age_2_to_1^2,
                    na.rm = TRUE
                )
        )
    },
    by = .(
        cohort_code,
        code_1,
        code_2
    )
]


# ==================================== #
#         VALIDATE PARTIAL STATS        #
# ==================================== #

partial_stats[
    ,
    count_validation := (
        code_1_before_code_2_count +
        code_2_before_code_1_count +
        same_date_count
    )
]

invalid_partial_counts <- partial_stats[
    count_validation !=
        participants_with_both
]

if (nrow(invalid_partial_counts) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_partial_counts)
            ),
            "partial pair summaries that failed count validation"
        ),
        call. = FALSE
    )
}

partial_stats[
    ,
    count_validation := NULL
]


# ==================================== #
#         CREATE DELAY OUTPUT           #
# ==================================== #

# This compact table preserves the values needed for exact quantiles
# during the final merge.
#
# Same-date observations are excluded because their delay is exactly
# zero and they do not represent a directed transition.

partial_delays <- temporal_pairs[
    same_date == FALSE,
    .(
        cohort_code,
        code_1,
        code_2,
        temporal_class,
        days_between,
        age_from,
        age_to,
        is_male
    )
]

setorder(
    partial_delays,
    code_1,
    code_2,
    temporal_class,
    days_between
)


# ==================================== #
#          OUTPUT VALIDATION            #
# ==================================== #

expected_delay_rows <- sum(
    temporal_pairs$code_1_before_code_2 |
    temporal_pairs$code_2_before_code_1
)

if (
    nrow(partial_delays) !=
    expected_delay_rows
) {
    stop(
        paste(
            "Delay output contains",
            format_count(
                nrow(partial_delays)
            ),
            "rows but expected",
            format_count(
                expected_delay_rows
            )
        ),
        call. = FALSE
    )
}

if (any(
    partial_delays$days_between <= 0L
)) {
    stop(
        "The directed delay output contains non-positive delays",
        call. = FALSE
    )
}


# ==================================== #
#              SORT OUTPUT              #
# ==================================== #

setorder(
    partial_stats,
    -participants_with_both,
    code_1,
    code_2
)


# ==================================== #
#              WRITE OUTPUTS            #
# ==================================== #

fwrite(
    partial_stats,
    file = partial_stats_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)

fwrite(
    partial_delays,
    file = partial_delays_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)


# ==================================== #
#               SUMMARY                 #
# ==================================== #

message(
    "Input participant-level pairs: ",
    format_count(
        nrow(temporal_pairs)
    )
)

message(
    "Unique code pairs in partial summary: ",
    format_count(
        nrow(partial_stats)
    )
)

message(
    "Directed delay rows retained: ",
    format_count(
        nrow(partial_delays)
    )
)

message(
    "Same-date rows represented only in aggregate counts: ",
    format_count(
        sum(temporal_pairs$same_date)
    )
)

message(
    "Wrote partial statistics: ",
    partial_stats_output
)

message(
    "Wrote partial directed delays: ",
    partial_delays_output
)