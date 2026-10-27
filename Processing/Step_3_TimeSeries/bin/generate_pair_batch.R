#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})


# ==================================== #
#              ARGUMENTS                #
# ==================================== #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
    stop(
        paste(
            "Expected 3 arguments:",
            "<participant_batch>",
            "<cohort_code>",
            "<output_file>"
        ),
        call. = FALSE
    )
}

participant_batch_file <- args[[1]]
cohort_code <- trimws(args[[2]])
output_file <- args[[3]]

if (!file.exists(participant_batch_file)) {
    stop(
        paste(
            "Participant batch file does not exist:",
            participant_batch_file
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
    "Processing participant batch: ",
    participant_batch_file
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


# ==================================== #
#              LOAD DATA                #
# ==================================== #

first_events <- fread(
    participant_batch_file,
    sep = ",",
    colClasses = c(
        eid = "character",
        code = "character",
        first_date = "character",
        last_date = "character",
        event_count = "integer",
        is_male = "integer",
        dob = "character",
        age_at_first_diagnosis = "numeric",
        age_at_last_diagnosis = "numeric"
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

message(
    "First-event rows loaded: ",
    format_count(nrow(first_events))
)

if (nrow(first_events) == 0L) {
    stop(
        "The participant batch is empty",
        call. = FALSE
    )
}


# ==================================== #
#          VALIDATE COLUMNS             #
# ==================================== #

required_columns <- c(
    "eid",
    "code",
    "first_date",
    "last_date",
    "event_count",
    "is_male",
    "dob",
    "age_at_first_diagnosis",
    "age_at_last_diagnosis"
)

missing_columns <- setdiff(
    required_columns,
    names(first_events)
)

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "Participant batch is missing required columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

# Ignore any additional columns that may have been added upstream.
first_events <- first_events[
    ,
    ..required_columns
]


# ==================================== #
#          STANDARDISE VALUES           #
# ==================================== #

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

first_events[
    ,
    last_date := as.IDate(
        last_date,
        format = "%Y-%m-%d"
    )
]

first_events[
    ,
    dob := as.IDate(
        dob,
        format = "%Y-%m-%d"
    )
]


# ==================================== #
#          VALIDATE VALUES              #
# ==================================== #

invalid_identifier_rows <- first_events[
    is.na(eid) |
    eid == "" |
    is.na(code) |
    code == ""
]

if (nrow(invalid_identifier_rows) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_identifier_rows)
            ),
            "rows with missing or empty eid or code values"
        ),
        call. = FALSE
    )
}

if (anyNA(first_events$first_date)) {
    stop(
        "Missing or invalid first_date values were found",
        call. = FALSE
    )
}

if (anyNA(first_events$last_date)) {
    stop(
        "Missing or invalid last_date values were found",
        call. = FALSE
    )
}

if (any(
    first_events$first_date >
        first_events$last_date
)) {
    stop(
        "Some first_date values occur after last_date",
        call. = FALSE
    )
}

if (any(
    !is.na(first_events$is_male) &
    !first_events$is_male %in% c(0L, 1L)
)) {
    stop(
        "is_male must contain only 0, 1, or NA",
        call. = FALSE
    )
}

if (any(
    first_events$event_count < 1L,
    na.rm = TRUE
)) {
    stop(
        "event_count contains values below one",
        call. = FALSE
    )
}


# ==================================== #
#       VALIDATE UNIQUE CODES           #
# ==================================== #

duplicate_participant_codes <- first_events[
    ,
    .N,
    by = .(
        eid,
        code
    )
][
    N > 1L
]

if (nrow(duplicate_participant_codes) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(duplicate_participant_codes)
            ),
            "duplicate eid-code combinations"
        ),
        call. = FALSE
    )
}


# ==================================== #
#     VALIDATE PARTICIPANT METADATA     #
# ==================================== #

participant_metadata_check <- first_events[
    ,
    .(
        distinct_nonmissing_sex_values =
            uniqueN(
                is_male[
                    !is.na(is_male)
                ]
            ),

        distinct_nonmissing_dob_values =
            uniqueN(
                dob[
                    !is.na(dob)
                ]
            )
    ),
    by = eid
]

inconsistent_sex <- participant_metadata_check[
    distinct_nonmissing_sex_values > 1L
]

if (nrow(inconsistent_sex) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(inconsistent_sex)
            ),
            "participants with inconsistent is_male values"
        ),
        call. = FALSE
    )
}

inconsistent_dob <- participant_metadata_check[
    distinct_nonmissing_dob_values > 1L
]

if (nrow(inconsistent_dob) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(inconsistent_dob)
            ),
            "participants with inconsistent dob values"
        ),
        call. = FALSE
    )
}


# ==================================== #
#       VALIDATE DIAGNOSIS AGES         #
# ==================================== #

invalid_negative_ages <- first_events[
    !is.na(age_at_first_diagnosis) &
    age_at_first_diagnosis < 0
]

if (nrow(invalid_negative_ages) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_negative_ages)
            ),
            "records with a negative age at first diagnosis"
        ),
        call. = FALSE
    )
}

implausible_ages <- first_events[
    !is.na(age_at_first_diagnosis) &
    age_at_first_diagnosis > 120
]

if (nrow(implausible_ages) > 0L) {
    warning(
        paste(
            "Found",
            format_count(
                nrow(implausible_ages)
            ),
            "records with age at first diagnosis above 120"
        ),
        call. = FALSE
    )
}


# ==================================== #
#              SORT DATA                #
# ==================================== #

setorder(
    first_events,
    eid,
    first_date,
    code
)


# ==================================== #
#      EXPECTED PAIR WORKLOAD           #
# ==================================== #

participant_work <- first_events[
    ,
    .(
        code_count = .N,

        expected_pair_count = (
            .N * (.N - 1L)
        ) / 2
    ),
    by = eid
]

participants_with_fewer_than_two_codes <- participant_work[
    code_count < 2L
]

if (
    nrow(participants_with_fewer_than_two_codes) >
    0L
) {
    warning(
        paste(
            "Found",
            format_count(
                nrow(
                    participants_with_fewer_than_two_codes
                )
            ),
            "participants with fewer than two codes.",
            "They will not produce temporal pairs."
        ),
        call. = FALSE
    )
}

expected_total_pairs <- sum(
    participant_work$expected_pair_count
)

if (expected_total_pairs < 1) {
    stop(
        "No participant in this batch can generate a code pair",
        call. = FALSE
    )
}

message(
    "Participants in batch: ",
    format_count(
        nrow(participant_work)
    )
)

message(
    "Expected participant-level pairs: ",
    format_count(
        expected_total_pairs
    )
)


# ==================================== #
#      GENERATE ONE PERSON'S PAIRS      #
# ==================================== #

generate_person_pairs <- function(
    codes,
    dates,
    ages,
    participant_sex
) {
    number_of_codes <- length(codes)

    if (number_of_codes < 2L) {
        return(NULL)
    }

    pair_indices <- combn(
        number_of_codes,
        2L
    )

    index_a <- pair_indices[1L, ]
    index_b <- pair_indices[2L, ]

    code_a <- codes[index_a]
    code_b <- codes[index_b]

    date_a <- dates[index_a]
    date_b <- dates[index_b]

    age_a <- ages[index_a]
    age_b <- ages[index_b]

    number_of_pairs <- length(index_a)

    # Canonical alphabetical code ordering.
    a_is_code_1 <- (
        code_a <= code_b
    )

    code_1 <- ifelse(
        a_is_code_1,
        code_a,
        code_b
    )

    code_2 <- ifelse(
        a_is_code_1,
        code_b,
        code_a
    )

    # Dates and ages corresponding to canonical code_1 and code_2.
    code_1_date <- as.IDate(
        rep(
            NA_integer_,
            number_of_pairs
        )
    )

    code_2_date <- as.IDate(
        rep(
            NA_integer_,
            number_of_pairs
        )
    )

    code_1_date[a_is_code_1] <-
        date_a[a_is_code_1]

    code_1_date[!a_is_code_1] <-
        date_b[!a_is_code_1]

    code_2_date[a_is_code_1] <-
        date_b[a_is_code_1]

    code_2_date[!a_is_code_1] <-
        date_a[!a_is_code_1]

    code_1_age <- rep(
        NA_real_,
        number_of_pairs
    )

    code_2_age <- rep(
        NA_real_,
        number_of_pairs
    )

    code_1_age[a_is_code_1] <-
        age_a[a_is_code_1]

    code_1_age[!a_is_code_1] <-
        age_b[!a_is_code_1]

    code_2_age[a_is_code_1] <-
        age_b[a_is_code_1]

    code_2_age[!a_is_code_1] <-
        age_a[!a_is_code_1]

    # Temporal ordering.
    same_date <- (
        date_a == date_b
    )

    a_before_b <- (
        date_a < date_b
    )

    b_before_a <- (
        date_b < date_a
    )

    code_from <- rep(
        NA_character_,
        number_of_pairs
    )

    code_to <- rep(
        NA_character_,
        number_of_pairs
    )

    code_from[a_before_b] <-
        code_a[a_before_b]

    code_to[a_before_b] <-
        code_b[a_before_b]

    code_from[b_before_a] <-
        code_b[b_before_a]

    code_to[b_before_a] <-
        code_a[b_before_a]

    date_from <- as.IDate(
        rep(
            NA_integer_,
            number_of_pairs
        )
    )

    date_to <- as.IDate(
        rep(
            NA_integer_,
            number_of_pairs
        )
    )

    date_from[a_before_b] <-
        date_a[a_before_b]

    date_to[a_before_b] <-
        date_b[a_before_b]

    date_from[b_before_a] <-
        date_b[b_before_a]

    date_to[b_before_a] <-
        date_a[b_before_a]

    age_from <- rep(
        NA_real_,
        number_of_pairs
    )

    age_to <- rep(
        NA_real_,
        number_of_pairs
    )

    age_from[a_before_b] <-
        age_a[a_before_b]

    age_to[a_before_b] <-
        age_b[a_before_b]

    age_from[b_before_a] <-
        age_b[b_before_a]

    age_to[b_before_a] <-
        age_a[b_before_a]

    days_between <- as.integer(
        abs(
            date_b -
            date_a
        )
    )

    direction <- rep(
        "same_date",
        number_of_pairs
    )

    direction[a_before_b] <- paste0(
        code_a[a_before_b],
        "_before_",
        code_b[a_before_b]
    )

    direction[b_before_a] <- paste0(
        code_b[b_before_a],
        "_before_",
        code_a[b_before_a]
    )

    data.table(
        is_male = rep(
            participant_sex[[1L]],
            number_of_pairs
        ),

        code_1 = code_1,
        code_2 = code_2,

        code_1_date = code_1_date,
        code_2_date = code_2_date,

        code_1_age = code_1_age,
        code_2_age = code_2_age,

        code_from = code_from,
        code_to = code_to,

        date_from = date_from,
        date_to = date_to,

        age_from = age_from,
        age_to = age_to,

        days_between = days_between,
        same_date = same_date,
        direction = direction
    )
}


# ==================================== #
#       GENERATE ALL BATCH PAIRS        #
# ==================================== #

message("Generating participant-level temporal pairs...")

all_temporal_pairs <- first_events[
    ,
    generate_person_pairs(
        codes = code,
        dates = first_date,
        ages = age_at_first_diagnosis,
        participant_sex = is_male
    ),
    by = eid
]

if (nrow(all_temporal_pairs) == 0L) {
    stop(
        "No temporal pairs were generated",
        call. = FALSE
    )
}

all_temporal_pairs[
    ,
    cohort_code := cohort_code
]

setcolorder(
    all_temporal_pairs,
    c(
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
)

setorder(
    all_temporal_pairs,
    eid,
    code_1,
    code_2
)


# ==================================== #
#          OUTPUT VALIDATION            #
# ==================================== #

if (
    nrow(all_temporal_pairs) !=
    expected_total_pairs
) {
    stop(
        paste(
            "Generated",
            format_count(
                nrow(all_temporal_pairs)
            ),
            "pairs but expected",
            format_count(
                expected_total_pairs
            )
        ),
        call. = FALSE
    )
}

duplicate_generated_pairs <- all_temporal_pairs[
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

if (nrow(duplicate_generated_pairs) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(duplicate_generated_pairs)
            ),
            "duplicate participant-level code pairs"
        ),
        call. = FALSE
    )
}

invalid_same_date_pairs <- all_temporal_pairs[
    same_date &
    (
        days_between != 0L |
        !is.na(code_from) |
        !is.na(code_to) |
        !is.na(date_from) |
        !is.na(date_to) |
        !is.na(age_from) |
        !is.na(age_to) |
        direction != "same_date"
    )
]

if (nrow(invalid_same_date_pairs) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_same_date_pairs)
            ),
            "invalid same-date pairs"
        ),
        call. = FALSE
    )
}

invalid_ordered_pairs <- all_temporal_pairs[
    !same_date &
    (
        days_between <= 0L |
        is.na(code_from) |
        is.na(code_to) |
        is.na(date_from) |
        is.na(date_to) |
        direction == "same_date"
    )
]

if (nrow(invalid_ordered_pairs) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(invalid_ordered_pairs)
            ),
            "invalid temporally ordered pairs"
        ),
        call. = FALSE
    )
}

if (any(
    all_temporal_pairs$code_1 >
        all_temporal_pairs$code_2
)) {
    stop(
        paste(
            "Canonical pair validation failed:",
            "code_1 is not always alphabetically before code_2"
        ),
        call. = FALSE
    )
}

date_order_error <- all_temporal_pairs[
    !same_date &
    date_from >= date_to
]

if (nrow(date_order_error) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(date_order_error)
            ),
            "ordered pairs where date_from is not before date_to"
        ),
        call. = FALSE
    )
}


# ==================================== #
#              WRITE OUTPUT             #
# ==================================== #

fwrite(
    all_temporal_pairs,
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
    "Participants represented: ",
    format_count(
        uniqueN(all_temporal_pairs$eid)
    )
)

message(
    "Participant-level pairs generated: ",
    format_count(
        nrow(all_temporal_pairs)
    )
)

message(
    "Unique batch-level code pairs: ",
    format_count(
        uniqueN(
            all_temporal_pairs,
            by = c(
                "code_1",
                "code_2"
            )
        )
    )
)

message(
    "Code 1 before code 2 pairs: ",
    format_count(
        sum(
            all_temporal_pairs$code_1_date <
            all_temporal_pairs$code_2_date
        )
    )
)

message(
    "Code 2 before code 1 pairs: ",
    format_count(
        sum(
            all_temporal_pairs$code_2_date <
            all_temporal_pairs$code_1_date
        )
    )
)

message(
    "Same-date pairs: ",
    format_count(
        sum(
            all_temporal_pairs$same_date
        )
    )
)

message(
    "Wrote temporal pair batch: ",
    output_file
)