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
            "<diagnosis_events>",
            "<participant_metadata>",
            "<cohort_code>",
            "<clean_events_output>",
            "<first_events_output>",
            "<code_summary_output>",
            "<cohort_summary_output>"
        ),
        call. = FALSE
    )
}

diagnosis_events_file <- args[[1]]
participant_metadata_file <- args[[2]]
cohort_code <- args[[3]]
clean_events_output <- args[[4]]
first_events_output <- args[[5]]
code_summary_output <- args[[6]]
cohort_summary_output <- args[[7]]

if (!file.exists(diagnosis_events_file)) {
    stop(
        paste(
            "Diagnosis event file does not exist:",
            diagnosis_events_file
        ),
        call. = FALSE
    )
}

if (!file.exists(participant_metadata_file)) {
    stop(
        paste(
            "Participant metadata file does not exist:",
            participant_metadata_file
        ),
        call. = FALSE
    )
}

if (
    is.na(cohort_code) ||
    !nzchar(trimws(cohort_code))
) {
    stop(
        "cohort_code cannot be empty",
        call. = FALSE
    )
}

cohort_code <- trimws(cohort_code)

message("Preparing cohort: ", cohort_code)
message("Diagnosis events: ", diagnosis_events_file)
message("Participant metadata: ", participant_metadata_file)


# ==================================== #
#         HELPER FUNCTIONS              #
# ==================================== #

add_calendar_months <- function(
    date_value,
    number_of_months = 1L
) {
    date_value <- as.Date(date_value)

    if (
        length(date_value) != 1L ||
        is.na(date_value)
    ) {
        return(as.IDate(NA_integer_))
    }

    original_year <- as.integer(
        format(
            date_value,
            "%Y"
        )
    )

    original_month <- as.integer(
        format(
            date_value,
            "%m"
        )
    )

    original_day <- as.integer(
        format(
            date_value,
            "%d"
        )
    )

    month_index <- (
        original_year * 12L +
        original_month -
        1L +
        as.integer(number_of_months)
    )

    target_year <- month_index %/% 12L
    target_month <- month_index %% 12L + 1L

    first_day_target_month <- as.Date(
        sprintf(
            "%04d-%02d-01",
            target_year,
            target_month
        )
    )

    first_day_next_month <- seq.Date(
        from = first_day_target_month,
        by = "month",
        length.out = 2L
    )[[2L]]

    last_day_target_month <- (
        first_day_next_month - 1L
    )

    last_day_number <- as.integer(
        format(
            last_day_target_month,
            "%d"
        )
    )

    target_day <- min(
        original_day,
        last_day_number
    )

    as.IDate(
        sprintf(
            "%04d-%02d-%02d",
            target_year,
            target_month,
            target_day
        ),
        format = "%Y-%m-%d"
    )
}


format_count <- function(value) {
    format(
        value,
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
    )
}


# ==================================== #
#      LOAD PARTICIPANT METADATA        #
# ==================================== #

participants <- fread(
    participant_metadata_file,
    sep = ",",
    colClasses = c(
        eid = "character",
        is_male = "integer",
        dob = "character"
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

required_participant_columns <- c(
    "eid",
    "is_male",
    "dob"
)

missing_participant_columns <- setdiff(
    required_participant_columns,
    names(participants)
)

if (length(missing_participant_columns) > 0L) {
    stop(
        paste(
            "Participant metadata is missing required columns:",
            paste(
                missing_participant_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

# Keep only columns used by this workflow.
participants <- participants[
    ,
    ..required_participant_columns
]

message(
    "Participant metadata rows loaded: ",
    format_count(nrow(participants))
)


# ==================================== #
#     CLEAN PARTICIPANT METADATA        #
# ==================================== #

participants[
    ,
    eid := trimws(eid)
]

participants[
    ,
    dob := trimws(dob)
]

participants[
    ,
    dob := as.IDate(
        dob,
        format = "%Y-%m-%d"
    )
]

invalid_eid_rows <- participants[
    is.na(eid) |
    eid == ""
]

if (nrow(invalid_eid_rows) > 0L) {
    stop(
        paste(
            "Participant metadata contains",
            format_count(nrow(invalid_eid_rows)),
            "rows with missing or empty eid values"
        ),
        call. = FALSE
    )
}

invalid_sex_rows <- participants[
    !is.na(is_male) &
    !is_male %in% c(0L, 1L)
]

if (nrow(invalid_sex_rows) > 0L) {
    stop(
        paste(
            "Participant metadata contains",
            format_count(nrow(invalid_sex_rows)),
            "rows where is_male is not 0, 1, or NA"
        ),
        call. = FALSE
    )
}

# Remove completely identical participant rows.
participant_rows_before_deduplication <- nrow(
    participants
)

participants <- unique(
    participants,
    by = c(
        "eid",
        "is_male",
        "dob"
    )
)

identical_participant_rows_removed <- (
    participant_rows_before_deduplication -
    nrow(participants)
)

message(
    "Identical participant rows removed: ",
    format_count(
        identical_participant_rows_removed
    )
)

# Any duplicate eid values remaining contain conflicting or repeated
# participant information and are therefore rejected.
duplicate_participant_ids <- participants[
    ,
    .N,
    by = eid
][
    N > 1L
]

if (nrow(duplicate_participant_ids) > 0L) {
    stop(
        paste(
            "Participant metadata contains",
            format_count(
                nrow(duplicate_participant_ids)
            ),
            "duplicate eid values after removing identical rows"
        ),
        call. = FALSE
    )
}

cohort_size <- uniqueN(
    participants$eid
)

male_count <- participants[
    is_male == 1L,
    .N
]

female_count <- participants[
    is_male == 0L,
    .N
]

missing_sex_count <- participants[
    is.na(is_male),
    .N
]

participants_with_dob <- participants[
    !is.na(dob),
    .N
]

participants_missing_dob <- participants[
    is.na(dob),
    .N
]

minimum_dob <- if (participants_with_dob > 0L) {
    min(
        participants$dob,
        na.rm = TRUE
    )
} else {
    as.IDate(NA_integer_)
}

maximum_dob <- if (participants_with_dob > 0L) {
    max(
        participants$dob,
        na.rm = TRUE
    )
} else {
    as.IDate(NA_integer_)
}

message(
    "Cohort size: ",
    format_count(cohort_size)
)

message(
    "Sex counts: male=",
    format_count(male_count),
    "; female=",
    format_count(female_count),
    "; missing=",
    format_count(missing_sex_count)
)

message(
    "Date of birth available: ",
    format_count(participants_with_dob),
    "; missing: ",
    format_count(participants_missing_dob)
)


# ==================================== #
#       LOAD DIAGNOSIS EVENTS           #
# ==================================== #

events <- fread(
    diagnosis_events_file,
    sep = ",",
    colClasses = c(
        eid = "character",
        code = "character",
        date = "character"
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

required_event_columns <- c(
    "eid",
    "code",
    "date"
)

missing_event_columns <- setdiff(
    required_event_columns,
    names(events)
)

if (length(missing_event_columns) > 0L) {
    stop(
        paste(
            "Diagnosis event file is missing required columns:",
            paste(
                missing_event_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

# Keep only the required event columns.
events <- events[
    ,
    ..required_event_columns
]

raw_event_count <- nrow(events)

message(
    "Raw diagnosis event rows: ",
    format_count(raw_event_count)
)


# ==================================== #
#       CLEAN DIAGNOSIS EVENTS          #
# ==================================== #

events[
    ,
    eid := trimws(eid)
]

events[
    ,
    code := toupper(
        trimws(code)
    )
]

events[
    ,
    date := trimws(date)
]

events[
    ,
    date := as.IDate(
        date,
        format = "%Y-%m-%d"
    )
]

invalid_event_rows <- events[
    is.na(eid) |
    eid == "" |
    is.na(code) |
    code == "" |
    is.na(date)
]

invalid_event_count <- nrow(
    invalid_event_rows
)

message(
    "Invalid diagnosis rows removed: ",
    format_count(invalid_event_count)
)

events <- events[
    !is.na(eid) &
    eid != "" &
    !is.na(code) &
    code != "" &
    !is.na(date)
]

if (nrow(events) == 0L) {
    stop(
        "No valid diagnosis events remain after cleaning",
        call. = FALSE
    )
}


# ==================================== #
#      REMOVE EXACT DUPLICATES          #
# ==================================== #

# Exact duplicate event rows are the same eid, code, and date.
# Removing them prevents duplicate source records from inflating
# event_count.

events_before_deduplication <- nrow(
    events
)

events <- unique(
    events,
    by = c(
        "eid",
        "code",
        "date"
    )
)

exact_duplicate_events_removed <- (
    events_before_deduplication -
    nrow(events)
)

message(
    "Exact duplicate diagnosis events removed: ",
    format_count(
        exact_duplicate_events_removed
    )
)


# ==================================== #
#       VALIDATE PARTICIPANT IDS        #
# ==================================== #

participant_ids <- participants$eid

event_eids_not_in_metadata <- setdiff(
    unique(events$eid),
    participant_ids
)

if (length(event_eids_not_in_metadata) > 0L) {
    stop(
        paste(
            "Diagnosis events contain",
            format_count(
                length(event_eids_not_in_metadata)
            ),
            "eid values not found in participant metadata"
        ),
        call. = FALSE
    )
}


# ==================================== #
#      ATTACH PARTICIPANT METADATA      #
# ==================================== #

events <- merge(
    events,
    participants,
    by = "eid",
    all.x = TRUE,
    sort = FALSE
)

if (anyNA(events$eid)) {
    stop(
        "Internal join error produced missing eid values",
        call. = FALSE
    )
}


# ==================================== #
#      VALIDATE DIAGNOSIS AGES          #
# ==================================== #

events[
    ,
    age_at_event := as.numeric(
        date - dob
    ) / 365.25
]

events[
    ,
    diagnosis_before_birth := (
        !is.na(dob) &
        date < dob
    )
]

diagnoses_before_birth <- events[
    diagnosis_before_birth == TRUE
]

if (nrow(diagnoses_before_birth) > 0L) {
    stop(
        paste(
            "Found",
            format_count(
                nrow(diagnoses_before_birth)
            ),
            "diagnosis events recorded before date of birth"
        ),
        call. = FALSE
    )
}

implausible_age_events <- events[
    !is.na(age_at_event) &
    age_at_event > 120
]

if (nrow(implausible_age_events) > 0L) {
    warning(
        paste(
            "Found",
            format_count(
                nrow(implausible_age_events)
            ),
            "diagnosis events at ages above 120 years"
        ),
        call. = FALSE
    )
}

events[
    ,
    diagnosis_before_birth := NULL
]


# ==================================== #
#        OBSERVATION PERIOD             #
# ==================================== #

# Temporary cohort-wide observation window:
#
# Start = earliest valid event date.
# End   = latest valid event date plus one calendar month.

observation_start <- min(
    events$date
)

latest_event_date <- max(
    events$date
)

observation_end <- add_calendar_months(
    latest_event_date,
    number_of_months = 1L
)

if (
    is.na(observation_end) ||
    observation_end <= observation_start
) {
    stop(
        "Could not derive a valid cohort observation period",
        call. = FALSE
    )
}

observation_days <- as.integer(
    observation_end -
    observation_start
)

observation_years <- (
    observation_days /
    365.25
)

message(
    "Observation period: ",
    observation_start,
    " to ",
    observation_end,
    " (",
    format_count(observation_days),
    " days; ",
    round(
        observation_years,
        2
    ),
    " years)"
)


# ==================================== #
#        CLEAN EVENT OUTPUT             #
# ==================================== #

clean_events <- events[
    ,
    .(
        eid,
        code,
        date
    )
]

setorder(
    clean_events,
    eid,
    date,
    code
)


# ==================================== #
#       COLLAPSE FIRST EVENTS           #
# ==================================== #

# One row per participant-code combination.
#
# event_count is the number of unique recorded dates after exact
# duplicate event rows have been removed.

first_events <- events[
    ,
    .(
        first_date = min(date),
        last_date = max(date),
        event_count = .N
    ),
    by = .(
        eid,
        code,
        is_male,
        dob
    )
]

first_events[
    ,
    age_at_first_diagnosis := as.numeric(
        first_date - dob
    ) / 365.25
]

first_events[
    ,
    age_at_last_diagnosis := as.numeric(
        last_date - dob
    ) / 365.25
]

setorder(
    first_events,
    eid,
    first_date,
    code
)


# ==================================== #
#       VALIDATE FIRST EVENTS           #
# ==================================== #

duplicate_first_events <- first_events[
    ,
    .N,
    by = .(
        eid,
        code
    )
][
    N > 1L
]

if (nrow(duplicate_first_events) > 0L) {
    stop(
        paste(
            "Internal validation failed:",
            format_count(
                nrow(duplicate_first_events)
            ),
            "duplicate eid-code combinations remain"
        ),
        call. = FALSE
    )
}

if (any(
    first_events$first_date >
        first_events$last_date
)) {
    stop(
        paste(
            "Internal validation failed:",
            "first_date occurs after last_date"
        ),
        call. = FALSE
    )
}

if (any(first_events$event_count < 1L)) {
    stop(
        paste(
            "Internal validation failed:",
            "event_count contains values below one"
        ),
        call. = FALSE
    )
}

participants_with_any_code <- uniqueN(
    first_events$eid
)

participants_without_codes <- (
    cohort_size -
    participants_with_any_code
)

if (participants_without_codes < 0L) {
    stop(
        paste(
            "Internal validation failed:",
            "more participants occur in first_events",
            "than in participant metadata"
        ),
        call. = FALSE
    )
}


# ==================================== #
#          CODE SUMMARY                 #
# ==================================== #

# first_events contains one row per participant and code.
# Therefore, .N gives participants_with_code.

code_summary <- first_events[
    ,
    .(
        participants_with_code = .N,

        total_recorded_events = sum(
            event_count
        ),

        first_recorded_date = min(
            first_date
        ),

        last_recorded_date = max(
            last_date
        ),

        mean_age_at_first_diagnosis =
            if (all(
                is.na(age_at_first_diagnosis)
            )) {
                NA_real_
            } else {
                mean(
                    age_at_first_diagnosis,
                    na.rm = TRUE
                )
            },

        median_age_at_first_diagnosis =
            if (all(
                is.na(age_at_first_diagnosis)
            )) {
                NA_real_
            } else {
                median(
                    age_at_first_diagnosis,
                    na.rm = TRUE
                )
            },

        male_participants_with_code =
            sum(
                is_male == 1L,
                na.rm = TRUE
            ),

        female_participants_with_code =
            sum(
                is_male == 0L,
                na.rm = TRUE
            ),

        participants_with_missing_sex =
            sum(
                is.na(is_male)
            )
    ),
    by = code
]

code_summary[
    ,
    cohort_code := cohort_code
]

code_summary[
    ,
    cohort_size := cohort_size
]

code_summary[
    ,
    participants_without_code := (
        cohort_size -
        participants_with_code
    )
]

code_summary[
    ,
    code_prevalence := (
        participants_with_code /
        cohort_size
    )
]

code_summary[
    ,
    code_prevalence_percent := (
        code_prevalence * 100
    )
]

setcolorder(
    code_summary,
    c(
        "cohort_code",
        "code",
        "cohort_size",
        "participants_with_code",
        "participants_without_code",
        "code_prevalence",
        "code_prevalence_percent",
        "total_recorded_events",
        "first_recorded_date",
        "last_recorded_date",
        "mean_age_at_first_diagnosis",
        "median_age_at_first_diagnosis",
        "male_participants_with_code",
        "female_participants_with_code",
        "participants_with_missing_sex"
    )
)

setorder(
    code_summary,
    -participants_with_code,
    code
)


# ==================================== #
#         COHORT SUMMARY                #
# ==================================== #

cohort_summary <- data.table(
    cohort_code = cohort_code,

    cohort_size = cohort_size,

    participants_with_any_code =
        participants_with_any_code,

    participants_without_codes =
        participants_without_codes,

    valid_event_rows =
        nrow(clean_events),

    invalid_event_rows_removed =
        invalid_event_count,

    exact_duplicate_events_removed =
        exact_duplicate_events_removed,

    unique_codes =
        uniqueN(first_events$code),

    participant_code_combinations =
        nrow(first_events),

    observation_start =
        observation_start,

    latest_event_date =
        latest_event_date,

    observation_end =
        observation_end,

    observation_days =
        observation_days,

    observation_years =
        observation_years,

    male_count =
        male_count,

    female_count =
        female_count,

    missing_sex_count =
        missing_sex_count,

    participants_with_dob =
        participants_with_dob,

    participants_missing_dob =
        participants_missing_dob,

    minimum_dob =
        minimum_dob,

    maximum_dob =
        maximum_dob
)


# ==================================== #
#       FINAL OUTPUT VALIDATION         #
# ==================================== #

if (any(
    code_summary$participants_with_code >
        cohort_size
)) {
    stop(
        paste(
            "Internal validation failed:",
            "a code contains more participants",
            "than the cohort size"
        ),
        call. = FALSE
    )
}

if (any(
    code_summary$participants_without_code <
        0L
)) {
    stop(
        paste(
            "Internal validation failed:",
            "participants_without_code contains",
            "negative values"
        ),
        call. = FALSE
    )
}

if (any(
    code_summary$code_prevalence < 0 |
    code_summary$code_prevalence > 1
)) {
    stop(
        paste(
            "Internal validation failed:",
            "code_prevalence is outside 0 to 1"
        ),
        call. = FALSE
    )
}


# ==================================== #
#              WRITE OUTPUTS            #
# ==================================== #

fwrite(
    clean_events,
    file = clean_events_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)

fwrite(
    first_events,
    file = first_events_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)

fwrite(
    code_summary,
    file = code_summary_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)

fwrite(
    cohort_summary,
    file = cohort_summary_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)


# ==================================== #
#                 SUMMARY               #
# ==================================== #

message(
    "Valid cleaned events: ",
    format_count(nrow(clean_events))
)

message(
    "Participants with at least one code: ",
    format_count(
        participants_with_any_code
    )
)

message(
    "Participants without recorded codes: ",
    format_count(
        participants_without_codes
    )
)

message(
    "Unique participant-code combinations: ",
    format_count(nrow(first_events))
)

message(
    "Unique codes: ",
    format_count(nrow(code_summary))
)

message(
    "Wrote cleaned events: ",
    clean_events_output
)

message(
    "Wrote first events: ",
    first_events_output
)

message(
    "Wrote code summary: ",
    code_summary_output
)

message(
    "Wrote cohort summary: ",
    cohort_summary_output
)