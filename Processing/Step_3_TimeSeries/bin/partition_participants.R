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
            "<first_events_file>",
            "<cohort_code>",
            "<target_pairs_per_batch>",
            "<manifest_output>"
        ),
        call. = FALSE
    )
}

first_events_file <- args[[1]]
cohort_code <- trimws(args[[2]])
target_pairs_per_batch <- as.numeric(args[[3]])
manifest_output <- args[[4]]

if (!file.exists(first_events_file)) {
    stop(
        paste(
            "First-events file does not exist:",
            first_events_file
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

if (
    is.na(target_pairs_per_batch) ||
    !is.finite(target_pairs_per_batch) ||
    target_pairs_per_batch < 1
) {
    stop(
        "target_pairs_per_batch must be a positive number",
        call. = FALSE
    )
}


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
    first_events_file,
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
            "First-events file is missing required columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

if (nrow(first_events) == 0L) {
    stop(
        "The first-events file is empty",
        call. = FALSE
    )
}


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

invalid_rows <- first_events[
    is.na(eid) |
    eid == "" |
    is.na(code) |
    code == ""
]

if (nrow(invalid_rows) > 0L) {
    stop(
        paste(
            "Found",
            format_count(nrow(invalid_rows)),
            "rows with missing or empty eid or code values"
        ),
        call. = FALSE
    )
}


# ==================================== #
#       VALIDATE PARTICIPANT CODES      #
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
#        CALCULATE PARTICIPANT WORK     #
# ==================================== #

participant_work <- first_events[
    ,
    .(
        code_count = .N,

        estimated_pair_count = (
            .N * (.N - 1)
        ) / 2,

        first_recorded_date = min(
            first_date
        ),

        last_recorded_date = max(
            last_date
        )
    ),
    by = eid
]

total_participants <- nrow(
    participant_work
)

participants_with_pairs <- participant_work[
    estimated_pair_count > 0
]

participants_without_pairs <- participant_work[
    estimated_pair_count == 0
]

message(
    "Participants in first-events file: ",
    format_count(total_participants)
)

message(
    "Participants with at least two codes: ",
    format_count(
        nrow(participants_with_pairs)
    )
)

message(
    "Participants with fewer than two codes: ",
    format_count(
        nrow(participants_without_pairs)
    )
)

if (nrow(participants_with_pairs) == 0L) {
    stop(
        paste(
            "No participants have at least two codes;",
            "no temporal pairs can be generated"
        ),
        call. = FALSE
    )
}

total_estimated_pairs <- sum(
    participants_with_pairs$estimated_pair_count
)

maximum_participant_work <- max(
    participants_with_pairs$estimated_pair_count
)

message(
    "Total estimated participant-level pairs: ",
    format_count(total_estimated_pairs)
)

message(
    "Largest individual participant workload: ",
    format_count(maximum_participant_work)
)


# ==================================== #
#       DETERMINE NUMBER OF BATCHES     #
# ==================================== #

number_of_batches <- max(
    1L,
    as.integer(
        ceiling(
            total_estimated_pairs /
            target_pairs_per_batch
        )
    )
)

# There cannot be more useful batches than participants.
number_of_batches <- min(
    number_of_batches,
    nrow(participants_with_pairs)
)

message(
    "Target pairs per batch: ",
    format_count(target_pairs_per_batch)
)

message(
    "Number of batches: ",
    format_count(number_of_batches)
)

if (
    maximum_participant_work >
    target_pairs_per_batch
) {
    warning(
        paste(
            "At least one participant has an estimated workload",
            "larger than target_pairs_per_batch.",
            "That participant will remain intact in one batch."
        ),
        call. = FALSE
    )
}


# ==================================== #
#        ASSIGN PARTICIPANTS            #
# ==================================== #

# Largest-workload-first greedy assignment:
#
# 1. Sort participants from highest to lowest workload.
# 2. Assign each participant to the batch with the current lowest
#    estimated workload.
#
# Sorting eid as the second key makes the assignment deterministic
# when participants have equal workloads.

setorder(
    participants_with_pairs,
    -estimated_pair_count,
    eid
)

batch_workloads <- numeric(
    number_of_batches
)

batch_participant_counts <- integer(
    number_of_batches
)

assigned_batches <- integer(
    nrow(participants_with_pairs)
)

for (
    participant_index in
    seq_len(nrow(participants_with_pairs))
) {
    selected_batch <- which.min(
        batch_workloads
    )

    assigned_batches[
        participant_index
    ] <- selected_batch

    participant_workload <- (
        participants_with_pairs$
            estimated_pair_count[
                participant_index
            ]
    )

    batch_workloads[
        selected_batch
    ] <- (
        batch_workloads[
            selected_batch
        ] +
        participant_workload
    )

    batch_participant_counts[
        selected_batch
    ] <- (
        batch_participant_counts[
            selected_batch
        ] +
        1L
    )
}

participants_with_pairs[
    ,
    batch_number := assigned_batches
]

participants_with_pairs[
    ,
    batch_id := sprintf(
        "%06d",
        batch_number
    )
]


# ==================================== #
#       ATTACH BATCH ASSIGNMENTS        #
# ==================================== #

batch_assignments <- participants_with_pairs[
    ,
    .(
        eid,
        batch_number,
        batch_id,
        code_count,
        estimated_pair_count
    )
]

partitioned_events <- merge(
    first_events,
    batch_assignments[
        ,
        .(
            eid,
            batch_number,
            batch_id
        )
    ],
    by = "eid",
    all = FALSE,
    sort = FALSE
)

expected_partition_rows <- first_events[
    eid %in% batch_assignments$eid,
    .N
]

if (nrow(partitioned_events) != expected_partition_rows) {
    stop(
        paste(
            "Internal validation failed:",
            "the number of partitioned rows does not match",
            "the expected number"
        ),
        call. = FALSE
    )
}




# ==================================== #
#          WRITE BATCH FILES            #
# ==================================== #

batch_files <- character(
    number_of_batches
)

for (
    current_batch_number in
    seq_len(number_of_batches)
) {
    current_batch_id <- sprintf(
        "%06d",
        current_batch_number
    )

    batch_file <- sprintf(
        "participant_batch_%s_%s.csv",
        current_batch_id,
        cohort_code
    )

    batch_files[
        current_batch_number
    ] <- batch_file

    # Explicitly use the loop variable from the calling environment.
    batch_data <- partitioned_events[
        batch_number ==
            current_batch_number
    ]

    if (nrow(batch_data) == 0L) {
        stop(
            paste(
                "Batch",
                current_batch_number,
                "contains no rows"
            ),
            call. = FALSE
        )
    }

    # Remove internal partitioning fields from output.
    batch_data[
        ,
        c(
            "batch_number",
            "batch_id"
        ) := NULL
    ]

    setorder(
        batch_data,
        eid,
        first_date,
        code
    )

    fwrite(
        batch_data,
        file = batch_file,
        sep = ",",
        quote = FALSE,
        na = "NA",
        nThread = 1L
    )

    message(
        "Wrote ",
        batch_file,
        ": ",
        format_count(
            uniqueN(batch_data$eid)
        ),
        " participants; ",
        format_count(
            nrow(batch_data)
        ),
        " rows; estimated ",
        format_count(
            batch_workloads[
                current_batch_number
            ]
        ),
        " pairs"
    )
}



# ==================================== #
#          CREATE MANIFEST              #
# ==================================== #

manifest <- batch_assignments[
    ,
    .(
        participant_count = .N,

        first_event_rows = sum(
            code_count
        ),

        estimated_pair_count = sum(
            estimated_pair_count
        ),

        minimum_codes_per_participant = min(
            code_count
        ),

        median_codes_per_participant = median(
            code_count
        ),

        mean_codes_per_participant = mean(
            code_count
        ),

        maximum_codes_per_participant = max(
            code_count
        ),

        minimum_pairs_per_participant = min(
            estimated_pair_count
        ),

        median_pairs_per_participant = median(
            estimated_pair_count
        ),

        mean_pairs_per_participant = mean(
            estimated_pair_count
        ),

        maximum_pairs_per_participant = max(
            estimated_pair_count
        )
    ),
    by = .(
        batch_number,
        batch_id
    )
]

manifest[
    ,
    cohort_code := cohort_code
]

manifest[
    ,
    batch_file := batch_files[
        batch_number
    ]
]

manifest[
    ,
    target_pairs_per_batch :=
        target_pairs_per_batch
]

manifest[
    ,
    workload_difference_from_target := (
        estimated_pair_count -
        target_pairs_per_batch
    )
]

manifest[
    ,
    workload_ratio_to_target := (
        estimated_pair_count /
        target_pairs_per_batch
    )
]

setcolorder(
    manifest,
    c(
        "cohort_code",
        "batch_number",
        "batch_id",
        "batch_file",
        "target_pairs_per_batch",
        "estimated_pair_count",
        "workload_difference_from_target",
        "workload_ratio_to_target",
        "participant_count",
        "first_event_rows",
        "minimum_codes_per_participant",
        "median_codes_per_participant",
        "mean_codes_per_participant",
        "maximum_codes_per_participant",
        "minimum_pairs_per_participant",
        "median_pairs_per_participant",
        "mean_pairs_per_participant",
        "maximum_pairs_per_participant"
    )
)

setorder(
    manifest,
    batch_number
)

fwrite(
    manifest,
    file = manifest_output,
    sep = ",",
    quote = FALSE,
    na = "NA",
    nThread = 1L
)

# ==================================== #
#       VALIDATE BATCH ASSIGNMENTS      #
# ==================================== #

participant_assignment_check <- batch_assignments[
    ,
    .N,
    by = eid
]

duplicated_assignments <- participant_assignment_check[
    N != 1L
]

if (nrow(duplicated_assignments) > 0L) {
    stop(
        paste(
            nrow(duplicated_assignments),
            "participants were assigned to more or fewer than one batch"
        ),
        call. = FALSE
    )
}

participants_in_batches <- sum(
    manifest$participant_count
)

if (
    participants_in_batches !=
    nrow(participants_with_pairs)
) {
    stop(
        paste(
            "Manifest contains",
            participants_in_batches,
            "participant assignments but expected",
            nrow(participants_with_pairs)
        ),
        call. = FALSE
    )
}

rows_in_batches <- sum(
    manifest$first_event_rows
)

if (
    rows_in_batches !=
    nrow(partitioned_events)
) {
    stop(
        paste(
            "Manifest contains",
            rows_in_batches,
            "first-event rows but expected",
            nrow(partitioned_events)
        ),
        call. = FALSE
    )
}

# ==================================== #
#              VALIDATION               #
# ==================================== #

if (
    sum(manifest$participant_count) !=
    nrow(participants_with_pairs)
) {
    stop(
        paste(
            "Internal validation failed:",
            "manifest participant counts do not match",
            "the assigned participant count"
        ),
        call. = FALSE
    )
}

if (
    sum(manifest$first_event_rows) !=
    nrow(partitioned_events)
) {
    stop(
        paste(
            "Internal validation failed:",
            "manifest event-row counts do not match",
            "the partitioned event rows"
        ),
        call. = FALSE
    )
}

if (
    sum(manifest$estimated_pair_count) !=
    total_estimated_pairs
) {
    stop(
        paste(
            "Internal validation failed:",
            "manifest pair workload does not match",
            "the total estimated pair workload"
        ),
        call. = FALSE
    )
}

assigned_eid_count <- uniqueN(
    partitioned_events$eid
)

if (
    assigned_eid_count !=
    nrow(participants_with_pairs)
) {
    stop(
        paste(
            "Internal validation failed:",
            "one or more eligible participants were lost",
            "or duplicated during partitioning"
        ),
        call. = FALSE
    )
}


# ==================================== #
#               SUMMARY                 #
# ==================================== #

message(
    "Participant batching completed"
)

message(
    "Eligible participants assigned: ",
    format_count(
        nrow(participants_with_pairs)
    )
)

message(
    "Participants excluded because they have fewer than two codes: ",
    format_count(
        nrow(participants_without_pairs)
    )
)

message(
    "Total estimated pairs assigned: ",
    format_count(
        total_estimated_pairs
    )
)

message(
    "Smallest batch workload: ",
    format_count(
        min(manifest$estimated_pair_count)
    )
)

message(
    "Largest batch workload: ",
    format_count(
        max(manifest$estimated_pair_count)
    )
)

message(
    "Mean batch workload: ",
    format_count(
        mean(manifest$estimated_pair_count)
    )
)

message(
    "Wrote batch manifest: ",
    manifest_output
)