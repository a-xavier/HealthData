#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
    library(survival)
})


# ==================================== #
#              ARGUMENTS                #
# ==================================== #

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5L) {
    stop(
        paste(
            "Expected 5 arguments:",
            "<transition_batch>",
            "<first_events>",
            "<cohort_summary>",
            "<cohort_code>",
            "<output_file>"
        ),
        call. = FALSE
    )
}

transition_batch_file <- args[[1]]
first_events_file <- args[[2]]
cohort_summary_file <- args[[3]]
cohort_code <- trimws(args[[4]])
output_file <- args[[5]]

input_files <- c(
    transition_batch_file,
    first_events_file,
    cohort_summary_file
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
#        REPORTING THRESHOLDS           #
# ==================================== #

# These thresholds determine whether a horizon-specific estimate is
# considered sufficiently supported for reporting.
#
# The estimate is still retained when unsupported, but the
# km_supported_* field is FALSE and km_status_* explains why.

minimum_eligible_count <- 50L
minimum_events_by_horizon <- 10L
minimum_number_at_risk <- 20L

confidence_level <- 0.95


# ==================================== #
#            ANALYSIS LABELS            #
# ==================================== #

# Death is not currently modelled as a competing event.
#
# Participants without the target diagnosis are administratively
# censored at the cohort-wide observation end.

km_analysis_type <- paste0(
    "Kaplan-Meier transition-free survival ",
    "with administrative censoring"
)

km_competing_death_handled <- FALSE


# ==================================== #
#              HORIZONS                 #
# ==================================== #

horizons <- c(
    "30_days" = 30,
    "90_days" = 90,
    "1_year" = 365,
    "2_years" = 731,
    "5_years" = 1826,
    "10_years" = 3653
)


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


# ------------------------------------ #
# Exact binomial interval               #
# ------------------------------------ #

exact_binomial_interval <- function(
    events,
    total,
    confidence_level = 0.95
) {
    result <- list(
        estimate = NA_real_,
        lower = NA_real_,
        upper = NA_real_
    )

    if (
        is.na(events) ||
        is.na(total) ||
        total <= 0L ||
        events < 0L ||
        events > total
    ) {
        return(result)
    }

    alpha <- 1 - confidence_level
    non_events <- total - events

    result$estimate <- events / total

    result$lower <- if (events == 0L) {
        0
    } else {
        qbeta(
            alpha / 2,
            events,
            non_events + 1L
        )
    }

    result$upper <- if (non_events == 0L) {
        1
    } else {
        qbeta(
            1 - alpha / 2,
            events + 1L,
            non_events
        )
    }

    result
}


# ------------------------------------ #
# Empty Kaplan-Meier result             #
# ------------------------------------ #

empty_km_result <- function() {

    result <- list(
        km_fit_ok = FALSE,
        km_fit_message = NA_character_,
        km_analysis_type = km_analysis_type,
        km_competing_death_handled =
            km_competing_death_handled,

        km_eligible_count = NA_integer_,
        km_event_count = NA_integer_,
        km_censored_count = NA_integer_,
        km_max_followup_days = NA_real_,

        km_median_estimable = FALSE,
        km_median_transition_days = NA_real_,
        km_median_transition_ci_low = NA_real_,
        km_median_transition_ci_high = NA_real_
    )

    for (horizon_name in names(horizons)) {

        result[[
            paste0(
                "km_probability_",
                horizon_name
            )
        ]] <- NA_real_

        result[[
            paste0(
                "km_probability_",
                horizon_name,
                "_ci_low"
            )
        ]] <- NA_real_

        result[[
            paste0(
                "km_probability_",
                horizon_name,
                "_ci_high"
            )
        ]] <- NA_real_

        result[[
            paste0(
                "km_number_at_risk_",
                horizon_name
            )
        ]] <- NA_integer_

        result[[
            paste0(
                "km_events_by_",
                horizon_name
            )
        ]] <- NA_integer_

        result[[
            paste0(
                "km_evaluable_at_",
                horizon_name
            )
        ]] <- NA_integer_

        result[[
            paste0(
                "km_supported_",
                horizon_name
            )
        ]] <- FALSE

        result[[
            paste0(
                "km_status_",
                horizon_name
            )
        ]] <- "not_calculated"

        result[[
            paste0(
                "km_ci_method_",
                horizon_name
            )
        ]] <- NA_character_
    }

    result
}


# ------------------------------------ #
# Extract median transition time        #
# ------------------------------------ #

extract_median_survival <- function(fit) {

    result <- list(
        estimable = FALSE,
        median = NA_real_,
        lower = NA_real_,
        upper = NA_real_
    )

    fit_summary <- tryCatch(
        summary(fit),
        error = function(e) {
            NULL
        }
    )

    if (
        is.null(fit_summary) ||
        is.null(fit_summary$table)
    ) {
        return(result)
    }

    fit_table <- fit_summary$table

    if (is.matrix(fit_table)) {
        fit_table <- fit_table[1L, ]
    }

    table_names <- names(fit_table)

    if ("median" %in% table_names) {
        result$median <- as.numeric(
            fit_table[["median"]]
        )
    }

    if ("0.95LCL" %in% table_names) {
        result$lower <- as.numeric(
            fit_table[["0.95LCL"]]
        )
    }

    if ("0.95UCL" %in% table_names) {
        result$upper <- as.numeric(
            fit_table[["0.95UCL"]]
        )
    }

    result$estimable <- (
        length(result$median) == 1L &&
        !is.na(result$median) &&
        is.finite(result$median)
    )

    if (!result$estimable) {
        result$median <- NA_real_
        result$lower <- NA_real_
        result$upper <- NA_real_
    }

    result
}


# ------------------------------------ #
# Extract one horizon                   #
# ------------------------------------ #

extract_survival_at_time <- function(
    fit,
    risk_set,
    horizon_days,
    maximum_followup,
    eligible_count
) {
    result <- list(
        probability = NA_real_,
        probability_ci_low = NA_real_,
        probability_ci_high = NA_real_,
        number_at_risk = NA_integer_,
        events_by_horizon = NA_integer_,
        evaluable_at_horizon = NA_integer_,
        supported = FALSE,
        status = "not_calculated",
        ci_method = NA_character_
    )

    if (
        is.na(maximum_followup) ||
        horizon_days > maximum_followup
    ) {
        result$status <- "horizon_beyond_followup"
        return(result)
    }

    events_by_horizon <- risk_set[
        event == 1L &
        time_days <= horizon_days,
        .N
    ]

    # A participant is evaluable at the horizon if:
    #
    # 1. the transition occurred by the horizon; or
    # 2. the participant remained observed through the horizon.
    #
    # Participants censored before the horizon without an event are
    # not counted as having a known horizon-specific outcome.

    evaluable_at_horizon <- risk_set[
        (
            event == 1L &
            time_days <= horizon_days
        ) |
        time_days >= horizon_days,
        .N
    ]

    fit_summary <- tryCatch(
        summary(
            fit,
            times = horizon_days,
            extend = FALSE
        ),
        error = function(e) {
            NULL
        }
    )

    if (
        is.null(fit_summary) ||
        length(fit_summary$surv) == 0L
    ) {
        result$status <- "survival_estimate_unavailable"
        return(result)
    }

    survival_probability <- as.numeric(
        fit_summary$surv[[1]]
    )

    survival_ci_low <- as.numeric(
        fit_summary$lower[[1]]
    )

    survival_ci_high <- as.numeric(
        fit_summary$upper[[1]]
    )

    number_at_risk <- as.integer(
        fit_summary$n.risk[[1]]
    )

    result$number_at_risk <- number_at_risk
    result$events_by_horizon <- events_by_horizon
    result$evaluable_at_horizon <-
        evaluable_at_horizon

    if (
        is.na(number_at_risk) ||
        number_at_risk <= 0L
    ) {
        result$status <- "empty_risk_set_at_horizon"
        return(result)
    }

    probability <- 1 - survival_probability

    probability_ci_low <- (
        1 - survival_ci_high
    )

    probability_ci_high <- (
        1 - survival_ci_low
    )

    ci_method <- "kaplan_meier_log_log"

    # The transformed Greenwood interval can be undefined or
    # degenerate at a boundary. When there are zero events by the
    # horizon, use an exact interval based on participants with a
    # known outcome through the horizon.
    #
    # This is explicitly labelled as a boundary fallback rather than
    # a standard Kaplan-Meier confidence interval.

    if (
        events_by_horizon == 0L &&
        evaluable_at_horizon > 0L
    ) {
        boundary_interval <- exact_binomial_interval(
            events = 0L,
            total = evaluable_at_horizon,
            confidence_level = confidence_level
        )

        probability <- 0
        probability_ci_low <-
            boundary_interval$lower
        probability_ci_high <-
            boundary_interval$upper

        ci_method <-
            "exact_binomial_boundary_fallback"
    }

    # Apply the same fallback if all evaluable participants have
    # transitioned by the horizon.

    if (
        evaluable_at_horizon > 0L &&
        events_by_horizon ==
            evaluable_at_horizon
    ) {
        boundary_interval <- exact_binomial_interval(
            events = events_by_horizon,
            total = evaluable_at_horizon,
            confidence_level = confidence_level
        )

        probability <- 1
        probability_ci_low <-
            boundary_interval$lower
        probability_ci_high <-
            boundary_interval$upper

        ci_method <-
            "exact_binomial_boundary_fallback"
    }

    # Clamp negligible floating-point excursions.

    probability <- min(
        max(probability, 0),
        1
    )

    if (!is.na(probability_ci_low)) {
        probability_ci_low <- min(
            max(probability_ci_low, 0),
            1
        )
    }

    if (!is.na(probability_ci_high)) {
        probability_ci_high <- min(
            max(probability_ci_high, 0),
            1
        )
    }

    result$probability <- probability
    result$probability_ci_low <-
        probability_ci_low
    result$probability_ci_high <-
        probability_ci_high
    result$ci_method <- ci_method

    # Determine whether the estimate is sufficiently supported for
    # reporting.

    if (
        eligible_count <
        minimum_eligible_count
    ) {
        result$status <- "insufficient_eligible"

    } else if (
        events_by_horizon <
        minimum_events_by_horizon
    ) {
        result$status <-
            "insufficient_events_by_horizon"

    } else if (
        number_at_risk <
        minimum_number_at_risk
    ) {
        result$status <-
            "insufficient_number_at_risk"

    } else if (
        is.na(result$probability)
    ) {
        result$status <-
            "probability_unavailable"

    } else if (
        is.na(result$probability_ci_low) ||
        is.na(result$probability_ci_high)
    ) {
        result$status <-
            "confidence_interval_unavailable"

    } else {
        result$status <- "supported"
        result$supported <- TRUE
    }

    result
}


# ------------------------------------ #
# Validate one KM horizon               #
# ------------------------------------ #

validate_horizon_result <- function(
    probability,
    lower,
    upper,
    number_at_risk,
    supported,
    status
) {
    if (
        is.na(status) ||
        status == ""
    ) {
        return(FALSE)
    }

    if (!is.na(probability)) {
        if (
            probability < 0 ||
            probability > 1
        ) {
            return(FALSE)
        }
    }

    if (!is.na(lower)) {
        if (
            lower < 0 ||
            lower > 1
        ) {
            return(FALSE)
        }
    }

    if (!is.na(upper)) {
        if (
            upper < 0 ||
            upper > 1
        ) {
            return(FALSE)
        }
    }

    if (
        !is.na(probability) &&
        !is.na(lower) &&
        lower > probability
    ) {
        return(FALSE)
    }

    if (
        !is.na(probability) &&
        !is.na(upper) &&
        upper < probability
    ) {
        return(FALSE)
    }

    if (
        !is.na(lower) &&
        !is.na(upper) &&
        lower > upper
    ) {
        return(FALSE)
    }

    if (isTRUE(supported)) {
        if (
            is.na(probability) ||
            is.na(lower) ||
            is.na(upper) ||
            is.na(number_at_risk) ||
            number_at_risk <
                minimum_number_at_risk ||
            status != "supported"
        ) {
            return(FALSE)
        }
    }

    TRUE
}


# ==================================== #
#          LOAD COHORT SUMMARY          #
# ==================================== #

cohort_summary <- fread(
    cohort_summary_file,
    sep = ",",
    colClasses = c(
        cohort_code = "character",
        observation_start = "character",
        observation_end = "character"
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

required_cohort_columns <- c(
    "cohort_code",
    "observation_start",
    "observation_end"
)

missing_columns <- setdiff(
    required_cohort_columns,
    names(cohort_summary)
)

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "Cohort summary is missing:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

if (nrow(cohort_summary) != 1L) {
    stop(
        "Cohort summary must contain exactly one row",
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
        "Cohort summary does not match cohort_code",
        call. = FALSE
    )
}

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


# ==================================== #
#          LOAD FIRST EVENTS            #
# ==================================== #

first_events <- fread(
    first_events_file,
    sep = ",",
    colClasses = c(
        eid = "character",
        code = "character",
        first_date = "character",
        is_male = "integer",
        age_at_first_diagnosis = "numeric"
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

required_event_columns <- c(
    "eid",
    "code",
    "first_date",
    "is_male",
    "age_at_first_diagnosis"
)

missing_columns <- setdiff(
    required_event_columns,
    names(first_events)
)

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "First-events file is missing:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

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

duplicate_events <- first_events[
    ,
    .N,
    by = .(
        eid,
        code
    )
][
    N > 1L
]

if (nrow(duplicate_events) > 0L) {
    stop(
        "first_events contains duplicate eid-code combinations",
        call. = FALSE
    )
}

if (any(
    first_events$first_date >
        observation_end
)) {
    stop(
        paste(
            "first_events contains diagnosis dates after",
            "the cohort observation end"
        ),
        call. = FALSE
    )
}

setkey(
    first_events,
    code,
    eid
)


# ==================================== #
#       LOAD TRANSITION BATCH           #
# ==================================== #

transition_batch <- fread(
    transition_batch_file,
    sep = ",",
    na.strings = c(
        "NA",
        ""
    ),
    showProgress = FALSE,
    nThread = 1L,
    check.names = FALSE
)

required_transition_columns <- c(
    "cohort_code",
    "source_code",
    "target_code"
)

missing_columns <- setdiff(
    required_transition_columns,
    names(transition_batch)
)

if (length(missing_columns) > 0L) {
    stop(
        paste(
            "Transition batch is missing:",
            paste(
                missing_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}

transition_batch[
    ,
    cohort_code := trimws(cohort_code)
]

transition_batch[
    ,
    source_code := toupper(
        trimws(source_code)
    )
]

transition_batch[
    ,
    target_code := toupper(
        trimws(target_code)
    )
]

if (any(
    transition_batch$cohort_code !=
        cohort_code
)) {
    stop(
        "Transition batch contains a different cohort_code",
        call. = FALSE
    )
}

if (any(
    transition_batch$source_code ==
        transition_batch$target_code
)) {
    stop(
        "Self-transitions were found",
        call. = FALSE
    )
}

if (anyDuplicated(
    transition_batch,
    by = c(
        "cohort_code",
        "source_code",
        "target_code"
    )
)) {
    stop(
        "Transition batch contains duplicate directed pairs",
        call. = FALSE
    )
}

message(
    "Directed transitions in batch: ",
    format_count(
        nrow(transition_batch)
    )
)


# ==================================== #
#       KAPLAN-MEIER FUNCTION           #
# ==================================== #

estimate_transition_km <- function(
    source_code,
    target_code
) {
    result <- empty_km_result()

    source_events <- first_events[
        .(source_code),
        on = .(
            code
        ),
        nomatch = 0L,
        .(
            eid,
            source_date =
                first_date,
            age_at_source =
                age_at_first_diagnosis,
            is_male
        )
    ]

    if (nrow(source_events) == 0L) {
        result$km_fit_message <-
            "source_code_not_found"
        return(result)
    }

    target_events <- first_events[
        .(target_code),
        on = .(
            code
        ),
        nomatch = 0L,
        .(
            eid,
            target_date =
                first_date
        )
    ]

    risk_set <- merge(
        source_events,
        target_events,
        by = "eid",
        all.x = TRUE,
        sort = FALSE
    )

    # Exclude participants whose target diagnosis was already
    # recorded before or on the source diagnosis date.
    risk_set <- risk_set[
        is.na(target_date) |
        target_date > source_date
    ]

    if (nrow(risk_set) == 0L) {
        result$km_fit_message <-
            "no_eligible_participants"
        return(result)
    }

    risk_set[
        ,
        event := as.integer(
            !is.na(target_date) &
            target_date <= observation_end
        )
    ]

    risk_set[
        ,
        end_date := observation_end
    ]

    risk_set[
        event == 1L,
        end_date := target_date
    ]

    risk_set[
        ,
        time_days := as.numeric(
            end_date -
            source_date
        )
    ]

    # Zero or negative follow-up cannot contribute to a conventional
    # right-sided survival analysis.
    risk_set <- risk_set[
        !is.na(time_days) &
        time_days > 0
    ]

    if (nrow(risk_set) == 0L) {
        result$km_fit_message <-
            "no_positive_followup"
        return(result)
    }

    result$km_eligible_count <- nrow(
        risk_set
    )

    result$km_event_count <- sum(
        risk_set$event
    )

    result$km_censored_count <- sum(
        risk_set$event == 0L
    )

    result$km_max_followup_days <- max(
        risk_set$time_days
    )

    fit_error_message <- NA_character_

    fit <- tryCatch(
        survfit(
            Surv(
                time_days,
                event
            ) ~ 1,
            data = risk_set,
            conf.int = confidence_level,
            conf.type = "log-log"
        ),
        error = function(e) {
            fit_error_message <<-
                conditionMessage(e)
            NULL
        }
    )

    if (is.null(fit)) {
        result$km_fit_message <- paste0(
            "survfit_failed: ",
            fit_error_message
        )
        return(result)
    }

    result$km_fit_ok <- TRUE
    result$km_fit_message <- "ok"

    median_result <- extract_median_survival(
        fit
    )

    result$km_median_estimable <-
        median_result$estimable

    result$km_median_transition_days <-
        median_result$median

    result$km_median_transition_ci_low <-
        median_result$lower

    result$km_median_transition_ci_high <-
        median_result$upper

    for (horizon_name in names(horizons)) {

        horizon_result <- extract_survival_at_time(
            fit = fit,
            risk_set = risk_set,
            horizon_days =
                horizons[[horizon_name]],
            maximum_followup =
                result$km_max_followup_days,
            eligible_count =
                result$km_eligible_count
        )

        result[[
            paste0(
                "km_probability_",
                horizon_name
            )
        ]] <- horizon_result$probability

        result[[
            paste0(
                "km_probability_",
                horizon_name,
                "_ci_low"
            )
        ]] <- horizon_result$
            probability_ci_low

        result[[
            paste0(
                "km_probability_",
                horizon_name,
                "_ci_high"
            )
        ]] <- horizon_result$
            probability_ci_high

        result[[
            paste0(
                "km_number_at_risk_",
                horizon_name
            )
        ]] <- horizon_result$
            number_at_risk

        result[[
            paste0(
                "km_events_by_",
                horizon_name
            )
        ]] <- horizon_result$
            events_by_horizon

        result[[
            paste0(
                "km_evaluable_at_",
                horizon_name
            )
        ]] <- horizon_result$
            evaluable_at_horizon

        result[[
            paste0(
                "km_supported_",
                horizon_name
            )
        ]] <- horizon_result$supported

        result[[
            paste0(
                "km_status_",
                horizon_name
            )
        ]] <- horizon_result$status

        result[[
            paste0(
                "km_ci_method_",
                horizon_name
            )
        ]] <- horizon_result$ci_method
    }

    result
}


# ==================================== #
#       PROCESS TRANSITION BATCH        #
# ==================================== #

km_results <- vector(
    mode = "list",
    length = nrow(transition_batch)
)

for (
    transition_index in
    seq_len(nrow(transition_batch))
) {
    source_code <- transition_batch$
        source_code[
            transition_index
        ]

    target_code <- transition_batch$
        target_code[
            transition_index
        ]

    if (
        transition_index == 1L ||
        transition_index %% 100L == 0L ||
        transition_index ==
            nrow(transition_batch)
    ) {
        message(
            "Processing transition ",
            transition_index,
            " of ",
            nrow(transition_batch),
            ": ",
            source_code,
            " -> ",
            target_code
        )
    }

    km_results[[
        transition_index
    ]] <- as.data.table(
        estimate_transition_km(
            source_code = source_code,
            target_code = target_code
        )
    )
}

km_results <- rbindlist(
    km_results,
    use.names = TRUE,
    fill = TRUE
)

result_table <- cbind(
    transition_batch,
    km_results
)


# ==================================== #
#       ADD REPORTING THRESHOLDS        #
# ==================================== #

result_table[
    ,
    km_minimum_eligible_required :=
        minimum_eligible_count
]

result_table[
    ,
    km_minimum_events_required :=
        minimum_events_by_horizon
]

result_table[
    ,
    km_minimum_at_risk_required :=
        minimum_number_at_risk
]


# ==================================== #
#              VALIDATION               #
# ==================================== #

if (nrow(result_table) != nrow(transition_batch)) {
    stop(
        "Kaplan-Meier output row count does not match input",
        call. = FALSE
    )
}

if (any(
    result_table$km_event_count >
        result_table$km_eligible_count,
    na.rm = TRUE
)) {
    stop(
        "KM event count exceeds eligible count",
        call. = FALSE
    )
}

if (any(
    result_table$km_censored_count +
        result_table$km_event_count !=
        result_table$km_eligible_count,
    na.rm = TRUE
)) {
    stop(
        paste(
            "KM event and censoring counts do not sum",
            "to the eligible population"
        ),
        call. = FALSE
    )
}

for (horizon_name in names(horizons)) {

    probability_column <- paste0(
        "km_probability_",
        horizon_name
    )

    lower_column <- paste0(
        probability_column,
        "_ci_low"
    )

    upper_column <- paste0(
        probability_column,
        "_ci_high"
    )

    at_risk_column <- paste0(
        "km_number_at_risk_",
        horizon_name
    )

    supported_column <- paste0(
        "km_supported_",
        horizon_name
    )

    status_column <- paste0(
        "km_status_",
        horizon_name
    )

    validation_result <- mapply(
        validate_horizon_result,

        probability =
            result_table[[
                probability_column
            ]],

        lower =
            result_table[[
                lower_column
            ]],

        upper =
            result_table[[
                upper_column
            ]],

        number_at_risk =
            result_table[[
                at_risk_column
            ]],

        supported =
            result_table[[
                supported_column
            ]],

        status =
            result_table[[
                status_column
            ]]
    )

    invalid_rows <- which(
        !validation_result
    )

    if (length(invalid_rows) > 0L) {

        diagnostic_columns <- c(
            "source_code",
            "target_code",
            probability_column,
            lower_column,
            upper_column,
            at_risk_column,
            supported_column,
            status_column
        )

        fwrite(
            result_table[
                invalid_rows,
                ..diagnostic_columns
            ],
            file = paste0(
                "invalid_km_",
                horizon_name,
                ".csv"
            ),
            sep = ",",
            quote = FALSE,
            na = "NA"
        )

        stop(
            paste(
                length(invalid_rows),
                "invalid Kaplan-Meier rows found for",
                horizon_name
            ),
            call. = FALSE
        )
    }
}


# ==================================== #
#              WRITE OUTPUT             #
# ==================================== #

fwrite(
    result_table,
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
    "Kaplan-Meier transitions processed: ",
    format_count(
        nrow(result_table)
    )
)

message(
    "Successful Kaplan-Meier fits: ",
    format_count(
        result_table[
            km_fit_ok == TRUE,
            .N
        ]
    )
)

message(
    "Transitions with at least one event: ",
    format_count(
        result_table[
            km_event_count > 0L,
            .N
        ]
    )
)

message(
    "Transitions with an estimable median transition time: ",
    format_count(
        result_table[
            km_median_estimable == TRUE,
            .N
        ]
    )
)

message(
    "Supported 30-day estimates: ",
    format_count(
        result_table[
            km_supported_30_days == TRUE,
            .N
        ]
    )
)

message(
    "Supported one-year estimates: ",
    format_count(
        result_table[
            km_supported_1_year == TRUE,
            .N
        ]
    )
)

message(
    "Supported five-year estimates: ",
    format_count(
        result_table[
            km_supported_5_years == TRUE,
            .N
        ]
    )
)

message(
    "Supported ten-year estimates: ",
    format_count(
        result_table[
            km_supported_10_years == TRUE,
            .N
        ]
    )
)

message(
    "Wrote Kaplan-Meier results: ",
    output_file
)