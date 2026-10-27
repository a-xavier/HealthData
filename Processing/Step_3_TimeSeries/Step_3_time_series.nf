nextflow.enable.dsl = 2

workflow {

    prepared = PrepareCohort(
        params.diagnosis_events,
        params.participant_metadata,
        params.cohort_code
    )

    participant_partition = PartitionParticipants(
        prepared.first_events,
        params.cohort_code,
        params.target_pair_work 
    )

    participant_batch_files = participant_partition
        .participant_batches
        .flatten()

    temporal_pair_batches = GeneratePairBatch(
        participant_batch_files,
        params.cohort_code
    )

    partial_stats = AggregateTemporalBatch(
        temporal_pair_batches.temporal_pairs,
        params.cohort_code
    )

    temporal_summary = MergeTemporalStats(
        partial_stats.partial_stats.collect(),
        partial_stats.partial_delays.collect(),
        params.cohort_code
    )

   transition_results = CalculateTransitionRisk(
        prepared.first_events,
        prepared.code_summary,
        prepared.cohort_summary,
        temporal_summary.temporal_summary,
        temporal_pair_batches.temporal_pairs.collect(),
        file(
            params.participant_metadata,
            checkIfExists: true
        ),
        params.cohort_code
    )

    transition_batches = SplitTransitionResults(
    transition_results.transition_results,
    params.cohort_code,
    params.km_pairs_per_batch as Integer
)

    transition_batch_files = transition_batches
        .km_pair_batches
        .flatten()

    km_batch_results = EstimateKaplanMeierBatch(
        transition_batch_files,
        prepared.first_events,
        prepared.cohort_summary,
        params.cohort_code
    )

    final_km_results = MergeKaplanMeierResults(
        km_batch_results.km_results.collect(),
        params.cohort_code
    )
}

process PrepareCohort {

    tag "${cohort_code}"

    input:
        path diagnosis_events
        path participant_metadata
        val cohort_code

    output:
        path "clean_events_${cohort_code}.csv",
            emit: clean_events

        path "first_events_${cohort_code}.csv",
            emit: first_events

        path "code_summary_${cohort_code}.csv",
            emit: code_summary

        path "cohort_summary_${cohort_code}.csv",
            emit: cohort_summary

    script:
    """
    prepare_cohort.R \
        "${diagnosis_events}" \
        "${participant_metadata}" \
        "${cohort_code}" \
        "clean_events_${cohort_code}.csv" \
        "first_events_${cohort_code}.csv" \
        "code_summary_${cohort_code}.csv" \
        "cohort_summary_${cohort_code}.csv"
    """
}

process PartitionParticipants {

    tag "${cohort_code}"

    input:
        path first_events
        val cohort_code
        val target_pairs_per_batch

    output:
        path "participant_batch_??????_${cohort_code}.csv",
            emit: participant_batches

        path "participant_batch_manifest_${cohort_code}.csv",
            emit: batch_manifest

    script:
    """
    partition_participants.R \
        "${first_events}" \
        "${cohort_code}" \
        "${target_pairs_per_batch}" \
        "participant_batch_manifest_${cohort_code}.csv"
    """
}

process GeneratePairBatch {

    tag "${participant_batch.simpleName}"

    maxForks 8
    cpus 1

    input:
        path participant_batch
        val cohort_code

    output:
        path "temporal_pairs_${participant_batch.simpleName}.csv",
            emit: temporal_pairs

    script:
        def output_file =
            "temporal_pairs_${participant_batch.simpleName}.csv"

    """
    generate_pair_batch.R \
        "${participant_batch}" \
        "${cohort_code}" \
        "${output_file}"
    """
}

process AggregateTemporalBatch {

    tag "${temporal_pair_batch.simpleName}"

    maxForks 8
    cpus 1

    input:
        path temporal_pair_batch
        val cohort_code

    output:
        path "partial_stats_${temporal_pair_batch.simpleName}.csv",
            emit: partial_stats

        path "partial_delays_${temporal_pair_batch.simpleName}.csv",
            emit: partial_delays

    script:
        def stats_output =
            "partial_stats_${temporal_pair_batch.simpleName}.csv"

        def delays_output =
            "partial_delays_${temporal_pair_batch.simpleName}.csv"

    """
    AggregateTemporalBatch.R \
        "${temporal_pair_batch}" \
        "${cohort_code}" \
        "${stats_output}" \
        "${delays_output}"
    """
}

process MergeTemporalStats {

    tag "${cohort_code}"

    input:
        path partial_stats_files
        path partial_delay_files
        val cohort_code

    output:
        path "temporal_pair_summary_${cohort_code}.csv",
            emit: temporal_summary

    script:
    """
    MergeTemporalStats.R \
        "${cohort_code}" \
        "temporal_pair_summary_${cohort_code}.csv"
    """
}

process CalculateTransitionRisk {

    tag "${cohort_code}"

    publishDir ".", mode: 'copy'

    input:
        path first_events
        path code_summary
        path cohort_summary
        path temporal_summary
        path temporal_pair_files
        path participant_metadata
        val cohort_code

    output:
        path "temporal_transition_results_${cohort_code}.csv",
            emit: transition_results

    script:
    """
    CalculateTransitionRisk.R \
        "${first_events}" \
        "${code_summary}" \
        "${cohort_summary}" \
        "${temporal_summary}" \
        "${participant_metadata}" \
        "${cohort_code}" \
        "temporal_transition_results_${cohort_code}.csv"
    """
}

process SplitTransitionResults {

    tag "${cohort_code}"

    input:
        path transition_results
        val cohort_code
        val km_pairs_per_batch

    output:
        path "km_pair_batch_??????_${cohort_code}.csv",
            emit: km_pair_batches

    shell:
    '''
    awk \
        -v batch_size="!{km_pairs_per_batch}" \
        -v cohort="!{cohort_code}" \
        '
        BEGIN {
            if (batch_size < 1) {
                print "batch_size must be a positive integer" > "/dev/stderr"
                exit 1
            }
        }

        NR == 1 {
            header = $0
            next
        }

        {
            batch_number = int((NR - 2) / batch_size) + 1
            output_file = sprintf("km_pair_batch_%06d_%s.csv", batch_number, cohort)

            if (!(output_file in initialized)) {
                print header > output_file
                initialized[output_file] = 1
            }

            print $0 >> output_file
        }

        END {
            if (NR < 2) {
                print "The transition results file contains no data rows" > "/dev/stderr"
                exit 1
            }
        }
        ' \
        "!{transition_results}"
    '''
}

process EstimateKaplanMeierBatch {

    tag "${transition_batch.simpleName}"

    maxForks 8
    cpus 1

    input:
        path transition_batch
        path first_events
        path cohort_summary
        val cohort_code

    output:
        path "km_results_${transition_batch.simpleName}.csv",
            emit: km_results

    script:
        def output_file =
            "km_results_${transition_batch.simpleName}.csv"

    """
    estimate_kaplan_meier_batch.R \
        "${transition_batch}" \
        "${first_events}" \
        "${cohort_summary}" \
        "${cohort_code}" \
        "${output_file}"
    """
}

process MergeKaplanMeierResults {

    tag "${cohort_code}"

    publishDir ".", mode: 'copy'

    input:
        path km_result_files
        val cohort_code

    output:
        path "temporal_transition_km_${cohort_code}.csv",
            emit: km_transition_results

    shell:
    '''
    output_file="temporal_transition_km_!{cohort_code}.csv"
    first_file=true

    for result_file in km_results_*.csv; do
        if [ "$first_file" = true ]; then
            cat "$result_file" > "$output_file"
            first_file=false
        else
            tail -n +2 "$result_file" >> "$output_file"
        fi
    done
    '''
}