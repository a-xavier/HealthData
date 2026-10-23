// this pipeline is be for processing icd10 paired codes
nextflow.enable.dsl=2

process ReadAndProcessPairs {

    tag "min_count=${min_both_count}, batch_size=${batch_size}"

    input:
        path icd10_pair_code_file
        val min_both_count
        val batch_size

    output:
        path "pair_batch_*", emit: pair_batches

    script:
    """
    awk -F '\\t' \
        -v min_count="${min_both_count}" \
        '\$3 >= min_count' \
        "${icd10_pair_code_file}" \
    | split \
        -l "${batch_size}" \
        -a 6 \
        - \
        pair_batch_
    """
}

process ProcessPairs {

    maxForks 5

    tag "${pair_batch.simpleName} - ${cohort_code}"

    input:
        path pair_batch
        val cohort_code
        path final_df

    output:
        path "results_${pair_batch.simpleName}_${cohort_code}.csv"

    script:
        def output_file = "results_${pair_batch.simpleName}_${cohort_code}.csv"

    """
    process.R \
        "${pair_batch}" \
        "${cohort_code}" \
        "${final_df}" \
        "${output_file}"
    """
}


process MergeReports {

    tag "${cohort_code}"

    publishDir ".", mode: 'copy'

    input:
        path report_files
        val cohort_code

    output:
        path "paired_codes_reports_${cohort_code}.csv"

    script:
    """
    output_file="paired_codes_reports_${cohort_code}.csv"

    first_file=1

    for report_file in results_*_${cohort_code}.csv; do
        if [ "\$first_file" -eq 1 ]; then
            cat "\$report_file" > "\$output_file"
            first_file=0
        else
            tail -n +2 "\$report_file" >> "\$output_file"
        fi
    done
    """
}

workflow {

    pair_batches_output = ReadAndProcessPairs(
        file(
            params.icd10_pair_code_file,
            checkIfExists: true
        ),
        params.min_both_count as Integer,
        params.batch_size as Integer
    )

    // Convert the list of batch files into one emission per batch file
    pair_batches_ch = pair_batches_output.pair_batches.flatten()

    report_files_ch = ProcessPairs(
        pair_batches_ch,
        params.cohort_code,
        file(
            params.final_df_file,
            checkIfExists: true
        )
    )

    MergeReports(
        report_files_ch.collect(),
        params.cohort_code
    )
}
