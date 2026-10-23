// This one is to test for each common covariate 
// age sex smokever bmi drinking moderate_activ vigorous_activity (also family for GS) PC 1 to 5

process ProcessSingleCode {

    input:
        tuple val(code), val(sex_flag), val(sex_prop), val(phenotype_count)
        val cohort_code
        path final_df_file

    output:
        path "processed_${code}_${cohort_code}.csv"

    script:
    """
    process.R $code $sex_flag $cohort_code $final_df_file
    """
}

process MergeResults {

    tag "${cohort_code}"

    publishDir ".", mode: 'copy'

    input:
        path processed_files
        val cohort_code

    output:
        path "merged_results_${cohort_code}.csv"

    script:
    """
    awk 'FNR == 1 && NR != 1 { next } { print }' \
        ${processed_files} \
        > merged_results_${cohort_code}.csv
    """
}




workflow {

    // read params.sex_code_file generated before 
    // split each line in 4 parts: code, sex_flag,  sex_prop, phenotype_count
    channel.fromPath(params.sex_code_file)
        .splitText()
        .map { line -> 
            def parts = line.split("\t")
            tuple(parts[0], parts[1], parts[2], parts[3])
        }
        .set { sex_code_channel }


    // Filter where parts[3] (phenotype_count) is greater than 5
    sex_code_channel = sex_code_channel.filter { tuple -> 
        tuple[3].toInteger() > 5
    }


    // Cohort code
    cohort_code_channel = channel.value(params.cohort_code)

    // Process each code 
    ProcessSingleCode(
        sex_code_channel,
        cohort_code_channel,
        file(params.final_df_file, checkIfExists: true)
    )

    // Merge all processed results for the cohort
    MergeResults(
        ProcessSingleCode.out.collect(),
        cohort_code_channel
    )
}