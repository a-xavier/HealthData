// Pipeline is simple 
// 2 inputs files to use in command line 
// - all_codes.txt which contains all codes 
// - the final_df.csv for both UKB and GS 
// also threshold to use for flagging sex imbalance (default 90)
// also cohort code to identify the cohort being processed and name output 

// Flow is:
// 1 - for each code check the all_levels columns and the is_male columns
// 2 - grep the code in all_levels and transform into 0 1 
// 3 -  Do a table phenotype 0 1 vs sex 0 1 and see see if imbalance, 
// if more then 90 percent is one sex, we flag as either F or M otherwise we leave it as NA



process Sex_Codes_Calculation {
    
    input:
        val code
        path final_df_file
        val threshold
        val cohort_code
    
    output:
        path "sex_codes_output_${code}_${cohort_code}.txt"

    script:

        """
        # Your processing code here
        process.py \
            --code ${code} \
            --final_df_file ${final_df_file} \
            --threshold ${threshold} \
            --cohort_code ${cohort_code}
        """
}

process Merge_Sex_Codes {

    publishDir "./" , mode: 'copy'

    input:
        path result_files
        val cohort_code

    output:
        path "combined_sex_${cohort_code}.tsv"

    script:
    """
    cat ${result_files.join(' ')} > combined_sex_${cohort_code}.tsv
    """
}

workflow {

    all_codes_ch = channel
            .fromPath(params.all_codes_file)
            .splitText()
            .map { it -> it.trim() }
            .filter { it -> it }

    // all_codes_ch = channel.fromList(["Chapter_I", "A00", "G35", "caca"])

    final_df_file_ch = channel.value(file(params.final_df_file))
    threshold_ch    = channel.value(params.threshold)
    cohort_code_ch  = channel.value(params.cohort_code)

    result_files = Sex_Codes_Calculation(
        all_codes_ch,
        final_df_file_ch,
        threshold_ch,
        cohort_code_ch
    )

    Merge_Sex_Codes(
        result_files.collect(),
        cohort_code_ch
    )
}

