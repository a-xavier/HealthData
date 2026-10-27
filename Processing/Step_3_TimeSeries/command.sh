nextflow run Step_3_time_series.nf -resume \
    --diagnosis_events /Users/alex/Documents/ICD10_Project_Official/Generation_Scotland/Codes_with_date/GS_code_date.csv \
    --participant_metadata /Users/alex/Documents/ICD10_Project_Official/Generation_Scotland/Codes_with_date/GS_metadata_age_dob.csv \
    --cohort_code GS \
    --cohort_size 20019

nextflow run Step_3_time_series.nf -resume \
    --diagnosis_events /Users/alex/Documents/ICD10_Project_Official/UKB/Codes_with_date/UKB_code_date.csv \
    --participant_metadata /Users/alex/Documents/ICD10_Project_Official/UKB/Codes_with_date/UKB_Metadata_age_dob.csv \
    --cohort_code UKB \
    --cohort_size 500000
