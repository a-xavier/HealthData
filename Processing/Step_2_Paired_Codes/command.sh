nextflow run Step_2_Paired_codes.nf -resume \
--icd10_pair_code_file /Users/alex/Downloads/slice_observed_pairs.txt \
--cohort_code UKB \
--final_df_file /Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/UKB/Final_Data_UKB.csv



Rscript /Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/Processing/Step_2_Paired_Codes/bin/process.R \
Chapter_I \
Chapter_II \
UKB \
/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/UKB/Final_Data_UKB.csv