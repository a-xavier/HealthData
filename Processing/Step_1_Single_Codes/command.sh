nextflow run Step_1_Single_Codes.nf \
--sex_code_file /Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/Processing/Step_0_Sex_Codes/combined_sex_UKB.tsv \
--final_df_file /Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/UKB/Final_Data_UKB.csv \
--cohort_code UKB


Rscript bin/process.R \
Chapter_I \
NA \
UKB \
/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/UKB/Final_Data_UKB.csv
