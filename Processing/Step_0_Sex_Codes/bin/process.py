#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description="Process sex codes")
    parser.add_argument("--code", required=True, help="the code to examine")
    parser.add_argument("--final_df_file", required=True, help="Path to final_df.csv")
    parser.add_argument("--threshold", type=int, default=90, help="Threshold for sex imbalance")
    parser.add_argument("--cohort_code", default="default_cohort", help="Cohort code")
    args = parser.parse_args()

    code = args.code
    final_df_file = args.final_df_file
    threshold = args.threshold
    cohort_code = args.cohort_code

    # Your processing code here
    print(f"Processing {code} and {final_df_file} with threshold {threshold} for cohort {cohort_code}")
    # Create dummy output
    # sex_codes_output_${cohort_code}.txt
    output_file = f"sex_codes_output_{code}_{cohort_code}.txt"

    # Read df 
    df = pd.read_csv(final_df_file)

    # Create 0 1 phenotype by grepping code from all_levels column
    df['phenotype'] = df['all_levels'].apply(lambda x: 1 if code in str(x) else 0)

    print(df["phenotype"])

    # phenotype count (otherwise not very informative)
    phenotype_count = df['phenotype'].sum()

    print(f"Phenotype count for {code}: {phenotype_count}")

    if phenotype_count >= 5:
        tmp_df = df[df['phenotype'] == 1]
        prop_male = tmp_df['is_male'].mean()
        if prop_male >= threshold / 100:
            sex_flag = "M"
        elif prop_male <= (1 - threshold / 100):
            sex_flag = "F"
        else:
            sex_flag = "NA"

    else:
        print(f"Phenotype count for {code} is less than 5")
        sex_flag = "NA"
        prop_male = "NA"

    # output format is one line separated by tab 
    # code sex flag and male proportion
    with open(output_file, "w") as f:
        f.write(f"{code}\t{sex_flag}\t{prop_male}\t{phenotype_count}\n")
    

if __name__ == "__main__":
    main()