import tqdm
import pandas as pd 
from collections import Counter
from itertools import combinations

# if all codes 14316 from all_codes.txt
# If observed codes 14200

# Files I have 
# Update when there are new cohorts added
files_i_have = [
# UKB
"/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/UKB/Final_Data_UKB.csv",
# GS
"/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/Generation_Scotland/Final_Data_GS.csv"
]

all_codes_observed = set()
for file in files_i_have:
    df = pd.read_csv(file)
    all_codes_in_this_cohort = df["all_levels"]
    # Codes are separated by pipes so split and add to the set
    for codes in all_codes_in_this_cohort.dropna():
        all_codes_observed.update(codes.split("|"))

print(f"Total unique codes observed: {len(all_codes_observed)}")

sorted_all_codes = sorted(all_codes_observed)

# Do combination but remove [a, a] pairs and order does not matter so [a, b] is the same as [b, a]
pairs = []
all_codes = sorted_all_codes
for i in tqdm.trange(len(all_codes)):
    for j in range(i + 1, len(all_codes)):
        pairs.append([all_codes[i], all_codes[j]])


print(f"Total pairs: {len(pairs)}")
# Write pairs to a file
with open("all_pairs.txt", "w") as f:
    for pair in pairs:
        f.write(f"{pair[0]}\t{pair[1]}\n")

# 102,466,770 total if using all codes 
# 100,812,900 total if using only observed codes

# Now only do observed pairs (meaning only codes that do co occur in real datasets)
observed_pair_counts = Counter()

for file_path in files_i_have:

    # Read only the required column
    all_levels = pd.read_csv(
        file_path,
        usecols=["all_levels"]
    )["all_levels"]

    for value in tqdm.tqdm(
        all_levels.dropna(),
        desc=f"Processing {file_path}"
    ):
        # Remove whitespace, empty values and duplicates within this row.
        # Sorting ensures that A|B and B|A both produce the pair (A, B).
        codes = sorted({
            code.strip()
            for code in value.split("|")
            if code.strip()
        })

        # A row with fewer than two distinct codes cannot produce a pair
        if len(codes) < 2:
            continue

        observed_pair_counts.update(combinations(codes, 2))


print(f"Total unique observed pairs: {len(observed_pair_counts):,}")
print(f"Total pair observations: {sum(observed_pair_counts.values()):,}")

with open("all_observed_pair.txt", "w") as output:
    for (code_1, code_2), count in observed_pair_counts.items():
        output.write(f"{code_1}\t{code_2}\t{count}\n")

# IF only observed
# Total unique observed pairs: 12,995,854
# This is like 10X less than just doing everything
# Total pair observations: 651,345,478