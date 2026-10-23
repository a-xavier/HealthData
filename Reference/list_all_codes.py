import pandas as pd 
df = pd.read_csv("icd_flat.csv")

columns_to_keep = ["chapter", "block", "category", "subcategory", "subclassification"]

all_codes = set()

for column in columns_to_keep:
    all_codes.update(df[column].dropna().to_list())


all_codes = sorted(all_codes)


# Write all codes
with open("all_codes.txt", "w") as f:
    for item in all_codes:
        f.write(f"{item}\n")