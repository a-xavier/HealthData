import os

files = [ 
    "/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/Reference/icd_flat.csv",
    "/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/Reference/all_codes.txt",
    "/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/UKB/Final_Data_UKB.csv",
    "/Users/alex/Library/CloudStorage/OneDrive-TheUniversityofNewcastle/Cairns_Group/ICD10_Project_Official/Generation_Scotland/Final_Data_GS.csv"
]

# For each file 
# grep "Chapter " and replace spaces with underscores
# replace Original file 

for file in files:
    print(f"Processing file: {file}", end="", flush=True)
    with open(file, "r") as f:
        content = f.read()
    content = content.replace("Chapter ", "Chapter_")
    print(" - Done", flush=True)
    print(f"Writing file: {file}", end="", flush=True)
    with open(file, "w") as f:
        f.write(content)
    print(" - Done", flush=True)