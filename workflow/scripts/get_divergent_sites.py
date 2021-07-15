# Import libraries
import os.path
import pandas as pd

# Set input/output files
input_files = snakemake.input
output_file = snakemake.output[0]

# Get sample names
sample_names = [os.path.splitext(os.path.basename(sample_file))[0] for sample_file in input_files]

# Create empty list to hold DFs
df_list = []

# Read CSVs into list
for i in range(len(input_files)):
    # set column names
    col_names = ["CHROM", "POS_1", "POS_2", "REF", "ALT", sample_names[i]]
    # read in DF
    temp_df = pd.read_csv(input_files[i], names = col_names)
    # append to list
    df_list.append(temp_df)

# Merge DFs
merged = pd.merge(df_list[0], df_list[1], on=["CHROM", "POS_1", "POS_2", "REF", "ALT"], how="inner")

# take last two columns and remove `/` and `|`
df_new = pd.DataFrame()
df_new[sample_names] = merged.loc[:, sample_names].replace('/|\|', '', regex = True)

# find which rows are divergent
out = merged.loc[df_new[sample_names[0]] != df_new[sample_names[1]]].copy()
# replace `|` with `/` in F0-1 (e.g. Cab) and F0-2 (e.g. Kaga) columns
out.loc[:, sample_names] = out.loc[:, sample_names].replace('\|', '/', regex = True)

# Write to file
out.to_csv(output_file, sep = '\t', header = False, index = False)