###This script takes as its input one of the output files from
###extract_insertion_sequences_and_report_indel_sizes_and_signatures.py
###It is used on nanopore data to sum all of the reads that are at or near
###the expected size for a given number of edits. This is necessary because
###nanopore sequencing often makes indel errors. (This work was all done using
###Oxford Nanopore Ligation Kit 14.)
###Written by ChatGPT3.5 and Theresa Loveless.

###This script was used for Figure 3 and Extended Data Figures 4-5.

import pandas as pd

# Create an empty DataFrame to store the results
result_df = pd.DataFrame()

# Specify the path pattern for your CSV files (e.g., 'file_{}.csv') and file names
file_path_pattern = 'PER104-{}_indel_histogram.csv'
sample_path_pattern = 'PER104-{}'

# Specify the column name you want to extract
column_name = 'counts'

# Create empty lists to store the sums and file names
samples = []
sums_wt = []
sums_one_insertion = []
sums_two_insertion = []
sums_three_insertion = []
sums_four_insertion = []
sums_five_insertion = []
sums_six_insertion = []
sums_seven_insertion = []
sums_eight_insertion = []

#Loop through files to make a list of file names
for file_number in range(1,45):
    sample = sample_path_pattern.format(file_number)
    samples.append(sample)

# Loop through files 1 to 64
for file_number in range(1, 45):
    # Construct the file path for the current file
    file_path = file_path_pattern.format(file_number)
    
    # Read the CSV file into a DataFrame
    df = pd.read_csv(file_path)
    
    # Extract rows 17 to 25 from the desired column and calculate the sum
    sum_result_wt = df.loc[17:26, column_name].sum()
    sums_wt.append(sum_result_wt)
    
    # Extract rows 37 to 45 from the desired column and calculate the sum
    sum_result_one = df.loc[37:46, column_name].sum()
    sums_one_insertion.append(sum_result_one)

    # Extract rows 57 to 65 from the desired column and calculate the sum
    sum_result_two = df.loc[57:66, column_name].sum()
    sums_two_insertion.append(sum_result_two)
    
    # Extract rows 77 to 85 from the desired column and calculate the sum
    sum_result_three = df.loc[77:86, column_name].sum()
    sums_three_insertion.append(sum_result_three)

    # Extract rows 97 to 105 from the desired column and calculate the sum
    sum_result_four = df.loc[97:106, column_name].sum()
    sums_four_insertion.append(sum_result_four)
    
    # Extract rows 117 to 125 from the desired column and calculate the sum
    sum_result_five = df.loc[117:126, column_name].sum()
    sums_five_insertion.append(sum_result_five)

    # Extract rows 137 to 145 from the desired column and calculate the sum
    sum_result_six = df.loc[137:146, column_name].sum()
    sums_six_insertion.append(sum_result_six)

    # Extract rows 157 to 165 from the desired column and calculate the sum
    sum_result_seven = df.loc[157:166, column_name].sum()
    sums_seven_insertion.append(sum_result_seven)
    
    # Extract rows 177 to 185 from the desired column and calculate the sum
    sum_result_eight = df.loc[177:186, column_name].sum()
    sums_eight_insertion.append(sum_result_eight)

# Create the result DataFrame with one column for each file
result_df['sample'] = samples
result_df['wt'] = sums_wt
result_df['one_insertion'] = sums_one_insertion
result_df['two_insertion'] = sums_two_insertion
result_df['three_insertion'] = sums_three_insertion
result_df['four_insertion'] = sums_four_insertion
result_df['five_insertion'] = sums_five_insertion
result_df['six_insertion'] = sums_six_insertion
result_df['seven_insertion'] = sums_seven_insertion
result_df['eight_insertion'] = sums_eight_insertion

result_df.to_csv('total_summary.csv', index=False)

# Display the resulting DataFrame
print(result_df)
