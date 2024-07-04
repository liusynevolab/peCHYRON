# This script was used to generate data for Figure 5c and Extended Data Figure 8c. It is part the "record stretching" analysis pipeline.

# This script takes as an input the [sample name]_signatures_dataframes.csv file that is generated for each sample by extract_insertion_sequences_and_report_indel_sizes_and_signatures.py.
# Records in the [sample name]_signatures_dataframes.csv file are filtered; if any record contains an unexpected signature (presumably from sequencing error), the entire record is thrown out.

# The output is simply a dataframe of the records after the filter has been applied: [sample name]_signatures_dataframe_removed_unexpected_sigs.csv

import sys
import csv
import pandas as pd

########## USER INPUTS ##########

expected_odd_sigs = ['C','L','T']
expected_even_sigs = ['c','l','t']
base_input_file_name = 'PER71-'   #what's the constant part at the beginning of the signatures_dataframe file names? For our data, this is the string immediately preceding the sample number. Here you are providing the names of the input files.
input_identifiers = [3,4]     #what sample numbers are you analyzing?  You can manually type in a list of sample numbers and comment out the For loop below, or leave this list empty and populate it with the For loop below.
'''for num in range(1,55):
    input_identifiers.append(num)'''

########## END OF USER INPUTS ##########

for file in input_identifiers:

    full_input_file = '{}{}_signatures_dataframe.csv'.format(base_input_file_name, file)
    raw_data = pd.read_csv(full_input_file)
    data = raw_data.dropna(axis=0) #removing rows that have NaN (it's just one row, because a wt record is denoted as an empty string at the top of the signatures_dataframe).  This avoids errors, plus we don't want to stretch the wt records.

    list_of_records = []
    for row in data.loc[:, 'signatures (1 letter)']: #making a list of the record types
        list_of_records.append(str(row))

    list_of_counts = []
    for counts in data.loc[:, 'counts']: #making a list of the counts that correspond to each record type
        list_of_counts.append(counts)

    df_filtered_records = []
    for record, freq in zip(list_of_records, list_of_counts):
        odd_or_even = 1
        all_expected = 'true'
        for char in record:
            if odd_or_even%2 == 1:
                if char not in expected_odd_sigs:
                    all_expected = 'false'
            if odd_or_even%2 == 0:
                if char not in expected_even_sigs:
                    all_expected = 'false'
            odd_or_even += 1
        
        if all_expected == 'true':
            df_filtered_records.append([record, len(record), freq])

    headers = ['record', 'length', 'counts']
    final_dataframe = pd.DataFrame(df_filtered_records, columns = headers)

    final_dataframe.to_csv('{}{}_signatures_dataframe_removed_unexpected_sigs.csv'.format(base_input_file_name, file))
