import csv
from collections import defaultdict

def compile_insertion_statistics(input_csv, output_csv, filtered_csv):
    """
    Reads in a CSV file with barcodes and insertions, compiles statistics for identical barcodes,
    writes the results to a new CSV file, and also creates a CSV with barcodes meeting specific criteria.

    :param input_csv: Path to the input CSV file with barcodes and insertions.
    :param output_csv: Path to the output CSV file to write the general statistics.
    :param filtered_csv: Path to the output CSV file to write barcodes with specific criteria.
    """
    barcode_dict = defaultdict(list)

    # Read the input CSV and group insertions by barcode
    with open(input_csv, 'r') as infile:
        reader = csv.reader(infile)
        next(reader)  # Skip the header
        for row in reader:
            barcode, insertion = row
            barcode_dict[barcode].append(insertion)

    # Prepare the statistics for each unique barcode
    stats = []
    filtered_stats = []
    for barcode, insertions in barcode_dict.items():
        unfound_count = insertions.count("unfound")
        xxnone_count = insertions.count("xxnone")
        len_20_count = sum(1 for i in insertions if len(i) == 20)
        len_40_count = sum(1 for i in insertions if len(i) == 40)
        len_60_count = sum(1 for i in insertions if len(i) == 60)
        len_80_count = sum(1 for i in insertions if len(i) == 80)
        insertion_list = ";".join(insertions)
        stats.append([barcode, unfound_count, xxnone_count, len_20_count, len_40_count, len_60_count, len_80_count, insertion_list])

        # Filter based on specific criteria
        if xxnone_count >= 2 and all(v == 0 for v in [unfound_count, len_20_count, len_40_count, len_60_count, len_80_count]):
            filtered_stats.append([barcode])

    # Write the statistics to the output CSV
    with open(output_csv, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Barcode', 'Unfound Count', 'xxnone Count', 'Length 20 Count', 'Length 40 Count', 'Length 60 Count', 'Length 80 Count', 'Insertions List'])
        writer.writerows(stats)

    # Write the filtered barcodes to the filtered CSV
    with open(filtered_csv, 'w', newline='') as filteredfile:
        writer = csv.writer(filteredfile)
        writer.writerow(['Barcode'])
        writer.writerows(filtered_stats)

##USER INPUTS
input_csv = 'iMEF_biorep1_15day-static_insertion_combo.csv'  # Replace with your input CSV file path from the previous script
output_csv = 'iMEF_biorep1_15day-compiled_stats.csv'  # Replace with your desired output CSV file path
filtered_csv = 'iMEF_biorep1_15day-filtered_barcodes.csv'  # Replace with your desired filtered output CSV file path

compile_insertion_statistics(input_csv, output_csv, filtered_csv)

print(f"CSV file '{output_csv}' has been created with the compiled statistics.")
print(f"CSV file '{filtered_csv}' has been created with filtered barcodes.")
