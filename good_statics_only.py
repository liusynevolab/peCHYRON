import csv
import sys

def filter_output_csv(output_csv, filtered_csv, updated_output_csv):
    """
    Updates the output CSV to only include barcodes that are in the filtered CSV list.
    # Increase the field size limit to handle larger fields
    :param output_csv: Path to the original CSV file with all statistics.
    :param filtered_csv: Path to the CSV file with barcodes to keep.
    :param updated_output_csv: Path to the updated CSV file that only includes filtered barcodes.
    """
    # Read the filtered barcodes into a set for quick lookup

    csv.field_size_limit(sys.maxsize)
    filtered_barcodes = set()
    with open(filtered_csv, 'r') as filteredfile:
        reader = csv.reader(filteredfile)
        next(reader)  # Skip the header
        for row in reader:
            barcode = row[0]
            filtered_barcodes.add(barcode)

    # Read the original output CSV and write only rows with barcodes in the filtered list
    with open(output_csv, 'r') as infile, open(updated_output_csv, 'w', newline='') as outfile:
        reader = csv.reader(infile)
        writer = csv.writer(outfile)
        
        header = next(reader)  # Read and write the header
        writer.writerow(header)
        
        for row in reader:
            barcode = row[0]  # Assuming barcode is in the first column
            if barcode in filtered_barcodes:
                writer.writerow(row)

##USER INPUTS
output_csv = 'iMEF_biorep1_15day-compiled_stats.csv'  # Replace with your original CSV file path
filtered_csv = 'iMEF_biorep1_zero_timepoint-filtered_barcodes.csv'  # Replace with your filtered CSV file path
updated_output_csv = 'iMEF_biorep1_15day-updated_compiled_stats.csv'  # Replace with your desired updated CSV file path

filter_output_csv(output_csv, filtered_csv, updated_output_csv)

print(f"CSV file '{updated_output_csv}' has been created with only the filtered barcodes.")
