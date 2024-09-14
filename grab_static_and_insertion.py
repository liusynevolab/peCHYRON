import re
import csv

def extract_sequences(fastq_file, start_sequence_insertion, end_sequence_insertion, start_sequence_barcode, end_sequence_barcode):
    """
    Extracts both the insertion and barcode from each sequence in the FASTQ file.
    If the upstream and downstream sequences are not found for the insertion or barcode, "unfound" is returned.
    If the sequences are found but there is nothing in between, "xxnone" is returned.

    :param fastq_file: Path to the input FASTQ file
    :param start_sequence_insertion: The known start sequence for the insertion (string)
    :param end_sequence_insertion: The known end sequence for the insertion (string)
    :param start_sequence_barcode: The known start sequence for the barcode (string)
    :param end_sequence_barcode: The known end sequence for the barcode (string)
    :return: List of tuples containing (barcode, insertion) for each entry
    """
    results = []

    with open(fastq_file, 'r') as file:
        line_number = 0
        for line in file:
            line_number += 1
            if line_number % 4 == 2:  # Sequence line in FASTQ format
                sequence = line.strip()

                # Extract the insertion
                insertion_match = re.search(f'{start_sequence_insertion}(.*?){end_sequence_insertion}', sequence)
                if insertion_match:
                    insertion = insertion_match.group(1)
                    insertion = "xxnone" if insertion == "" else insertion
                else:
                    insertion = "unfound"

                # Extract the barcode
                barcode_match = re.search(f'{start_sequence_barcode}(.*?){end_sequence_barcode}', sequence)
                if barcode_match:
                    barcode = barcode_match.group(1)
                    barcode = "xxnone" if barcode == "" else barcode
                else:
                    barcode = "unfound"

                results.append((barcode, insertion))

    return results

def write_to_csv(output_file, data):
    """
    Writes the barcode and insertion data to a CSV file.

    :param output_file: Path to the output CSV file
    :param data: List of tuples containing (barcode, insertion)
    """
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Barcode', 'Insertion'])  # Write the header
        writer.writerows(data)  # Write the data rows


##USER INPUTS
fastq_file = 'iMEF_biorep1_15day.fastq'  # Replace with your FASTQ file path
output_file = 'iMEF_biorep1_15day-static_insertion_combo.csv'  # Replace with your desired output CSV file path

start_sequence_insertion = 'CAAGCCATCATCACCATCAT'  # Replace with your known start sequence for the insertion
end_sequence_insertion = 'TGATGGCAGAGGAAAGGAAG'  # Replace with your known end sequence for the insertion

start_sequence_barcode = 'AGGAAGCCCTGCTTCCTCCA'  # Replace with your known start sequence for the barcode
end_sequence_barcode = 'CTGAAATTCTCGACCTCGAG'  # Replace with your known end sequence for the barcode

data = extract_sequences(fastq_file, start_sequence_insertion, end_sequence_insertion, start_sequence_barcode, end_sequence_barcode)
write_to_csv(output_file, data)

print(f"CSV file '{output_file}' has been created with barcode and insertion data.")
