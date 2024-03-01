import argparse
import os
from Bio import SeqIO

def split_file(input_file, output_dir, n):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Determine file format (FASTA or FASTQ)
    file_format = "fasta" if input_file.lower().endswith(".fasta") else "fastq"

    # Read input file and split records into n subsets
    records = list(SeqIO.parse(input_file, file_format))
    subset_size = len(records) // n

    for i in range(n):
        start_idx = i * subset_size
        end_idx = (i + 1) * subset_size if i < n - 1 else len(records)
        subset_records = records[start_idx:end_idx]
        
        # Write subset to output file
        subset_filename = os.path.splitext(os.path.basename(input_file))[0] + f"_sub{i+1}.{file_format}"
        subset_filepath = os.path.join(output_dir, subset_filename)
        with open(subset_filepath, "w") as output_handle:
            SeqIO.write(subset_records, output_handle, file_format)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split FASTA or FASTQ file into n subsets")
    parser.add_argument("-i", "--input", help="Input FASTA or FASTQ file", required=True)
    parser.add_argument("-n", "--number", help="Number of subsets to create", type=int, required=True)
    parser.add_argument("-o", "--output", help="Output directory", required=True)
    args = parser.parse_args()

    split_file(args.input, args.output, args.number)
