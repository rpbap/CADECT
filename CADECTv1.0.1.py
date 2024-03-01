import argparse
import os
import time
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

def parse_args():
    parser = argparse.ArgumentParser(description="Parse FASTA/FASTQ files and analyze sliding windows.")
    parser.add_argument("-R", "--reads", required=True, help="Input FASTA/FASTQ file")
    parser.add_argument("-w", "--window_size", type=int, default=500, help="Window size (default: 500)")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    return parser.parse_args()

def parse_reads(input_file):
    """
    Parse input FASTA/FASTQ file and return a list of SeqRecord objects.
    """
    records = []
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, determine_format(input_file)):
            records.append(record)
    return records

def determine_format(file):
    """
    Determine the file format based on the file extension.
    """
    _, file_extension = os.path.splitext(file)
    if file_extension.lower() == ".fastq":
        return "fastq"
    else:
        return "fasta"

def classify_sequences(records, window_size, log_file, table_file):
    """
    Classify sequences as short, putative concatemers, or non-concatemers based on the number of windows.
    """
    short_sequences = []
    concatemers = []
    non_concatemers = []
    
    num_parsed = 0
    total_reads = len(records)
    
    for record in records:
        num_parsed += 1
        windows = split_sequence(record, window_size)
        num_windows = len(windows)
        num_overlaps = get_overlap_count(windows, window_size)
        if num_windows == 1:
            classification = "short"
            short_sequences.append(record)
        elif num_overlaps > 0:
            classification = "putative_concatemers"
            concatemers.append(record)
        else:
            classification = "non_concatemers"
            non_concatemers.append(record)
        
        # Write classification to output file
        table_file.write(f"{record.id}\t{classification}\t{num_windows}\t{num_overlaps}\n")
        
        # Write progress to log file
        log_file.write(f"Classifying: {num_parsed}/{total_reads}\n")
        log_file.flush()  # Ensure the message is written immediately
        
    return short_sequences, concatemers, non_concatemers

def split_sequence(record, window_size):
    """
    Split a sequence record into sliding windows of a given size.
    """
    windows = []
    for i in range(0, len(record.seq), window_size):
        window = record.seq[i:i+window_size]
        windows.append(window)
    return windows

def get_overlap_count(windows, window_size):
    """
    Get the number of overlaps between sliding windows.
    """
    overlap_count = 0
    for i in range(len(windows)):
        for j in range(i+1, len(windows)):
            alignments = pairwise2.align.globalxx(windows[i], windows[j], score_only=True)
            align_length = max(len(windows[i]), len(windows[j]))
            if alignments / align_length >= 0.7:
                overlap_count += 1
    return overlap_count

def main():
    args = parse_args()
    output_dir = args.output_dir
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Parse input file
    records = parse_reads(args.reads)
    
    # Open log file for writing progress
    with open(os.path.join(output_dir, "progress.log"), "w") as log_file, open(os.path.join(output_dir, "classification_table.txt"), "w") as table_file:
        start_time = time.time()
        log_file.write("Start time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
        
        # Write table header
        table_file.write("Read ID\tClassification\tNum Windows\tNum Overlaps\n")
        
        # Classify sequences
        short_seqs, concatemers, non_concatemers = classify_sequences(records, args.window_size, log_file, table_file)
        
        # Write sequences to output files
        for seq_list, classification in [(short_seqs, "short"), (concatemers, "putative_concatemers"), (non_concatemers, "non_concatemers")]:
            output_file = os.path.join(output_dir, f"{classification}.{determine_format(args.reads)}")
            SeqIO.write(seq_list, output_file, determine_format(args.reads))
        
        end_time = time.time()
        log_file.write("End time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
        log_file.write("Total time: {:.2f} seconds\n".format(end_time - start_time))

if __name__ == "__main__":
    main()

