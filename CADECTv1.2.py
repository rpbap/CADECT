import argparse
import os
import time
import gzip
from Bio import SeqIO
from Bio import pairwise2


def parse_args():
    parser = argparse.ArgumentParser(description="Parse FASTA/FASTQ files and analyze sliding windows.")
    parser.add_argument("-R", "--reads", required=True, help="Input FASTA/FASTQ file (optionally gzipped)")
    parser.add_argument("-w", "--window_size", type=int, default=500, help="Window size (default: 500)")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory")
    return parser.parse_args()


def determine_format(file):
    """
    Determine the file format based on the file extension, supporting gzip extensions.
    """
    base, ext = os.path.splitext(file)
    if ext == ".gz":  # Check if it's gzipped
        _, ext = os.path.splitext(base)  # Get the extension before .gz

    if ext.lower() in [".fastq", ".fq"]:
        return "fastq"
    else:
        return "fasta"


def open_reads(input_file):
    """
    Open input FASTA/FASTQ file (possibly gzipped) and return (handle, file_format).
    Caller is responsible for closing the handle.
    """
    open_func = gzip.open if input_file.endswith(".gz") else open
    file_format = determine_format(input_file)
    handle = open_func(input_file, "rt")
    return handle, file_format


def count_reads(input_file, file_format):
    """
    Count total number of reads in the input file (streaming, low memory).
    """
    open_func = gzip.open if input_file.endswith(".gz") else open
    total = 0
    with open_func(input_file, "rt") as handle:
        for _ in SeqIO.parse(handle, file_format):
            total += 1
    return total


def split_sequence(record, window_size):
    """
    Split a sequence record into windows of a given size (non-overlapping).
    """
    windows = []
    seq = record.seq
    for i in range(0, len(seq), window_size):
        window = seq[i: i + window_size]
        windows.append(window)
    return windows


def get_overlap_count(windows, window_size):
    """
    Get the number of overlaps between windows based on pairwise similarity.
    NOTE: This is O(n^2) in the number of windows and can be slow for very long reads.
    """
    overlap_count = 0
    n = len(windows)
    for i in range(n):
        for j in range(i + 1, n):
            score = pairwise2.align.globalxx(windows[i], windows[j], score_only=True)
            align_length = max(len(windows[i]), len(windows[j]))
            if align_length == 0:
                continue
            if score / align_length >= 0.7:
                overlap_count += 1
    return overlap_count


def classify_sequences(handle,
                       file_format,
                       window_size,
                       log_file,
                       table_file,
                       short_handle,
                       concat_handle,
                       non_concat_handle,
                       total_reads):
    """
    Stream through the input reads, classify each, and write directly to
    the appropriate output file. No large lists are kept in memory.
    """
    num_parsed = 0

    for record in SeqIO.parse(handle, file_format):
        num_parsed += 1

        windows = split_sequence(record, window_size)
        num_windows = len(windows)
        num_overlaps = get_overlap_count(windows, window_size)

        if num_windows == 1:
            classification = "short"
            SeqIO.write(record, short_handle, file_format)
        elif num_overlaps > 0:
            classification = "putative_concatemers"
            SeqIO.write(record, concat_handle, file_format)
        else:
            classification = "non_concatemers"
            SeqIO.write(record, non_concat_handle, file_format)

        # Classification table
        table_file.write(
            f"{record.id}\t{classification}\t{num_windows}\t{num_overlaps}\n"
        )

        # Progress log with percentage
        if total_reads > 0:
            pct = (num_parsed / total_reads) * 100.0
            log_file.write(
                f"Classifying: {num_parsed} / {total_reads} ({pct:.2f}%)\n"
            )
        else:
            log_file.write(f"Classifying: {num_parsed}\n")

        log_file.flush()


def main():
    args = parse_args()
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # Determine file format once
    file_format = determine_format(args.reads)
    out_ext = "fq" if file_format == "fastq" else "fasta"

    # Open all output files immediately so they appear in the folder
    with open(os.path.join(output_dir, "progress.log"), "w") as log_file, \
         open(os.path.join(output_dir, "classification_table.txt"), "w") as table_file, \
         open(os.path.join(output_dir, f"short.{out_ext}"), "w") as short_handle, \
         open(os.path.join(output_dir, f"putative_concatemers.{out_ext}"), "w") as concat_handle, \
         open(os.path.join(output_dir, f"non_concatemers.{out_ext}"), "w") as non_concat_handle:

        start_time = time.time()
        log_file.write("Start time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
        log_file.write("Counting reads...\n")
        log_file.flush()

        # Pre-count total reads (streaming, low memory)
        total_reads = count_reads(args.reads, file_format)
        log_file.write(f"Total reads: {total_reads}\n")
        log_file.flush()

        # Header for the classification table
        table_file.write("Read ID\tClassification\tNum Windows\tNum Overlaps\n")

        # Now open input for streaming classification
        handle, _ = open_reads(args.reads)

        # Stream and classify
        classify_sequences(handle,
                           file_format,
                           args.window_size,
                           log_file,
                           table_file,
                           short_handle,
                           concat_handle,
                           non_concat_handle,
                           total_reads)

        handle.close()

        end_time = time.time()
        log_file.write("End time: {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S")))
        log_file.write("Total time: {:.2f} seconds\n".format(end_time - start_time))


if __name__ == "__main__":
    main()
