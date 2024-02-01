import argparse
import random
from Bio import SeqIO

def generate_chimeric_read(sequence):
    length = len(sequence)
    split_point = random.randint(1, length - 1)
    chimeric_read = sequence[:split_point] + sequence[split_point:][::-1]
    return chimeric_read

def generate_random_chimeric_reads(sequences, num_reads):
    chimeric_reads = []
    for _ in range(num_reads):
        random_sequence = random.choice(sequences)
        chimeric_reads.append(generate_chimeric_read(random_sequence))
    return chimeric_reads

def read_fastq(file_path):
    sequences = []
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            sequences.append(str(record.seq))
    return sequences

def main():
    parser = argparse.ArgumentParser(description="Generate random chimeric reads from a FastQ file.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FastQ file")
    parser.add_argument("-n", "--num_reads", type=int, required=True, help="Number of chimeric reads to generate")
    parser.add_argument("-o", "--output", required=True, help="Path to the output file for generated chimeric reads")
    args = parser.parse_args()

    input_sequences = read_fastq(args.input)
    generated_chimeric_reads = generate_random_chimeric_reads(input_sequences, args.num_reads)

    with open(args.output, "w") as output_file:
        for idx, chimeric_read in enumerate(generated_chimeric_reads, start=1):
            output_file.write(f">Chimeric_Read_{idx}\n{chimeric_read}\n")

if __name__ == "__main__":
    main()
