# Support scripts usages

## Chimera_generator.py

* It gets your input fastq file and generate a fasta file with simulated inverted chimeras generated to mimicry template switching events. This script just generates chimeras witha aminimum of 1kb length. The purpose of this script was to test the capability of CADECT to detect chimeras in our simulated datasets.
* Usage:
   `python script.py -i path/to/your/input.fastq -n 100 -o path/to/your/output.fasta`
* Parameters:

  **-n** number of chimeric reads to be generated
