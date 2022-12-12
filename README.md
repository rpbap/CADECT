# CADECT v.0.2.1
Concatemer by Amplification DEteCtion Tool V.0.2.1


Whole Genome Amplification using multiple displacement amplification (MDA) sometimes can introduce potential false concatemer sequences that can affect whole genome assembly assays. Here we propose a Concatemer detection tool for those WGA assays.

<p align="center">
  
<img width="400" height="300" alt="CADET_AI" src="https://user-images.githubusercontent.com/28576450/206807841-2de5a0b3-4e00-460a-aaf1-34576894bf85.png">


</p>

**Figure. Concatemer-Mediated Multiple Displacement Amplification.** Annealing of random hexamer primers and addition of phi29-DNA polymerase leads to concatemers-mediated multiple displacement amplification from *(A)* linear and *(B)* circular concatemers respectively.  [Modified from Shoaib et al., 2008](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-9-415).

## How it works?
It splits all reads in separate files to perform sliding windows with the user prefered size and the gap between these windows. For ONT amplified reads, we suggest windows >= 500bp with no overlaps (e.g. `-w 500` and `-s 500`). If the read is not able to generate more than one window (< 1kb in size in the 500bp window example) the read is skipped and its ID stored in the `short.txt` output file. Reads with more than two windows, will have their fragment windows aligned using nucmer and reads with overlaps are reported in stats and their IDs stored in the `concat_ID` output file. Statistics with the total number of reads, number of putative concatemers, number with no concatemer detection and the overlap frequency can be found in the `stats.txt` output.
Fastq/Fasta files with the characterized reads will also be generated.

  
### Workflow

  <p align="center">

<img width="500" height="500" alt="CADET_0 2 1" src="https://user-images.githubusercontent.com/28576450/206768044-3ed65ee6-b119-470c-bed7-a361d497efb4.png">

    
  </p>

## Instalation

**Software Requirements:**

- Mummer v3.23

- Seqtk v1.3

**Easy install unisng conda/mamba**

```
mamba create -n cadect -c bioconda -c conda-forge mummer seqtk 
conda activate cadect
git clone https://github.com/rpbap/CADECT.git
cd CADECT
chmod a+x CADECT_v.0.2.1.sh
```

## Usage
```
./CADECT_v.0.2.1.sh [OPTIONS] -R <Reads.fastq> -w <window size> -s <slide size> -p <your_prefix>

Flag description:

Required:
  -R  --reads     fastq (or fasta*) file with reads generated by WGA sequencing using ONT (required)

Options: 
    -w  --window    length of desired window sequences in bp (required) (default = 500)
    -s  --slide     length to slide each window over in bp (required) (default = 500)
    -p  --prefix    Prefix name for your output folder (default = "CADECT_output")  
    -h  --help      display this message

*Note: if using a FASTA file as input, the Concat.fastq output will be in fasta format
```

## Output Files
| Output File | Description |
| --- | --- |
|`concat_IDs`|sequence ID of putative concatemeric reads|
|`coords (folder)`| Folder with all nucmer coord alignments|
|`short.txt`| IDs of reads detected as short and skipped by the pipeline|
|`stats.txt`| File statistics of the CADECT pipeline|
|`Non_conc.fastq`|fastq/fasta file containing non-concatemeric reads|
|`Conc.fastq`|fastq/fasta file containing putative concatemeric reads|
|`Short.fastq`|fastq/fasta file containing short reads|

### stats output example

Short sequences:                          2   

Number of non-concatemer detected:        2   

Number of putative concatemer detected:   5   

Total number of Reads:                    9   

`###` putative concatemers `###`

| read_number *(useful for coords)* | read_ID | self_alignemts |
| --- | --- | ---|
|window_4|d159b5a3-ee3b-4cc4-92ad-1422bf7a5a28|17.5|  
|window_5|159ffb63-2583-4a7d-88a5-639111d4fe99|104|
|window_6|6d5ce662-395e-4af2-a68c-37015af5913b|54|
|window_8|b8194fa6-aa7b-4017-bd55-5538b8f31039|142.5|
|window_9|a6b76c03-832a-47a1-bb80-0a57b862118a|9|


## Acknowledgements

The sliding window logic was made as a modification from [Kim Dill-McFarland, UW 2020](https://github.com/kdillmcfarland/sliding_windows).
