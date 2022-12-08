# CADET v.0.1
Concatemer by Amplification DEtection Tool V.0.1


Whole Genome Amplification using multiple displacement amplification (MDA) sometimes can introduce potential false concatemer sequences that can affect whole genome assembly assays. Here we propose a Concatemer detection tool for those WGA assays.



## How it works?
It splits all reads in separate files to perform sliding windows with the user prefered size and the gap between these windows. For ONT amplified reads, we suggest windows >= 500bp with no overlaps (e.g. <window size> = 500 and <slide size> = 500). If the read is not able to generate more than one window (< 1kb in size in the 500bp window example) the read is skipped and its ID stored in the short.txt output file. Reads with more than two windows, will have their fragment windows aligned using nucmer and reads with overlaps are reported in stats and their IDs stored in the concat_ID output file. Statistics with the total number of reads, number of putative concatemers, number with no concatemer detection and the overlap frequency can be found in the stats.txt output.
<p align="center">
<img width="350" height="350" alt="concatemer" src="https://user-images.githubusercontent.com/28576450/206510156-341b185c-b284-41a2-9937-484a46f24266.png">
</p>
Figure. Concatamer-Mediated Multiple Displacement Amplification. The principle of concatamermediated multiple displacement amplification. 1-Religation of DNA fragments with T4 DNA ligase, which leads to two types of products, 2.1-Linear Concatamers and 2.2-Circular Concatamers. 3.1 and 3.2-Annealing of random hexamer primers and addition of phi29-DNA polymerase leads to concatamers-mediated multiple displacement amplification from linear and circular concatamers respectively. (Shoaib et al., 2008)
  
### workflow

![CADET](https://user-images.githubusercontent.com/28576450/206311935-1010d792-fd71-467f-89ba-8f55dabafe9c.png)

## Installation

```
mamba create -n cadet -c bioconda -c conda-forge mummer seqtk 
```

## Usage
```
./CADET.sh <fasta> <window size> <slide size>
```

## Output Files
| Output File | Description |
| --- | --- |
|`concat_IDs`|sequence ID of putative concatemeric reads (useful for `seqtk subseq file.fastq **concat_IDs** > conc.fastq` )|
|`concat_count.tab`||
|`coords (folder)`||
|`short.txt`||
|`stats.txt`||
