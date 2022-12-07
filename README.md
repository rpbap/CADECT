# CADET v.0.1
Concatemer by Amplification DEtection Tool V.0.1


Whole Genome Amplification using multiple displacement amplification (MDA) sometimes can introduce potential false concatemer sequences that can affect whole genome assembly assays. Here we propose a Concatemer detection tool for those WGA assays.

## How it works?
![CADET](https://user-images.githubusercontent.com/28576450/206311935-1010d792-fd71-467f-89ba-8f55dabafe9c.png)

## Installation

```
mamba create -n cadet -c bioconda -c conda-forge mummer seqtk 
```

## Usage
```
./CADET.sh <fasta> <window size> <slide size>
```
