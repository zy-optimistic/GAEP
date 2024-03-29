# GAEP
A genome assembly evaluating pipeline.

## Introduction

GAEP is a pipeline to assess genome assembly. 

## Installation
```shell
git clone https://github.com/zy-optimistic/GAEP.git
cd GAEP  
./gaep  
```
## Usage
```
gaep <command> [options]

pipe   (NGS,TGS,trans)  let GAEP to determine the module to be executed based on the input data   
stat                    report genome basic information  
macc   (NGS)            base accuracy based on reads mapping  
kacc   (NGS)            base accuracy based on K-mer  
bkp    (TGS)            misassembly breakpoints detected  
snvcov (NGS,TGS)        SNV-coverage dot plot  
busco                   run busco v5  
```

## Running example
```bash
gaep pipe -r genome.fasta --lr TGS.fastq -x pb --sr1 NGS_1.fastq --sr2 NGS_2.fastq -t 3 -c config.txt
#You can list your data and dependancies in the config.txt. The template of config.txt is in GAEP/config/.
```

## Dependencies

The specified versions were used for testing GAEP.

### Misassembliy breakpoint detection (bkp)

* Bio::DB::HTS (perl module). Can be installed by conda with command "conda install -c bioconda perl-bio-db-hts".
* minimap2 v2.17-r974-dirty
* samtools v1.9 

### SNV-Cov dot plot (snvcov)

* minimap2 v2.17-r974-dirty (for TGS reads mapping)
* bwa v0.7.17-r1198-dirty (for NGS reads mapping)
* bcftools v1.9 (for SNV calling)
* Rscript v4.0.3 (for plotting)
* bedtools
* ggplot2 (R module)
* ggExtra (R module)

### Base accuracy (macc)

* bwa v0.7.17-r1198-dirty
* samtools v1.9
* bcftools v1.9

### Base accuracy (kacc)

* meryl v1.3
* merqury v1.3

### busco v5

* busco v5
* metaeuk (for eukaryotes)
* prodigal (for non-eukaryotes)
* hmmsearch
* bbtools
