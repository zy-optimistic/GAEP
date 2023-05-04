# GAAP
A genome assembly assessment pipeline.

## Introduction

GAAP is a pipeline to assess genome assembly. 

## Installation

git clone https://github.com/zy-optimistic/GAAP.git
cd GAAP
./gaap

## Usage

gaap <command> [options]

pipe       let GAAP to determine the module to be executed based on the input data
stat       report genome basic information
macc       base accuracy based on reads mapping
kacc       base accuracy based on K-mer
bkp        misassembly breakpoints detected
snvcov     SNV-coverage dot plot
busco      run busco v3

## Dependencies

### Misassembliy breakpoint detection (bkp)

* Bio::DB::HTS (perl module). Can be installed by conda with command "conda install -c bioconda perl-bio-db-hts".
* minimap2 v2.17-r974-dirty
* samtools v1.9 

### SNV-depth dot plot (snvcov)

* minimap2 v1.9 (for TGS reads mapping)
* bwa v0.7.17-r1198-dirty (for NGS reads mapping)
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

### busco

* Augustus (> 3.2.1)
* blastn
* hmmer v3.1b2


