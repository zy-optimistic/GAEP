# GAAP
A genome assembly assessment pipeline.

## Introduction

GAAP is a pipeline to assess genome assembly. 

## Installation


## Dependencies

### Base accuracy (mapping)

* bwa v0.7.17-r1198-dirty
* samtools v1.9
* bcftools v1.9

### Base accuracy (k-mer)

* meryl v1.3
* merqury v1.3

### breakpoint detection

* Bio::DB::Sam v1.43 (perl module)
* minimap2 v2.17-r974-dirty
* samtools v1.9 

### busco

* Augustus (> 3.2.1)
* blastn
* hmmer v3.1b2

### SNV-depth dot plot

* minimap2 v1.9 (for TGS reads mapping). Can be install by conda with command "conda install perl-bio-samtools -c bioconda".
* bwa v0.7.17-r1198-dirty (for NGS reads mapping)
* Rscript v4.0.3 (for plotting)
* bedtools
* ggplot2 (R module)
* ggExtra (R module)

## Usage

GAAP

