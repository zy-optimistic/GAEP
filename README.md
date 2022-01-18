# GAAP
A genome assembly assessment pipeline.

## Introduction

GAAP is a pipeline to assess genome assembly. 

## Installation


## Dependencies

### Base accuracy (mapping)

* bwa
* samtools 
* bcftools

### Base accuracy (k-mer)

* meryl
* merqury

### breakpoint detection

* Bio::DB::Sam (perl module)
* minimap2
* samtools 

### busco

* Augustus
* blastn
* hmmer

### SNV-depth dot plot

* minimap2 (for TGS reads mapping)
* bwa (for NGS reads mapping)
* Rscript (for plotting)
* ggplot2 (R module)
* ggExtra (R module)

## Usage

GAAP

