#! /usr/bin/perl -w

use strict;
use Getopt::Long;

my $task = "run_merqury";
my $mode = 0;



##! /bin/bash
#
#OMP_NUM_THREADS=$1
#k=21
#genome="NA19240"
#
#if [ ! $OMP_NUM_THREADS ]
#        then
#        echo [Error]Input OMP_NUM_THREADS number!
#        exit -1
#fi
#
#echo $OMP_NUM_THREADS
#
#date
#echo Now running Meryl
#
#for i in GM19240_NoIndex_L007_R1 GM19240_NoIndex_L007_R2 GM19240_NoIndex_L008_R1 GM19240_NoIndex_L008_R2
#do
#        echo meryl k=$k count output meryl/$i.meryl ../NGS/${i}_001.fastq.gz
#        meryl k=$k count output meryl/$i.meryl ../NGS/${i}_001.fastq.gz
#done
#
### 合并两个reads 的db
#echo meryl union-sum output $genome.meryl meryl/*.meryl
#meryl union-sum output $genome.meryl meryl/*.meryl
#
#echo /home/zhangyong/software/merqury-1.3/merqury.sh $genome.meryl NA19240_pri.mat.fasta merqury_output
#/home/zhangyong/software/merqury-1.3/merqury.sh $genome.meryl NA19240_pri.mat.fasta merqury_output
