#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $task = 'run_merqury.pl';
my $mode = 1;

my $usage = <<USAGE;
Program: $task
Version: v1.0
Author:  Zhang Yong

Command: 
    $task [option] assembly fastq1 fastq2

Options: 
    -d <PATH>       output directory
    -o <STR>        output prefix
    -t <INT>        threads use in meryl
    --meryl         path to meryl
    --merqury       path to merqury.sh
    -h              print help massage
USAGE

my ($threads,$dir,$out_prefix,$meryl,$merqury,$help);
GetOptions(
	"t:i"            => \$threads,
	"d:s"            => \$dir,
	"o:s"            => \$out_prefix,
	"meryl:s"        => \$meryl,
	"merqury:s"      => \$merqury,
	"h"              => \$help,
);

my $assembly = shift;
my $fq1      = shift;
my $fq2      = shift;

die $usage unless $assembly && $fq1 && $fq2;

#---------------------------------variates-------------------------------------#

$out_prefix   = "output" unless $out_prefix;
$out_prefix   = "$dir/$out_prefix" if $dir;
my $merqury_dir = dirname($merqury);
my $genome_size = (-s $assembly);
my $qv_file = $out_prefix.

#-----------------------------------main---------------------------------------#

##best kmer
my $best_k_path = $merqury_dir.'/best_k.sh';
my $best_kmer   = (split /\n/, (`$best_k_path $genome_size`))[2];
$best_kmer = int($best_kmer);
print STDERR "Use kmer length of $best_kmer.\n";

##run meryl
unless ($threads) {
	print STDERR "No threads has been set, by default, meryl uses up to 64 cpus and all memories available.\n";
}else {
	$ENV{OMP_NUM_THREADS}=$threads;
	print STDERR "Use $threads CPUs to run meryl.\n";
}

my $i = 1;
for ($fq1,$fq1) {
	my $meryl_cmd = "meryl k=$best_kmer count output $out_prefix.$i.meryl $_";
	_system($meryl_cmd, $mode);
	$i ++;
}
#merge the output of meryl
my $merge_cmd = "meryl union-sum output $out_prefix.merge.meryl $out_prefix.1.meryl $out_prefix.2.meryl";
_system($merge_cmd, $mode);


##run merqury
#/home/zhangyong/software/merqury-1.3/merqury.sh 02428.meryl ../02428.genome.fasta merqury_output
my $merqury_cmd = "$merqury $out_prefix.merge.meryl $assembly $out_prefix";
_system($merqury_cmd, $mode);


##-------------------------------subroutines----------------------------------##



sub _system {
	my $cmd = shift;
	my $mode = @_ ? shift : 0;
	if ( $mode == 1 ) {
		print $cmd,"\n";
	}else {
		print $cmd,"\n";
		die "[$task] Can't run \"$cmd\".\n" if ( system $cmd );
	}
	return;
}

sub check_software {
	my $software = shift;
	my $path = shift if @_;
	if ( $path ) {
		if ( basename($path) eq $software && -X $path ) {
			$software = $path;
		}else {
			return "-1";
		}
	}else {
		my @path = split /:/, $ENV{PATH};
		foreach ( @path ) {
			$_ =~ s/\/$//; 
			return $software = "$_/$software" if -X "$_/$software";
		}
		return "-1";
	}
}