#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
Getopt::Long::Configure qw( bundling no_ignore_case );

my $task = "stat_accuracy_base";
my $mode = 0;

my $usage = qq(
Program: $task
Version: v1.0

Usage:   $task 
Options:
         -r <FILE>           genome assembly in a fasta file
         -i <FILE>           fastq1 file. Can be set more than one time
         -I <FILE>           fastq2 file. Can be set more than one time
         -l <FILE>           file listing fastq files. Read1 and read2 delimited with tab
         -b <FILE>           bam file, will ignore -i, -I and -l
         -v <FILE>
         -d <INT>            output directory
         -o <INT>            output prefix
         --bcftools <PATH>   path to bcftools
         --bedtools <PATH>   path to bedtools
         -h                  print help massage
);

##----------------------------------options-----------------------------------##

my ($assembly,@reads1,@reads2,$file_list,$threads,$vcf,$bam,$dir,$prefix,$bcftools,$bedtools,$help);
GetOptions (
	"r:s"        => \$assembly,
	"i:s"        => \@reads1,
	"I:s"        => \@reads2,
	"l|list:s"   => \$file_list,
	"t:i"        => \$threads,
	"v:s"        => \$vcf,
	"b:s"        => \$bam,
	"d:s"        => \$dir,
	"o:s"        => \$prefix,
	"h"          => \$help,
	"bcftools:s" => \$bcftools, 
	"bedtools:s" => \$bedtools,
);
die $usage if $help;
##Check input files.
##Check input files and directory.
die $usage if !@reads1 && !@reads2 && !$file_list && !$bam && !$vcf;

$dir = "gaap_${task}_$$" unless $dir;
$dir =~ s/\/$//;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix = "accu_output_$$" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;


##Check software.
$bcftools = $bcftools ? check_software("bcftools", $bcftools) : check_software("bcftools");
die "[$task] Error! Bcftools is not found.
[$task] Check your config file or set it in your environment.\n" if $bcftools eq "-1";
$bedtools = $bedtools ? check_software("bedtools", $bedtools) : check_software("bedtools");
die "[$task] Error! Bedtools is not found.
[$task] Check your config file or set it in your environment.\n" if $bedtools eq "-1";

if (!$vcf && !$bam) {
	my $map_cmd = "$RealBin/run_bwa.pl -r $assembly ";
	for (my $i = 0; $i <= $#reads1 && $i <= $#reads2; $i++) {
		$map_cmd .= "-i $reads1[$i] -I $reads2[$i] ";
	}
	$map_cmd .= "-l $file_list " if $file_list;
	$map_cmd .= "-t $threads "   if $threads;
	$map_cmd .= "-d $dir/mapping/NGS_mapping ";
	$map_cmd .= "-o $prefix ";
	_system($map_cmd, $mode);
	$bam = "$dir/mapping/NGS_mapping/$prefix.paired.sorted.bam";
	
	my $var_cmd = "$RealBin/var_calling.pl -r $assembly ";
	$var_cmd .= "-b $bam ";
	$var_cmd .= "-t $threads "   if $threads;
	$var_cmd .= "-d $dir/variants/ ";
	$var_cmd .= "-o $prefix ";
	_system($var_cmd, $mode);
	$vcf = "$dir/variants/$prefix.flt.vcf"
}

##---------------------------------variants-----------------------------------##

my $mean_depth = 0;
my $total_len  = 0;
my $map_len    = 0;
my $var_num    = 0;
my $err_base   = 0;
my $accuracy   = 0;
my $QV         = 0;

my $out_file = "$prefix_out.base.qv";
##-----------------------------------main-------------------------------------##

##mappable length
my $genomecov = "$prefix_out.genomecov";
my $bedtools_cmd = "$bedtools genomecov -ibam $bam -split > $genomecov";
_system($bedtools_cmd, $mode);
open COV, '<', $genomecov or die "[$task] Can't open such file: $genomecov.\n";
while (<COV>) {
	chomp;
	my @line = split;
	if ($line[0] eq 'genome') {
		$mean_depth += $line[1] * $line[2];
		$total_len  += $line[2];
	} 
}
$mean_depth = $mean_depth / $total_len;
seek COV, 0, 0;
while (<COV>) {
	chomp;
	my @line = split;
	$map_len += $line[2] if ($line[0] eq 'genome' && $line[1] >= 3 && $line[1] <= $mean_depth*12 );
}
close COV;

##var number
open VCF, '-|', "$bcftools view -Hi '%QUAL>=10 && DP>=3 && DP<=$mean_depth*12 && (GT=\"AA\" || GT=\"Aa\")' $vcf" || 
    die "Fail to open pipe: $bcftools view -Hi '%QUAL>=10 && DP>=3 && DP<=$mean_depth*12 && (GT=\"AA\" || GT=\"Aa\")' $vcf.\n";
while (<VCF>) {
	#print $_;
	my ($ref,$alt) = (split)[3,4];
	my $alt_len = 0;
	if ( $alt =~ /,/ ) {
		my @alt_s = split /,/, $alt;
		$alt_len = $alt_len > length $_ ? $alt_len : length $_ for @alt_s;
	}else {
		$alt_len = length $alt;
	}
	$var_num += (length $ref) == $alt_len ? $alt_len : abs((length $ref) - $alt_len);
}
close VCF;

$err_base = $var_num/$map_len;
if ($err_base==0) {
	$QV = 0;
	print "[$task] No error base detected. QV can't be caculated.\n";
}else {
	$QV = sprintf("%.2f", -10*log($err_base)/log(10));
}
$accuracy = 1-$err_base;

open OUT, '>', $out_file or die "[$task] Can't open such file: $out_file.\n";
if ($QV == 0) {
	print OUT "No error base has been detected. QV can't be caculated.\n";
}else {
	print OUT "QV\t$QV\n";
}
close OUT;
##---------------------------------subroutine---------------------------------##

sub _system {
	my $cmd = shift;
	my $mode = shift if @_;
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
		if ( basename($path) eq "software" && -X $path ) {
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