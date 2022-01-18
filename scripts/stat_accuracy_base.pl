#! /usr/bin/perl -w

use strict;
use Getopt::Long;

my $task = "stat_accuracy_base";
my $mode = 0;

my $usage = qq(
Program: $task
Version: v1.0

Usage: $task 
Options:
         -v <FILE>           VCF/BCF files
         -b <FILE>           bam file
         -d <INT>            output directory
         -o <INT>            output prefix
         --bcftools <PATH>   path to bcftools
         --bedtools <PATH>   path to bedtools
         -h                  print help massage
);

##----------------------------------options-----------------------------------##

my ($vcf,$bam,$dir,$prefix_out,$bcftools,$bedtools,$help);
GetOptions (
	"v:s"        => \$vcf,
	"b:s"        => \$bam,
	"d:s"        => \$dir,
	"o:s"        => \$prefix_out,
	"h"          => \$help,
	"bcftools:s" => \$bcftools, 
	"bedtools:s" => \$bedtools,
);
die ($usage) if $help;
##Check input files.
##Check input files and directory.
die ($usage) unless $vcf and $bam;
print STDERR "[$task] No such file: $vcf." unless -e $vcf;
print STDERR "[$task] No such file: $bam." unless -e $bam;
$prefix_out = "${task}_output" unless $prefix_out;
if ( $dir ){
	unless ( -e $dir ) {
		die ("[$task]Error! Can't make directory:\"$dir\"\n") if ( system "mkdir $dir" );
	}else {
		$dir =~ s/\/$//;
	}
	$prefix_out = "$dir/$prefix_out";
}

##Check software.
$bcftools = $bcftools ? check_software("bcftools", $bcftools) : check_software("bcftools");
die "[$task] Error! Bcftools is not found.
[$task] Check your config file or set it in your environment.\n" if $bcftools eq "-1";
$bedtools = $bedtools ? check_software("bedtools", $bedtools) : check_software("bedtools");
die "[$task] Error! Bedtools is not found.
[$task] Check your config file or set it in your environment.\n" if $bedtools eq "-1";


##---------------------------------variants-----------------------------------##

my $mean_depth = 0;
my $total_len  = 0;
my $map_len = 0;
my $var_num = 0;
my $accuracy = 0;
my $QV = 0;

my $out_file = "$prefix_out.base.qv";
##-----------------------------------main-------------------------------------##

##mappable length
my $genomecov = "$prefix_out.genomecov";
my $bedtools_cmd = "$bedtools genomecov -ibam $bam -split > $genomecov";
#_system($bedtools_cmd, $mode);
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

$accuracy = $var_num/$map_len;
$QV = sprintf("%.2f", -10*log($accuracy)/log(10));
$accuracy = 1-$accuracy;

open OUT, '>', $out_file or die "[$task] Can't open such file: $out_file.\n";
print OUT "QV\t$QV\n";
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