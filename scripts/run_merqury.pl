#! /usr/bin/perl -w

use strict;
use File::Basename;
use FindBin qw($RealBin);
use Getopt::Long;

my $task = "run_merqury";
my $mode = 0;

my $usage = qq(
Program: $task
Version: v1.0

Usage:   $task <-r fasta> <-i fastq> [<-i fastq>] [options]
Options:
         -r        genome assembly in a fasta file
         -g        total length of assembly
         -i        fastq file. Can be set more than one time
         -t        max threads to use
         -d        output directory
         -o        output prefix
         -h        print help massage
         --meryl:s"    => \$meryl, 
         --merqury:s"  => \$merqury,
);

##----------------------------------options-----------------------------------##

my ($assembly, @reads, $threads,  $dir, $prefix, $meryl, $merqury, $help);
my $genome_len = 0;
GetOptions (
    "r:s"        => \$assembly,
    "g:i"        => \$genome_len,
    "i:s"        => \@reads,
    "t:i"        => \$threads,
    "d:s"        => \$dir,
    "o:s"        => \$prefix,
    "h"          => \$help,
    "meryl:s"    => \$meryl, 
    "merqury:s"  => \$merqury,
);
die $usage if !$assembly || !@reads;
die $usage if $help;

$dir = "gaep_${task}_$$" unless $dir;
if (! -e $dir){
    if (system "mkdir -p $dir"){
        die "[$task] Error! Can't make directory:\"$dir\"\n";
    }
}else {
	my @f = <$dir>;
	if (@f) {
			die "[$task] Directory $dir is not empty. Please specify another directory or a empty directory.\n";
	}
}
if (! -e "$dir/meryl"){
    if (system "mkdir -p $dir/meryl"){
        die "[$task] Error! Can't make directory:\"$dir/meryl\"\n";
    }
}

$prefix = "merqury_$$" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;
print "\n[$task] Output directory is $dir.\n";

##---------------------------------variants-----------------------------------##
my $qv = "$prefix_out.qv";

##-----------------------------------main-------------------------------------##
##Check software.
$meryl = $meryl ? check_software("meryl", $meryl) : check_software("meryl");
die "[$task] Error! Meryl is not found.
[$task] Check your config file or set it in your environment.\n" if $meryl eq "-1";
$merqury = $merqury ? check_software("merqury.sh", $merqury) : check_software("merqury.sh");
die "[$task] Error! Merqury is not found.
[$task] Check your config file or set it in your environment.\n" if $merqury eq "-1";

##calculate genome length
if (! $genome_len) {
    open GENO, '<', $assembly or die "[$task] Can't open $assembly for reading.\n";
    while (<GENO>) {
       if (/^>/) {
           next;
       }else {
           chomp;
           $genome_len += length $_;
       }
    }
    close GENO;
}

##calculate k length
my $best_k = "$RealBin/../third_party/merqury/best_k.sh";
die "[$task] Can't find best_k.sh at $RealBin/../third_party/merqury/.\n" if (! -e $best_k);
my $klen = `$best_k $genome_len`;
$klen = int $klen;
$klen += 1 if $klen % 2 == 0;

##set threads
$ENV{OMP_NUM_THREADS} = $threads;

##count K-mer in reads
for my $read (@reads) {
    my $meryl_cmd = "$meryl k=$klen count output $dir/meryl $read";
    _system($meryl_cmd, $mode);
}

##change directory and create symbolic link.
my $assembly_a = `pwd`;
$assembly_a =~ s/\n?$//;
$assembly_a =~ s/\/?$//;
$assembly_a = $assembly_a."/".$assembly;
my $assembly_s = basename($assembly);
chdir("$dir") or die "[$task] Can't open directory: $dir.\n";
symlink $assembly_a, $assembly_s or die "[$task] Can't create symbolic link: $assembly_s.\n"; 
##merge
if (@reads > 1) {
    my $merge_cmd = "$meryl union-sum output $prefix.meryl meryl";
    _system($merge_cmd, $mode);
}

##run merqury
my $merqury_cmd = "$merqury $prefix.meryl $assembly_s $prefix";
_system($merqury_cmd, $mode);
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