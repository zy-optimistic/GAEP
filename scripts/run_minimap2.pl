#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $task = "run_minimap2";
my $sys_run = 0;

my $usage = qq (
Program: $task
Version: v1.0

Usage : $task -r <assembly> (-i <fastq> | -l <fastq list>) [options]

Options:
         -r <FILE>    Genome assembly in a fasta file.
         -i <FILE>    Fastq file. Can be set more than one time.
         -l <FILE>    File listing fastq files. One per line.
         -x <STR>     Reads type. (ont,pb)
         -t <INT>     Threads use in minimap2 and samtools.
         -o <STR>     Output prefix.
         -d <DIR>     Output directory.
         --minimap2   Path to minimap2.
         --samtools   Path to samtools.
); 

my ($assembly,$file_list,$reads_type,$threads,$dir,$prefix,$help,$minimap2,$samtools);
my @sequences = ();
GetOptions(
	"r:s"        => \$assembly,
	"i:s"        => \@sequences,     #FASTQ or FASTA
	"l|list:s"   => \$file_list,
	"x:s"        => \$reads_type,	   #(ont,pb)
	"t:i"        => \$threads,
	"d:s"        => \$dir,
	"o:s"        => \$prefix,
	"h"          => \$help,
	"minimap2:s" => \$minimap2,
	"samtools:s" => \$samtools,
);

##Initial setting
unless ($assembly && (@sequences || $file_list)) {
	die $usage;
}

$dir = "gaep_${task}_$$" unless $dir;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix = "tgs_mapping_$$" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;

print "\n[$task] Output directory is $dir.\n";

my $out_bam = "${prefix_out}_TGS_mapping.bam";
my @file_list = ();

##Check software.
$minimap2 = $minimap2 ? check_software("minimap2", $minimap2) : check_software("minimap2");
die "[$task] Error! Minimap2 is not found.
[$task] Check your config file or set it in your environment.\n" if $minimap2 eq "-1";

$samtools = $samtools ? check_software("samtools", $samtools) : check_software("samtools");
die "[$task] Error! Samtools is not found.
[$task] Check your config file or set it in your environment.\n" if $samtools eq "-1";

##reading file list
if (@sequences){
	@file_list = @sequences;
}
if ($file_list){
	open CONFIG , $file_list;
	while (<CONFIG>){
		chomp;
		push @file_list , $_;
	}
}

my $minimap2_cmd = "$minimap2 --secondary=no ";
$minimap2_cmd .= "-t $threads " if $threads;
my $temp = $minimap2_cmd;
my $count = 0;
my @sam_list;

foreach (@file_list){
	$count ++;
	#run minimap2
	if (! $reads_type) {
		$minimap2_cmd .= "-L -a $assembly $_ -o ${prefix_out}_${count}.sam"; 
	}elsif ($reads_type	eq 'ont'){
		$minimap2_cmd .= "-L -ax map-ont $assembly $_ -o ${prefix_out}_${count}.sam"; 
	}elsif ($reads_type	eq 'pb'){
		$minimap2_cmd .= "-L -ax map-pb $assembly $_ -o ${prefix_out}_${count}.sam"; 
	}elsif ($reads_type	eq 'ccs'){
		$minimap2_cmd .= "-L -ax asm20 $assembly $_ -o ${prefix_out}_${count}.sam"; 
	}else {
		$minimap2_cmd .= "$assembly $_ -o ${prefix_out}_${count}.sam";
	}
	_system($minimap2_cmd, $sys_run);
	push @sam_list, "${prefix_out}_${count}";
	$minimap2_cmd = $temp;
}

##sort
my $sort_cmd = "samtools sort -O BAM ";
$sort_cmd .= "-@ $threads " if $threads;
$temp = $sort_cmd;
foreach (@sam_list){
	$sort_cmd .= "$_.sam -o $_.sorted.bam";
	_system($sort_cmd, $sys_run);
	$sort_cmd = $temp;
}

##merge
if (@file_list > 1){
	my $merge_cmd = "samtools merge ";
	$merge_cmd .= "-@ $threads " if $threads;
	$merge_cmd .= "-h $sam_list[0].sorted.bam $out_bam ";
	$merge_cmd .= join (" ",map {"$_.sorted.bam"} @sam_list);
	_system($merge_cmd, $sys_run);
}else {
	rename "${prefix_out}_1.sorted.bam", $out_bam or die "[$task] Can't rename ${prefix_out}_1.sorted.bam.\n";
}

##-------------------------------subroutines----------------------------------##

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
