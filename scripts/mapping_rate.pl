#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;

my $task = "mapping_rate";

my ($bam, $threads, $out_file, $samtools);
GetOptions(
	"b:s"            => \$bam,
	"o:s"            => \$out_file,
	"t:i"            => \$threads,
	"--samtools:s"   => \$samtools,
);
die "No bam file was input.\n" if !$bam;

$samtools = $samtools ? check_software("samtools", $samtools) : check_software("samtools");
die "[$task] Error! Samtools was not found.
[$task] Check your config file or set it in your environment.\n" if $samtools eq "-1";

$out_file = "$bam.flagstat_$$.txt" if !$out_file;

my $cmd = "$samtools flagstat ";
$cmd .= "-@ $threads " if $threads;
$cmd .= "$bam ";
$cmd .= "> $out_file";

if (system $cmd) {
	die "[$task] Running error with command line:\n $cmd";
}

##-------------------------------subroutines----------------------------------##

sub check_software {
	my $software = shift;
	my $path = shift if @_;
	if ( $path ) {
		if ( basename($path) eq $software && -X $path ) {
			return $software = $path;
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