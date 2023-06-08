#! /usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;
Getopt::Long::Configure qw(bundling no_ignore_case);

my $task = "run_hisat2";
my $sys_run = 0;

my ($assembly,$file_list,@reads1,@reads2,
    $dir,$pre_out,$prefix,$threads,
    $hisat2,$hisat2_build,$samtools);
GetOptions(
	"r:s"            => \$assembly,
	"i:s"            => \@reads1,
	"I:s"            => \@reads2,
	"l|list:s"       => \$file_list,
	"d:s"            => \$dir,
	"o:s"            => \$prefix,
	"t|threads:i"    => \$threads,
	"hisat2:s"       => \$hisat2,
	"hisat2_build:s" => \$hisat2_build,
	"samtools:s"     => \$samtools
);

##Check input files and directory.
die ("[$task]Error!Please input your assembly file.\n") unless $assembly;
$dir = "gaep_${task}_$$" unless $dir;
$dir =~ s/\/$//;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix = "${task}_output" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;

##-------------------------------variants---------------------------------##
my %flist = ();
my $index = "${prefix_out}_index";

##Check software.
$hisat2 = $hisat2 ? check_software("hisat2", $hisat2) : check_software("hisat2");
die "[$task] Error! Hisat2 is not found. 
[$task] Check your config file or set it in your environment.\n" if $hisat2 eq "-1";

$hisat2_build = $hisat2_build ? check_software("hisat2-build", $hisat2_build) : check_software("hisat2-build");
die "[$task] Error! Hisat2_build is not found. 
[$task] Check your config file or set it in your environment.\n" if $hisat2_build eq "-1";

$samtools = $samtools ? check_software("samtools", $samtools) : check_software("samtools");
die "[$task] Error! Samtools is not found. 
[$task] Check your config file or set it in your environment.\n" if $samtools eq "-1";

##reading file list
if ( $file_list ) {
	%flist = read_file_list($file_list);
	die "[$task] Empty file $file_list.\n" unless %flist; 
}

if ( @reads1 && @reads2 ) {
	if ( @reads1 != @reads2 ) {
		die( "[$task] Invalid input files.\n[$task] Different number of files in reads1 and reads2.\n" );
	}else {
		for ( 0..$#reads1 ) {
			push @{ $flist{paired} }, [$reads1[$_], $reads2[$_]];
		}
	}
}elsif ( @reads1 ) {
	push @{ $flist{single} }, @reads1;
}elsif ( @reads2 ) {
	push @{ $flist{single} }, @reads2;
}

my $out_file = "${prefix_out}_TRANS_mapping.bam";

##Creating a HISAT2 index
my $index_cmd = "$hisat2_build ";
$index_cmd .= "-p $threads " if $threads;
$index_cmd .= "$assembly $index";
_system($index_cmd,$sys_run);

if ( defined $flist{paired} ) {
	my $hisat2_paired_out = "$prefix_out.paired.sam";
	my $hisat2_paired_bam = "$prefix_out.paired.bam";
	#my $hisat2_paired_fixmate = "$prefix_out.paired.fixmate.bam";
	my $hisat2_paired_sort_bam = "$prefix_out.paired.sorted.bam";
	#my $hisat2_paired_markdup = "$prefix_out.paired.sorted.mkdup.bam";
	##run hisat2
	my $hisat2_cmd = "$hisat2 ";
	$hisat2_cmd .= "-p $threads " if $threads;
	$hisat2_cmd .= "-x $index ";
	for my $i ( 0..$#{ $flist{paired} } ) {
		my $hisat2_cmd = "$hisat2 ";
		$hisat2_cmd .= "-p $threads " if $threads;
		$hisat2_cmd .= "-x $index ";
		my $reads1 = $flist{paired}[$i][0];
		my $reads2 = $flist{paired}[$i][1];
		$hisat2_cmd .= "-1 $reads1 -2 $reads2 -S ${prefix_out}_paired_$i.sam";
		_system($hisat2_cmd,$sys_run);
	}
	
	##samtools merge
	if ( scalar @{ $flist{paired} }  > 1 ) {
		my $merge_cmd = "$samtools merge ";
		$merge_cmd .= "-@ $threads " if $threads;
		$merge_cmd .= "$hisat2_paired_out ";
		$merge_cmd .= "${prefix_out}_paired_${_}.sam " for ( 0..$#{ $flist{paired} } );
		_system($merge_cmd, $sys_run);
	}elsif ( scalar @{ $flist{paired} } == 1 ) {
		rename "${prefix_out}_paired_0.sam", $hisat2_paired_out or die "[$task] Can't rename ${prefix_out}_paired_0.sam to $hisat2_paired_out.\n";
	}
	##samtools view
	my $view_cmd = "$samtools view ";
	$view_cmd .= "-@ $threads " if $threads;
	$view_cmd .= "-bS $hisat2_paired_out > $hisat2_paired_bam";
	_system($view_cmd, $sys_run);
	##samtools sort
	my $sort_cmd = "$samtools sort ";
	$sort_cmd .= "-@ $threads " if $threads;
	$sort_cmd .= "$hisat2_paired_bam -o $hisat2_paired_sort_bam";
	_system ($sort_cmd,$sys_run);
	
	rename $hisat2_paired_sort_bam, $out_file or die "[$task] Can't rename $hisat2_paired_sort_bam to $out_file.\n";
}

if ( defined $flist{single} ) {
	my $hisat2_single_out = "$prefix_out.single.sam";
	my $hisat2_single_bam = "$prefix_out.single.bam";
	#my $hisat2_single_fixmate = "$prefix_out.single.fixmate.bam";
	my $hisat2_single_sort_bam = "$prefix_out.single.sorted.bam";
	#my $hisat2_single_markdup = "$prefix_out.single.sorted.mkdup.bam";
	
	##run hisat2
	for my $i ( 0..$#{ $flist{single} } ){
		my $hisat2_cmd = "$hisat2 ";
		$hisat2_cmd .= "-p $threads " if $threads;
		$hisat2_cmd .= "-x $index ";
		my $reads1 = $flist{single}[0];
		$hisat2_cmd .= "-U $reads1 -S ${prefix_out}_single_$i.sam";
		_system($hisat2_cmd,$sys_run);
	}
	
	##samtools merge
	if ( scalar @{ $flist{single} }  > 1 ) {
		my $merge_cmd = "$samtools merge ";
		$merge_cmd .= "-@ $threads " if $threads;
		$merge_cmd .= "$hisat2_single_out ";
		$merge_cmd .= "${prefix_out}_single_${_}.sam " for ( 0..$#{ $flist{single} } );
		_system($merge_cmd, $sys_run);
	}elsif ( scalar @{ $flist{single} } == 1 ) {
		rename "${prefix_out}_single_0.sam", $hisat2_single_out or die "[$task] Can't rename ${prefix_out}_single_0.sam to $hisat2_single_out.\n";
	}
	##samtools view
	my $view_cmd = "$samtools view ";
	$view_cmd .= "-@ $threads " if $threads;
	$view_cmd .= "-bS $hisat2_single_out > $hisat2_single_bam";
	_system($view_cmd, $sys_run);
	##samtools sort
	my $sort_cmd = "$samtools sort ";
	$sort_cmd .= "-@ $threads " if $threads;
	$sort_cmd .= "$hisat2_single_bam -o $hisat2_single_sort_bam";
	_system ($sort_cmd,$sys_run);
	rename $hisat2_single_sort_bam, $out_file or die "[$task] Can't rename $hisat2_single_sort_bam to $out_file.\n";
}


##samtools flagstats
#my $flagstats = "samtools flagstats ";
#$flagstats .= "-@ $threads " if $threads;
#$flagstats .= "$prefix_out.merged.bam";
#$mapping_rate = `$flagstats | perl -lne 'print \$1 if /mapped\\s+\\((.*%)/'`;
#
#printf ("Mapping rate:\t%d", $mapping_rate0);

##---------------------------------subroutine---------------------------------##

sub read_file_list {
	my $flist = shift;
	my %file_list;
	open FLIST, '<', "$flist" or die "Can't open file: $flist.\n";
	while (<FLIST>) {
		chomp;
		my @line = split;
		if ( @line == 1 ) {
			push @{ $file_list{single} }, $_;
		}elsif ( @line == 2 ) {
			push @{ $file_list{paired} }, [@line];
		}else {
			die "[$task] Found more than two elements in $flist:$..\n";
		}
	}
	close FLIST;
	return %file_list;
}

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