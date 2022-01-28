#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

Getopt::Long::Configure qw( bundling no_ignore_case );

my $task = "run_bwa";
my $sys_run = 0;

my ($assembly,@reads1,@reads2,$threads,$dir,$prefix,$index,$file_list,$bwa,$samtools);
GetOptions(
	"r:s"        => \$assembly,
	"i:s"        => \@reads1,
	"I:s"        => \@reads2,
	"l|list:s"   => \$file_list,
	"t:i"        => \$threads,
	"d:s"        => \$dir,
	"o|output:s" => \$prefix,    #prefix of output files
	"bwa:s"      => \$bwa,
	"samtools:s" => \$samtools
);

##Check input files and directory.
die ("[$task]Error!Please input your assembly file.\n") unless $assembly;

$dir = "gaap_${task}_$$" unless $dir;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix = "ngs_mapping_$$" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;

##check software
$bwa = $bwa ? check_software("bwa", $bwa) : check_software("bwa");
die "[$task] Error! Bwa is not found. 
[$task] Check your config file or set it in your environment.\n" if $bwa eq "-1";
$samtools = $samtools ? check_software("samtools", $samtools) : check_software("samtools");
die "[$task] Error! Samtools is not found. 
[$task] Check your config file or set it in your environment.\n" if $samtools eq "-1";

##-------------------------------variants---------------------------------##
my %flist = ();


##Make a file list
if ( $file_list ){
	%flist = read_file_list($file_list);
	die "[$task] Empty file $file_list.\n" unless %flist; 
}

if ( @reads1 && @reads2 ){
	if ( @reads1 != @reads2 ){
		die( "[$task] Invalid input files.\n[$task] Different number of files in reads1 and reads2.\n" );
	}else {
		for ( 0..$#reads1 ){
			push @{ $flist{paired} }, [$reads1[$_], $reads2[$_]];
		}
	}
}elsif ( @reads1 ) {
	push @{ $flist{single} }, @reads1;
}elsif ( @reads2 ) {
	push @{ $flist{single} }, @reads2;
}
#print $#{$flist{paired}};
#print $#{$flist{single}};
##Run bwa to align reads to assembly.
#index
unless ( -f "$assembly.bwt" && -f "$assembly.pac" && -f "$assembly.sa" && -f "$assembly.ann" && -f "$assembly.amb" ) {
	print STDERR "[$task] Failed to locate the BWA index. Run bwa index.\n";
	my $fsize = -s $assembly;
	my $index_cmd = "$bwa index ";
	
	$index_cmd .= "$assembly";
	_system($index_cmd,$sys_run);
}

if ( defined $flist{paired} ) {
	my $bwa_out = "$prefix_out.paired.sam";
	my $bam = "$prefix_out.paired.bam";
	my $fixmate_out = "$prefix_out.paired.fixmate.bam";
	my $sort_bam = "$prefix_out.paired.sorted.bam";
	#my $markdup = "$prefix_out.paired.sorted.mkdup.bam";
	#run bwa mem
	for my $i ( 0..$#{ $flist{paired} } ){
		my $reads1 = $flist{paired}[$i][0];
		my $reads2 = $flist{paired}[$i][1];
		my $bwa_mem_cmd = "$bwa mem -M ";
		$bwa_mem_cmd .= "-t $threads " if $threads;
		$bwa_mem_cmd .= "$assembly ";
		$bwa_mem_cmd .= "$reads1 $reads2 -o ${prefix_out}_paired_${i}.sam";
		_system($bwa_mem_cmd, $sys_run);
	}
	
	
	##Run samtools
	#samtools merge
	if ( scalar @{ $flist{paired} }  > 1 ) {
		my $merge_cmd = "$samtools merge ";
		$merge_cmd .= "-@ $threads " if $threads;
		$merge_cmd .= "$bwa_out ";
		$merge_cmd .= "${prefix_out}_paired_${_}.sam " for ( 0..$#{ $flist{paired} } );
		_system($merge_cmd, $sys_run);
	}elsif ( scalar @{ $flist{paired} } == 1 ) {
		rename "${prefix_out}_paired_0.sam", $bwa_out or die "[$task] Can't rename ${prefix_out}_paired_0.sam to $bwa_out.\n";
	}
	
	#sam -> bam 
	##my $view_cmd = "$samtools view ";
	##$view_cmd .= "-@ $threads " if $threads;
	##$view_cmd .= "-bS $bwa_out > $bam";
	##_system($view_cmd, $sys_run);
	#Run samtools fixmate to mark duplications later.
	my $fixmate_cmd = "$samtools fixmate ";
	$fixmate_cmd .= "-@ $threads " if $threads;
	$fixmate_cmd .= "-m $bwa_out $fixmate_out";
	_system($fixmate_cmd, $sys_run);
	#sort bam
	my $sort_cmd = "$samtools sort ";
	$sort_cmd .= "-@ $threads " if $threads;
	$sort_cmd .= "$fixmate_out -o $sort_bam";
	_system($sort_cmd, $sys_run);
	##bam markdup 
	#my $markdup_cmd = "$samtools markdup ";
	#$markdup_cmd .= "-@ $threads " if $threads;
	#$markdup_cmd .= "$sort_bam $markdup";
	#_system($markdup_cmd, $sys_run);
}

if ( defined $flist{single} ) {
	my $bwa_out = "$prefix_out.single.sam";
	my $bam = "$prefix_out.single.bam";
	my $fixmate_out = "$prefix_out.single.fixmate.bam";
	my $sort_bam = "$prefix_out.single.sorted.bam";
	#my $markdup = "$prefix_out.single.sorted.mkdup.bam";
	
	#run bwa mem
	for my $i ( 0..$#{ $flist{single} } ){
		my $reads1 = $flist{single}[$i];
		my $bwa_mem_cmd = "$bwa mem -M ";
		$bwa_mem_cmd .= "-t $threads " if $threads;
		$bwa_mem_cmd .= "$assembly ";
		$bwa_mem_cmd .= "$reads1 > ${prefix_out}_single_${i}.sam";
		_system($bwa_mem_cmd, $sys_run);
	}
	
	##Run samtools
	#samtools merge
	if ( scalar @{ $flist{single} }  > 1 ) {
		my $merge_cmd = "$samtools merge ";
		$merge_cmd .= "-@ $threads " if $threads;
		$merge_cmd .= "$bwa_out ";
		$merge_cmd .= "${prefix_out}_single_${_}.sam " for ( 0..$#{ $flist{single} } );
		_system($merge_cmd, $sys_run);
	}elsif ( scalar @{ $flist{single} } == 1 ) {
		rename "${prefix_out}_single_0.sam", $bwa_out or die "[$task] Can't rename ${prefix_out}_single_0.sam to $bwa_out.\n";
	}
	#sam -> bam 
	my $view_cmd = "$samtools view ";
	$view_cmd .= "-@ $threads " if $threads;
	$view_cmd .= "-bS $bwa_out > $bam";
	_system($view_cmd, $sys_run);
	#Run samtools fixmate to mark duplications later.
	my $fixmate_cmd = "$samtools fixmate ";
	$fixmate_cmd .= "-@ $threads " if $threads;
	$fixmate_cmd .= "-m $bam $fixmate_out";
	_system($fixmate_cmd, $sys_run);
	#sort bam
	my $sort_cmd = "$samtools sort ";
	$sort_cmd .= "-@ $threads " if $threads;
	$sort_cmd .= "$fixmate_out -o $sort_bam";
	_system($sort_cmd, $sys_run);
	#bam markdup 
	#my $markdup_cmd = "$samtools markdup ";
	#$markdup_cmd .= "-@ $threads " if $threads;
	#$markdup_cmd .= "$sort_bam $markdup";
	#_system($markdup_cmd, $sys_run);
}

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
		if ( system $cmd ) {
			die "[$task] Can't run \"$cmd\".\n";
		} 
	}
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

