#! /usr/bin/perl -w 

=head1 Name

  cov_snp_dot.pl

=head1 Description

  Generate a dotplot about reads coverage and SNP density of windows of a genome
  sequences.

=head1 Usage
 
    gaep snvcov <-r assembly.fasta> -i sr1.fastq -I sr2.fastq -l lr.fastq -x (ont, pb, ccs) [options]
  	-r         <FILE>    genome assembly in a fasta file
	-i         <FILE>    NGS fastq1 file. Can be set more than one time
	-I         <FILE>    NGS fastq2 file. Can be set more than one time
	--nl       <FILE>    file listing NGS fastq files, Read1 and read2 delimited with tab
	-l         <FILE>    TGS fastq file. Can be set more than one time
	--tl       <FILE>    file listing TGS fastq files
	-x         <STR>     TGS reads type (ont, pb, ccs)
	-t         <INT>     number of threads
	-m         <INT>     contig short than INT will be ignored
	-s         <INT>     windows size
	-f         <STR>     output plot format (pdf, png)
	-o         <STR>     out prefix
	-d         <STR>     out directory
	--bedtools <PATH>    path to bedtools
	--bcftools <PATH>    path to bcftools
	--Rscript  <PATH>    path to Rscript
	-h                   print help information

=head1 Usage
 
  GAEP uses the short reads to called SNVs between reads and contigs, while the 
  long reads are utilized for calculating the read coverage. Both types of reads
  will be aligned to the assembly.
  
=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use FindBin qw($RealBin);
Getopt::Long::Configure qw( bundling no_ignore_case );

my $task = "cov_snp_dot";
my $mode = 0;

my ($assembly, @reads1, @reads2, $ngs_list,
	@long_reads, $tgs_list, $reads_type, $ngs_bam,
	$threads, $fwindows, $bam, $vcf, $min_ctg, 
	$prefix, $dir, $bedtools, $bcftools, $win_size,
	$Rscript, $help);
	
my $format = 'pdf';
GetOptions(
	"r:s"        => \$assembly,
	"i:s"        => \@reads1,
	"I:s"        => \@reads2,
	"nl:s"       => \$ngs_list,
	"l:s"        => \@long_reads,
	"tl:s"       => \$tgs_list,
	"x:s"        => \$reads_type,	   #(ont,pb)
	"t:i"        => \$threads,
	"b:s"        => \$bam,
	"n:s"        => \$ngs_bam,
	"v:s"        => \$vcf,
	"m:i"        => \$min_ctg,
	"w:s"        => \$fwindows,
	"s:i"        => \$win_size,
	"f:s"        => \$format,
	"o:s"        => \$prefix,
	"d:s"        => \$dir,
	"bedtools:s" => \$bedtools,
	"bcftools:s" => \$bcftools,
	"Rscript:s"  => \$Rscript,
	"h"          => \$help
);

##Initial setting
die `pod2text $0` if $help;
if (!$vcf && !$bam) {
	if (!$assembly) {
		die "[$task] Error! Can't find the assembly file.\n", `pod2text $0`;
	}
	if (!$vcf) {
		if (!$ngs_bam) {
			if (!@reads1 && !@reads2 && !$ngs_list) {
				print STDERR "[$task] Error! No vcf file!";
				die `pod2text $0`;
			}
		}
	}
	if (!$bam) {
		if (!$ngs_bam) {
			if (!@reads1 && !@reads2 && !$ngs_list) {
				if (!@long_reads && !$tgs_list) {
					print STDERR "[$task] Error! No bam file!";
					die `pod2text $0`;	
				}
			}
		}
			
	}
}
$dir = "gaep_${task}_$$" unless $dir;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix = "snv_depth_$$" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;
print "\n[$task] Output directory is $dir.\n";

my $hist = "$prefix_out.${task}_$$.hist";

##Check software.
$bcftools = $bcftools ? check_software("bcftools", $bcftools) : check_software("bcftools");
die "[$task] Error! bcftools is not found.
[$task] Check your config file or set it in your environment.\n" if $bcftools eq "-1";
$bedtools = $bedtools ? check_software("bedtools", $bedtools) : check_software("bedtools");
die "[$task] Error! bedtools is not found.
[$task] Check your config file or set it in your environment.\n" if $bedtools eq "-1";

if (!$vcf) {
	if (!$ngs_bam) {
		system "mkdir -p $dir/mapping/NGS_mapping/";
		my $map_cmd = "perl $RealBin/run_bwa.pl -r $assembly ";
		for (my $i = 0; $i <= $#reads1 && $i <= $#reads2; $i++) {
			$map_cmd .= "-i $reads1[$i] -I $reads2[$i] ";
		}
		$map_cmd .= "-l $ngs_list " if $ngs_list;
		$map_cmd .= "-t $threads "   if $threads;
		$map_cmd .= "-d $dir/mapping/NGS_mapping ";
		$map_cmd .= "-o $prefix ";
		_system($map_cmd, $mode);
		$ngs_bam = "$dir/mapping/NGS_mapping/$prefix.paired.sorted.bam";
	}
	system "mkdir -p $dir/variants/";
	my $var_cmd = "perl $RealBin/var_calling.pl -r $assembly ";
	$var_cmd .= "-b $ngs_bam ";
	$var_cmd .= "-t $threads "   if $threads;
	$var_cmd .= "-d $dir/variants/ ";
	$var_cmd .= "-o $prefix ";
	_system($var_cmd, $mode);
	$vcf = "$dir/variants/$prefix.flt.vcf";
}

if (!$bam) {
	system "mkdir -p $dir/mapping/TGS_mapping/";
	if (@long_reads) {
		my $tgs_cmd = "perl $RealBin/run_minimap2.pl -r $assembly ";
		$tgs_cmd .= "-i $_ " for @long_reads;
		$tgs_cmd .= "-x $reads_type " if $reads_type;
		$tgs_cmd .= "-t $threads "   if $threads;
		$tgs_cmd .= "-d $dir/mapping/TGS_mapping/ ";
		$tgs_cmd .= "-o $prefix ";
		_system($tgs_cmd, $mode);
		$bam = "$dir/mapping/TGS_mapping/${prefix}_TGS_mapping.bam";
	}elsif (-e $ngs_bam) {
		$bam = $ngs_bam;
	}else {
		print "[$task] Can't find a bam file.\n";
		die `pod2text $0`;
	}
}

##--get length--##
$min_ctg = $min_ctg ? $min_ctg : 0;
my $seq_len = read_header($bam, $min_ctg);
my $tolen = 0;
$tolen += $_->[1] for @$seq_len;
if (!$win_size) {
	if ($tolen < 500000000) {
		$win_size = int($tolen/5000);
	}elsif ($tolen > 4000000000) {
		$win_size = 200000;
	}else {
		$win_size = 100000;
	}
}

print STDERR "[$task] Make windows in length $win_size.\n";
##--make windiows--##
if ( ! $fwindows ) {
	$fwindows = "$dir/${prefix}_${win_size}_windows.txt";
	open WIN, '>', $fwindows || die "Can't open such file: $fwindows.\n";
	for my $a ( @{$seq_len} ) {
		my $start = 0;
		my $end = 0;
		for  ( my $i = 1; $i <= int($a->[1]/$win_size); $i++ ) {
			$end = $i * $win_size;
			print WIN join("\t",$a->[0],$start,$end),"\n";
			$start = $end;
		}
		print WIN join("\t",$a->[0],$start,$a->[1]),"\n" if ( $end < $a->[1] );
	}
	close WIN;
}
##--check windiows--##
##--read windows--##
my %windows = ();
open WIN, '<', $fwindows || die "Can't open such file: $fwindows.\n";
while (<WIN>) {
	chomp;
	my @line = split;
	$windows{$line[0]}{$line[1]} = 0;
}
close WIN;
#print Dumper(\%windows);

##--SNP density--##
window2SNPden(\%windows);

##--Coverage--##
my $cmd = "$bedtools coverage -mean -sorted -bed -nobuf -split ";
   $cmd .= "-a $fwindows ";
   $cmd .= "-b $bam";
print STDERR $cmd;
open COV, '-|', $cmd || die "[$task] Fail to open pipe: $cmd.\n";
open OUT, '>', $hist || die "[$task] Fail to open pipe: $hist.\n";
while (<COV>) {
	chomp;
	my @line = split;
	print OUT "$_\t", $windows{$line[0]}{$line[1]}/($line[2]-$line[1])*100, "\n";
}
close COV;
close OUT;

##plot
$Rscript = $Rscript ? check_software("Rscript", $Rscript) : check_software("Rscript");
die "[$task] Error! Rscript is not found.
[$task] Check your config file or set it in your environment.\n" if $Rscript eq "-1";
my $plot_cmd = "$Rscript $RealBin/snp_coverage_dotplot.R $hist $prefix_out $format";
_system($plot_cmd, $mode);

#---------------------------------Subroutine-----------------------------------#
sub read_header {
        my $bam = shift;
        my $length = [];
        my $pip = "samtools view -H $bam ";
        open my $BAM, '-|', $pip || die "Fail to read header of $bam.\n";
        while (<$BAM>) {
            if (/^\@SQ\s+SN:(\S+)\s+LN:(\d+)/){
				push @$length, [$1,$2] if ($2 >= $min_ctg);
				
			}
        }
        die "Fail to read header of $bam.\n" unless (@{$length});
        return $length;
}

sub window2SNPden {
	my $win = shift;
	open VCF, '-|', "$bcftools view -i 'GT=\"het\" && TYPE=\"snp\"' $vcf" || die "Fail to open pipe: $bcftools view $vcf.\n";
	while (<VCF>) {
		next if ( /^#/ );
		my ($chr,$posi) = (split/\s+/,$_,3)[0,1];
		my $index = (int($posi/$win_size)) * $win_size;
		#print join("\t",$posi$win_size,$win_size),"\n";
		#print $index,"\n";
		$win->{$chr}{$index} += 1;
	}
	close VCF;
	#print Dumper($win);
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
