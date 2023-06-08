#! /usr/bin/perl -w

use strict;
use File::Basename;
use Getopt::Long;
use threads;
use Thread::Semaphore;
use Thread::Queue;

my $task = "var_calling";
my $mode = 0;

my $threads = 4;
my ($assembly,$bam,$dir,$prefix_out,$markduped,$bcftools,$samtools);
GetOptions (
	"r:s"           => \$assembly,
	"b:s"           => \$bam,
	"t:i"           => \$threads,
	"d:s"           => \$dir,
	"o:s"           => \$prefix_out,
	"skip_markdup"  => \$markduped, #bam file has been marked duplications
	"bcftools:s"    => \$bcftools, 
	"samtools:s"    => \$samtools
);

my $usage = qq/[$task]
usage:$task -r <file.fa> -i <file.bam> 
            -t INT	threads 
/;
##Check input files.
##Check input files and directory.
die ($usage) unless $assembly and $bam;
print STDERR "[$task] No such file: $assembly." unless -e $assembly;
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

##test software
$bcftools = $bcftools ? check_software("bcftools", $bcftools) : check_software("bcftools");
die "[$task] Error! Bcftools is not found. 
[$task] Check your config file or set it in your environment.\n" if $bcftools eq "-1";
$samtools = $samtools ? check_software("samtools", $samtools) : check_software("samtools");
die "[$task] Error! Samtools is not found. 
[$task] Check your config file or set it in your environment.\n" if $samtools eq "-1";

##---------------------------------variants-----------------------------------##
my $markdup  = "$prefix_out.markdup.bam";
my $var_result  = "$prefix_out.flt.vcf";
my %seq_len = read_header($bam);
my $ctg_queue = Thread::Queue->new();
my $semaphore = Thread::Semaphore->new($threads);

#bam markdup
if (! $markduped) {
	print STDERR "[$task] Run samtools markdup to mark duplications.\n";
	my $markdup_cmd = "$samtools markdup ";
	$markdup_cmd .= "-@ $threads " if $threads;
	$markdup_cmd .= "$bam $markdup";
	_system($markdup_cmd, $mode);
	$bam = $markdup;
	#index
	my $index_cmd = "$samtools index ";
	$index_cmd .= "-@ $threads " if $threads;
	$index_cmd .= "$bam";
	_system($index_cmd, $mode);
}

if (!-s "${bam}.bai" && !-s "${bam}.csi"){
	my $cmd = "$samtools index ";
	$cmd   .= "-@ $threads " if $threads;
	$cmd   .= "$bam";
	_system($cmd, $mode);
}

for (sort keys %seq_len) {
	$ctg_queue->enqueue($_);
}
$ctg_queue->end();

##varants calling
for (my $i = 0; $i < $threads; $i++) {
	$semaphore->down(1);
	threads->create(\&process) -> detach();
}
$semaphore->down($threads);
$semaphore->up($threads);

##merge all vcfs into one
my $concat = "$bcftools concat $prefix_out.*.flt.vcf > $var_result";
_system($concat, $mode);


##---------------------------------subroutine---------------------------------##

sub process {
	while (defined(my $ctg = $ctg_queue->dequeue())) {
		snv_calling($ctg);
	}
	$semaphore->up(1);
}

sub snv_calling {
	my $ctg = shift;
	my $var_call = "$prefix_out.$ctg.vcf";
	my $var_flt  = "$prefix_out.$ctg.flt.vcf";
	
	##snp calling
	my $call_snp = "$bcftools mpileup -r $ctg ";
	$call_snp .= " -a SP -f $assembly $bam | $bcftools call ";
	$call_snp .= "-vm -f GQ -o $var_call";
	_system($call_snp, $mode);
	##filter
	my $filter_cmd = "$bcftools filter ";
	$filter_cmd .= " -e '%QUAL<10 || DP <= 3' $var_call -o $var_flt";
	_system($filter_cmd, $mode);
}

sub read_header {
	my $bam = shift;
	my %length = ();
	my $pip = "$samtools view -H $bam ";
	open my $BAM, '-|', $pip || die "[$task 1]Fail to read header of $bam.\n";
	while (<$BAM>) {
		/^\@SQ\s+SN:(\S+)\s+LN:(\d+)/ ? $length{$1} = $2 : next;
	}
	die "[$task 2]Fail to read header of $bam.\n" unless (%length);
	return %length;
}

sub _system {
	my $cmd = shift;
	my $mode = shift if @_;
	if ( $mode == 1 ) {
		print $cmd, "\n";
	}else {
		print $cmd, "\n";
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
