#! /usr/bin/perl -w

=head1 Description

 GAAP: Genome assembly assessment pipeline
 Version: V1.0

=head1 Usage

 gaap pipe <-r fasta_file> [options]

=head1 Options

=head2 Necessary

 -r assembly
 
=head2 TGS data

 --lr         TGS long reads, can be specified more than once.
 -x/--type    Type of TGS long reads(pb, ont, ccs), used for mapping.
 
=head2 NGS data

 --sr1        Fisrt mate of paired-end NGS reads, can be specified more than once.
 --sr2        Second mate of paired-end NGS reads, can be specified more than once.

=head2 TRANS data
	
 --tr1        Fisrt mate of paired-end trascripts reads, can be specified more than once. 
 --tr2        Second mate of paired-end trascripts reads, can be specified more than once.
 
=head2 Global setting

 -d           Output directory.
 -o           Output prefix.
 -c           Config file, specify the path of dependencies, data file can also be specified.
 -t           Max CPU number to use, default: 4.

=head2 Busco setting

 -l           Busco lineage dataset, odb9(https://busco-archive.ezlab.org/v3/) is supported, busco will not run without lineage dataset.
 -s           Species used for augustus training in busco.

=head1 Note: 

 Gaap will find a dependency in enviromnent if its path has not been specified 
 in command line or config file. If a dependency can't be found, the corresponding 
 assessments will running with error. It won't affect the assessments without 
 this dependency.

=cut

use strict;
use Getopt::Long;
use FindBin qw($RealBin);
#use Data::Dumper;

my $task = 'gaap pipe';
my $mode = 0;

my %data = (
	assembly   => '',
	TGS        => {},
	NGS        => {},
	TRANS      => {},
);
my ($busco_lineage, $species, $conf,  $dir, $prefix, $samtools, $help);
my $threads = 4;
GetOptions(
	"r:s"             => \$data{assembly},
	"lr:s"            => \@{$data{TGS}->{data}},
	"x|type:s"        => \$data{TGS}->{type},
	"sr1:s"           => \@{$data{NGS}->{data1}},
	"sr2:s"           => \@{$data{NGS}->{data2}},
	"tr1:s"           => \@{$data{TRANS}->{data1}},
	"tr2:s"           => \@{$data{TRANS}->{data2}},
	
	"l:s"             => \$busco_lineage,
	"s:s"             => \$species,
	"c:s"             => \$conf,
 	"t:i"             => \$threads,
	"d:s"             => \$dir,
	"o:s"             => \$prefix,
	"samtools"        => \$samtools,
	"h"               => \$help
);
die `pod2text $0` if $help;
die "[$task] Error! No assembly file has input.\n", `pod2text $0` unless $data{assembly};

#print Dumper(\%data);

$dir = "./gaap_output_$$" unless $dir;
$dir =~ s/\/$//;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix = "gaap_output" unless $prefix;
my $prefix_out = "$dir/$prefix" if $dir;

my %software = ();
my %result = ();
my @cmd = ();
my $checkpoint = 0x0000;

if ($conf) {
	parse_conf($conf, \%data, \%software);
}

#print Dumper(\%data);
#print Dumper(\%software);
##path
my $geno_stat     = "${prefix_out}_genome_stat.txt";

my $mapping_tgs   = "$dir/mapping/tgs_mapping";
my $bam_tgs       = "$mapping_tgs/${prefix}_TGS_mapping.bam";
my $maprate_tgs   = "$dir/mapping/tgs_mapping/flagstat_tgs.txt";
my $bkp_dir       = "$dir/breakpoint_output/";

my $mapping_ngs   = "$dir/mapping/ngs_mapping";
my $bam_ngs       = "$mapping_ngs/$prefix.paired.sorted.bam";
my $maprate_ngs   = "$dir/mapping/ngs_mapping/flagstat_ngs.txt";
my $var_dir       = "$dir/variants_calling/";
my $vcf           = "$var_dir/$prefix.flt.vcf";

my $mapping_trans = "$dir/mapping/trans_mapping";
my $bam_trans     = "$mapping_trans/${prefix_out}_TRANS_mapping.bam";
my $maprate_trans = "$dir/mapping/trans_mapping/flagstat_trans.txt";

my $snvcov_tgs    = "$dir/snv_coverage_tgs";
my $snvcov_ngs    = "$dir/snv_coverage_ngs";

my $busco_dir     = "$dir/busco_out"; 

my $bkp_out       = "$dir/breakpoint_output";
my $accu_out      = "$dir/accuracy_output";
my $snvcov_out    = "$dir/snvcov_output";

print "\n";
print join("\n", $mapping_ngs, $mapping_tgs, $mapping_trans, $geno_stat , $bkp_out , $accu_out , $snvcov_out);
print "\n";

##cmd list
#genome stat
push @cmd, [0x0000, 0x0001, "[$task] Now getting genome assembly informations.\n",
            "$RealBin/basic.pl $data{assembly} > ${prefix_out}_genome_stat.txt"];
#TGS mapping
if (@{$data{TGS}->{data}}) {
	my $tgs_cmd = "$RealBin/run_minimap2.pl -r $data{assembly} ";
	$tgs_cmd .= "-i $_ " for @{$data{TGS}->{data}};
	$tgs_cmd .= "-x $data{TGS}->{type} " if $data{TGS}->{type};
	$tgs_cmd .= "-t $threads "   if $threads;
	$tgs_cmd .= "-d $mapping_tgs ";
	$tgs_cmd .= "-o $prefix ";
	$tgs_cmd .= "--minimap2 $software{minimap2}->{minimap2} "if $software{minimap2}->{minimap2};
    $tgs_cmd .= "--samtools $software{samtools}->{samtools} "if $software{samtools}->{samtools};
	push @cmd, [0x0000, 0x8000, "[$task] Now mapping TGS reads.\n", $tgs_cmd];
	
	#mapping rate
	my $flg_cmd = "$RealBin/mapping_rate.pl -b $bam_tgs ";
	$flg_cmd .= "-t $threads "   if $threads;
	$flg_cmd .= "-o $maprate_tgs ";
	$flg_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x8000, 0x0002, "[$task] Now running samtools flagstat to gain mapping rate of tgs data.\n", $flg_cmd];
	
	#breakpoint detect
	my $bkp_cmd = "$RealBin/breakpoint_detected.pl -b $bam_tgs ";
	$bkp_cmd .= "-t $threads " if $threads;
	$bkp_cmd .= "-d $bkp_dir ";
	$bkp_cmd .= "-o $prefix ";
	$bkp_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x8000, 0x0004, "[$task] Now detecting breakpoints.\n", $bkp_cmd];
}

#NGS mapping
if (@{$data{NGS}->{data1}} or @{$data{NGS}->{data2}}) {
	my $ngs_cmd = "$RealBin/run_bwa.pl -r $data{assembly} ";
	if (@{$data{NGS}->{data1}}) {
		$ngs_cmd .= "-i $_ " for @{$data{NGS}->{data1}};
	}
	if (@{$data{NGS}->{data2}}) {
		$ngs_cmd .= "-I $_ " for @{$data{NGS}->{data2}};
	}
	$ngs_cmd .= "-t $threads "   if $threads;
	$ngs_cmd .= "-d $mapping_ngs ";
	$ngs_cmd .= "-o $prefix ";
	$ngs_cmd .= "--bwa $software{bwa}->{bwa} " if $software{bwa}->{bwa};
	$ngs_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x0000, 0x4000, "[$task] Now mapping NGS reads.\n", $ngs_cmd];
	
	#mapping rate
	my $flg_cmd = "$RealBin/mapping_rate.pl -b $bam_ngs ";
	$flg_cmd .= "-t $threads "   if $threads;
	$flg_cmd .= "-o $maprate_ngs ";
	$flg_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x4000, 0x0008, "[$task] Now running samtools flagstat to gain mapping rate of tgs data.\n", $flg_cmd];
	
	my $var_cmd = "$RealBin/var_calling.pl -r $data{assembly} ";
	$var_cmd .= "-b $bam_ngs ";
	$var_cmd .= "-t $threads "   if $threads;
	$var_cmd .= "-d $var_dir ";
	$var_cmd .= "-o $prefix ";
	$var_cmd .= "--bcftools $software{bcftools}->{bcftools} " if $software{bcftools}->{bcftools};
	$var_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x4000, 0x1000, "[$task] Now calling variants.\n", $var_cmd];
}

#snv_coverage dot-plot
if (-e $bam_tgs) {
	my $sc_tgs_cmd = "$RealBin/cov_snp_dot.pl -r $data{assembly} -b $bam_tgs -v $vcf ";
	$sc_tgs_cmd   .= "-d $snvcov_tgs ";
	$sc_tgs_cmd   .= "-o $prefix ";
	$sc_tgs_cmd   .= "--bedtools $software{bedtools}->{bedtools} " if $software{bedtools}->{bedtools};
	$sc_tgs_cmd   .= "--bcftools $software{bcftools}->{bcftools} " if $software{bcftools}->{bcftools};
	$sc_tgs_cmd   .= "--Rscript $software{Rscript}->{Rscript} " if $software{Rscript}->{Rscript};
	push @cmd, [0x9000, 0x0010, "[$task] Now plot snv_coverage dot-plot using tgs bam file.\n", $sc_tgs_cmd];
}elsif (-e $bam_ngs) {
	my $sc_ngs_cmd = "$RealBin/cov_snp_dot.pl -r $data{assembly} -b $bam_ngs -v $vcf ";
	$sc_ngs_cmd   .= "-d $snvcov_ngs ";
	$sc_ngs_cmd   .= "-o $prefix ";
	$sc_ngs_cmd   .= "--bedtools $software{bedtools}->{bedtools} " if $software{bedtools}->{bedtools};
	$sc_ngs_cmd   .= "--bcftools $software{bcftools}->{bcftools} " if $software{bcftools}->{bcftools};
	$sc_ngs_cmd   .= "--Rscript $software{Rscript}->{Rscript} " if $software{Rscript}->{Rscript};
	push @cmd, [0x5000, 0x0020, "[$task] Now plot snv_coverage dot-plot using ngs bam file.\n", $sc_ngs_cmd];
}




#Trans mapping
if (@{$data{TRANS}->{data1}} or @{$data{TRANS}->{data2}}) {
	my $tr_cmd = "$RealBin/run_hisat2.pl -r $data{assembly} ";
	if (@{$data{TRANS}->{data1}}) {
		$tr_cmd .= "-i $_ " for @{$data{TRANS}->{data1}};
	}
	if (@{$data{TRANS}->{data2}}) {
		$tr_cmd .= "-I $_ " for @{$data{TRANS}->{data2}};
	}
	$tr_cmd .= "-t $threads "   if $threads;
	$tr_cmd .= "-d $mapping_trans ";
	$tr_cmd .= "-o $prefix ";
	$tr_cmd .= "-o $prefix ";
	$tr_cmd .= "--hisat2 $software{hisat2}->{hisat2} " if $software{hisat2}->{hisat2};
	$tr_cmd .= "--hisat2_build $software{hisat2_build}->{hisat2_build} " if $software{hisat2_build}->{hisat2_build};
	$tr_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x0000, 0x2000, "[$task] Now mapping transcripts reads.\n", $tr_cmd];
	
	#mapping rate
	my $flg_cmd = "$RealBin/mapping_rate.pl -b $bam_trans ";
	$flg_cmd .= "-t $threads "   if $threads;
	$flg_cmd .= "-o $maprate_trans ";
	$flg_cmd .= "--samtools $software{samtools}->{samtools} " if $software{samtools}->{samtools};
	push @cmd, [0x2000, 0x0040, "[$task] Now running samtools flagstat to gain mapping rate of tgs data.\n", $flg_cmd];
}

#run busco
if ($busco_lineage) {
	my $busco_cmd = "$RealBin/run_busco.pl -r $data{assembly} ";
	$busco_cmd .= "-l $busco_lineage ";
	$busco_cmd .= "-t $threads "   if $threads;
	$busco_cmd .= "-d $busco_dir ";
	$busco_cmd .= "-o $prefix ";
	$busco_cmd .= "-s $species " if $species;
	$busco_cmd .= "--tblastn   $software{tblastn}->{tblastn}     " if $software{tblastn}->{tblastn};    
	$busco_cmd .= "--augustus  $software{augustus}->{augustus}   " if $software{augustus}->{augustus};
	$busco_cmd .= "--hmmsearch $software{hmmsearch}->{hmmsearch} " if $software{hmmsearch}->{hmmsearch};
	push @cmd, [0x0000, 0x0080, "[$task] Now running busco.\n", $busco_cmd];
}
	
##running
if ($mode == 1) {
	for (@cmd) {
		if (($checkpoint & $_->[0]) == $_->[0]) {
			print $_->[2];
			print $_->[3],"\n";
			print "[$task] Running completed!\n\n";
			$checkpoint |= $_->[1];
		}
	}
}else {
	for (@cmd) {
		if (($checkpoint & $_->[0]) == $_->[0]) {
			print $_->[2];
			if (system $_->[3]) {
				print "[$task] Running with error!\n\n";
			}else {
				print "[$task] Running completed!\n\n";
				$checkpoint |= $_->[1];
			}
		}
	}
}

##---------------------------------subroutine---------------------------------##

sub parse_conf {
	my $conf = shift;
	my $data = shift;
	my $sofw = shift;
	
	my $phrase = '';
	my $temp = {};
	open IN, '<', $conf or die "[$task] Can't open such file: $conf.\n";
	while (<IN>) {
		next if /^#/;
		next if /^\s*$/;
		if (/^\[(\S+)\]/) {
			$phrase = $1;
			$temp = $data;
		}elsif (/^\{(\S+)\}/) {
			$phrase = $1;
			$temp = $sofw;
		}else {
			chomp;
			my ($opt, $value) = split '=', $_, 2;
			next if !$opt or !$value;
			$opt =~ s/^\s+//;$opt =~ s/\s+$//;
			$value =~ s/^\s+//;$value =~ s/\s+$//;
			if ($opt =~ /^data/) {
				push @{$temp->{$phrase}->{$opt}}, $value;
			}else {
				$temp->{$phrase}->{$opt} = $value;
			}
			
		}
	}
	close IN;
	return;
}

sub check_dependencies {
	my @dep = qw ( minimap2 samtools bwa bcftools hisat2 Rscript 
	               busco tblastn augustus hmmsearch );
	for my $software (@dep) {
		$software = $software ? check_software("$software", $software) : check_software("software");
		die "[$task] Error! $software is not found.
		[$task] Check your config file or set it in your environment.\n" if $software eq "-1";
	}
}

sub check_software {
	my $software = shift;
	my $path = shift if @_;
	if ( $path ) {
		if ( basename($path) eq "software" && -X $path ) {
			$software = $path;
		}else {
			return "1";
		}
	}else {
		my @path = split /:/, $ENV{PATH};
		foreach ( @path ) {
			$_ =~ s/\/$//; 
			return $software = "$_/$software" if -X "$_/$software";
		}
		return "1";
	}
}




__END__

##DATA (delete the # or add a new line if you have more than one sample)
[TGS]
data=
#data=
type=

[NGS]
data1=
data2=
#data1=
#data2=

[TRANS]
data1=
data2=
#data1=
#data2=
 
##Dependencies(path to executable file)
{samtools}
samtools=

{minimap2}
minimap2=

{bwa}
bwa=

{bcftools}
bcftools=

{hisat2}
hisat2=
{hisat2-build}
hisat2-build=

{Rscript}
Rscript=

#busco depencies
{tblastn}
tblastn=
{augustus}
augustus=
{hmmsearch}
hmmsearch=