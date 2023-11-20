#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);


my $task = "run_busco5";
my $mode = 0;

my $usage = <<USAGE;
Program: gaep busco

Command: 
    gaep busco <-r assembly> [options]

Options: 
    -d <PATH>       output directory
    -o <STR>        output prefix
    -t <INT>        number of threads to use
    -i <path>       busco config file
    -l <PATh>       busco lineage path
    -h              print help massage
    --busco         path to busco v5
    --prodigal      path to prodigal binary, for non-eukaryotes
    --metaeuk       path to metaeuk binary, for eukaryotes
    --hmmsearch     path to hmmsearch binary
    --bbtools       path to bbtools stats.sh

Note:
    GAEP calls BUSCO in offline mode, so please provide an unpacked lineage dataset, 
	which can be downloaded from https://busco-data.ezlab.org/v5/data/lineages/.
    If the paths of the dependencies are not specified in the command line, GAEP 
    will search for them in the environment.
USAGE

my $config = {};
my ($assembly,$dir,$prefix_out,$threads,$fconfig,$lineage,$busco,$help);#`$species,
my $domain = '';
GetOptions(
	"r:s"            => \$assembly,
	"d:s"            => \$dir,
	"o:s"            => \$prefix_out,
	"t:i"            => \$threads,
	"i:s"            => \$fconfig,
	"l:s"            => \$lineage,
	"busco:s"        => \$busco,
	"prodigal:s"     => \$config->{prodigal},
	"metaeuk:s"      => \$config->{metaeuk},
	"hmmsearch:s"    => \$config->{hmmsearch},
	"bbtools:s"      => \$config->{bbtools},
	"h"              => \$help,
);

die $usage if $help || (!$assembly || !$lineage);

$lineage = "./$lineage" if $lineage !~ /\//;
if (-e "$lineage/dataset.cfg") {
	open CFG, '<', "$lineage/dataset.cfg" or die "[$task]Error! Can't find $lineage/dataset.cfg, make sure the lineage dataset is complete.\n";
	while (<CFG>) {
		if ($_ =~ /^domain=(\S+)/) {
			$domain = $1;
			last;
		}
	}
	close CFG;
}else {
	die "[$task]Error! Can't find $lineage/dataset.cfg, make sure the lineage dataset is complete.\n";
}

$busco = $busco ? check_software("busco", $busco) : check_software("busco");
die "[$task] Error! Busco is not found. 
[$task] Check your config file or set it in your environment.\n" if $busco eq "-1";

if ($dir && !-e $dir){
	if (system "mkdir $dir"){
		die "[$task]Error! Can't make directory:\"$dir\"\n";
	}
}
#elsif(!$dir) {
#	$dir='./run_busco5_output';
#}
$prefix_out = "busco_output" unless $prefix_out;

#$ENV{PATH} = "$ENV{PATH}".":".dirname($config->{'augustus'});
#my $AUGUSTUS_CONFIG_PATH = dirname($config->{'augustus'})."/../config";
#$ENV{AUGUSTUS_CONFIG_PATH} = $AUGUSTUS_CONFIG_PATH;
if (!$fconfig) {
	my $default_config = "$RealBin/../third_party/busco/busco5_config/config.ini";
	$fconfig = create_config($config, $default_config, $domain);
}
my $BUSCO_CONFIG_FILE    = "$fconfig";
$ENV{BUSCO_CONFIG_FILE}    = $BUSCO_CONFIG_FILE;


my $busco_cmd = "$busco -m genome --offline -f ";
$busco_cmd   .= "--config $fconfig ";
$busco_cmd   .= "-i $assembly ";
$busco_cmd   .= "-l $lineage ";
#$busco_cmd   .= "-sp $species " if $species;;
$busco_cmd   .= "-o $prefix_out ";
$busco_cmd   .= "-c $threads " if $threads;

_system($busco_cmd, $mode);



##-------------------------------subroutines----------------------------------##

sub create_config {
	my $config  = shift;
	my $def_con = shift;
	if (lc($domain) eq 'eukaryota') {
		$config->{metaeuk} = $config->{metaeuk} ? check_software("metaeuk", $config->{metaeuk}) : check_software("metaeuk");
		die "[$task] Error! Can't find metaeuk.\nPlease input a correct path or add it to environment.\n" if $config->{'metaeuk'} eq "-1"; 
		$config->{metaeuk}   = dirname($config->{metaeuk});
	}else {
		$config->{prodigal} = $config->{prodigal} ? check_software("prodigal", $config->{prodigal}) : check_software("prodigal");
		die "[$task] Error! Can't find prodigal.\nPlease input a correct path or add it to environment.\n" if $config->{'prodigal'} eq "-1"; 
		$config->{prodigal}  = dirname($config->{prodigal});
	}

	$config->{bbtools} = $config->{bbtools} ? check_software("stats.sh", $config->{bbtools}) : check_software("stats.sh");
	die "[$task] Error! Can't find bbtools.\nPlease input a correct path to \"stats.sh\" or add it to environment.\n" if $config->{'bbtools'} eq "-1"; 
	$config->{bbtools} = dirname($config->{bbtools});

	$config->{hmmsearch} = $config->{hmmsearch} ? check_software("hmmsearch", $config->{hmmsearch}) : check_software("hmmsearch");
	die "[$task] Error! Can't find hmmsearch.\nPlease input a correct path or add it to environment.\n" if $config->{'hmmsearch'} eq "-1"; 
	$config->{hmmsearch} = dirname($config->{hmmsearch});
	
	open DEFT, '<', $def_con or die "Can't open such file: $def_con.\n";
	open CONFIG, '>', "$dir/busco_config.ini" or die "Can't open such file: $dir/busco_config.ini.\n";
	my $mask = '';
	while (<DEFT>) {
		if (/^\[(\S+)\]$/) {
			if (exists $config->{$1}) {
				$mask = $config->{$1};
			}
			print CONFIG $_;
		}elsif (/^path/) {
			print CONFIG "path = $mask\n";
			$mask = ''
		}elsif (/^;?out_path/) {
			print CONFIG "out_path = $dir\n";
		}else {
			print CONFIG $_;
		}
	}
	close DEFT;
	close CONFIG;
	return "$dir/busco_config.ini";
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