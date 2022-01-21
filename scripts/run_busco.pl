#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use FindBin qw($RealBin);


my $task = "run_busco";
my $mode = 0;

my $usage = <<USAGE;
Program: $task
Version: v1.0

Command: 
    $task -r assembly [option]

Options: 
    -d <PATH>       output directory
    -o <STR>        output prefix
    -t <INT>        number of threads to use
    -i <path>       busco config file
    -l <PATh>       busco lineage path
    -s <STR>        augustus species
    --tblastn       path to tblastn binary
    --augustus      path to augustus binary
    --hmmsearch     path to hmmsearch binary
    -h              print help massage
USAGE

my $config = {};
my ($assembly,$dir,$prefix_out,$threads,$fconfig,$lineage,$species,$help);
GetOptions(
	"r:s"            => \$assembly,
	"d:s"            => \$dir,
	"o:s"            => \$prefix_out,
	"t:i"            => \$threads,
	"i:s"            => \$fconfig,
	"l:s"            => \$lineage,
	"s:s"            => \$species,
	"tblastn:s"      => \$config->{tblastn},
	"augustus:s"     => \$config->{augustus},
	"hmmsearch:s"    => \$config->{hmmsearch},
	"h"              => \$help,
);

die $usage if $help || (!$assembly && !$lineage);

if ($dir && !-e $dir){
	if (system "mkdir $dir"){
		die "[$task]Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix_out   = "output" unless $prefix_out;

$ENV{PATH} = "$ENV{PATH}".":".dirname($config->{'augustus'});
my $AUGUSTUS_CONFIG_PATH = dirname($config->{'augustus'})."/../config";
$ENV{AUGUSTUS_CONFIG_PATH} = $AUGUSTUS_CONFIG_PATH;
if (!$fconfig) {
	my $default_config = "$RealBin/../third_party/busco/config/config.ini.default";
	$fconfig = create_config($config, $default_config);
}
my $BUSCO_CONFIG_FILE    = "$fconfig";
$ENV{BUSCO_CONFIG_FILE}    = $BUSCO_CONFIG_FILE;


my $busco_cmd = "$RealBin/../third_party/busco/build/run_BUSCO.py -m genome ";
$busco_cmd   .= "-i $assembly ";
$busco_cmd   .= "-l $lineage ";
$busco_cmd   .= "-o $prefix_out ";
$busco_cmd   .= "-c $threads " if $threads;

_system($busco_cmd, $mode);



##-------------------------------subroutines----------------------------------##

sub create_config {
	my $config  = shift;
	my $def_con = shift;
	$config->{'tblastn'}              = dirname($config->{'tblastn'});
	$config->{'makeblastdb'}          = $config->{'tblastn'};
	$config->{'augustus'}             = dirname($config->{'augustus'});
	$config->{'etraining'}            = $config->{'augustus'};
	$config->{'gff2gbSmallDNA.pl'}    = $config->{'augustus'};
	$config->{'new_species.pl'}       = $config->{'augustus'};
	$config->{'optimize_augustus.pl'} = $config->{'augustus'};
	$config->{'hmmsearch'}            = dirname($config->{'hmmsearch'});
	
	for (keys %{$config}) {
		my $check = $config->{$_} ? check_software($_, "$config->{$_}/$_") : check_software($_);
		die "[$task] Error! Can't find $_.\nPlease input a correct path or add it to environment.\n" if $check eq "-1";
	}
	
	print "$def_con\n";
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