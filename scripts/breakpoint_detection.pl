#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Bio::DB::HTS;
use threads;
use Thread::Semaphore;
use threads::shared;
use Thread::Queue;
use FindBin qw($RealBin);

my $task = "bkp_detection";
my $mode = 0;

my $usage = <<USAGE;
Program: gaep bkp

Usage:   gaep bkp <-r fasta_file> [options]
Options:
         -r <FILE>           genome assembly in a fasta file
         -i <FILE>           Fastq file. Can be set more than one time.
         -l <FILE>           File listing fastq files. One per line.
         -x <STR>            Reads type. (ont,pb,ccs)
         -b <FILE>           bam alignment by TGS data, will ignore -i, -l and -x
         -e <INT>            give a estimated depth
         -c <INT>            a read with more than INT bp cliped as cliped read, default 1000
         -g <INT>            ignore i bp at both end of contigs
         -d <PATH>           output directory
         -o <STR>            output prefix
         -t <INT>            parallel number
         -h                  print help massage
         --samtools <PATH>   path to samtools
USAGE

my $threads;
my ($assembly,$file_list,$read_type,$bam,$bam_depth,$reads_type,$dir,$prefix,$no_index,$help,$samtools);
my $site_distance  = 500;
my $ignore_end     = 5000;
my $clip_threshold = 1000;
my @sequences;
GetOptions(
	"r:s"      => \$assembly,
	"i:s"      => \@sequences,      #FASTQ or FASTA
	"l|list:s" => \$file_list,
	"x:s"      => \$read_type,	    #(ont,pb)
	"b:s"      => \$bam,		   	#bam alignment by TGS data
	"e:i"      => \$bam_depth,     	#give a estimated depth 
	"c:i"      => \$clip_threshold,	
	"g:i"      => \$ignore_end,    	#ignore i bp at both end of contigslength regarded as a assumptive misjoin
 	"t:i"      => \$threads,
	"d:s"      => \$dir,
	"o:s"      => \$prefix,
	"samtools:s" => \$samtools,
	"h"        => \$help
);

##Initial setting
die $usage if $help;
if (!$bam && !$assembly) {
	print "[$task] Make sure to give a bam file or give the assembly and TGS reads for alignment.\n";
	die $usage;
	
}elsif (!$bam && $assembly && !(@sequences || $file_list)) {
	print "[$task] Error! No reads file.\n";
	die $usage;
}
unless ($threads) {
	print "[$task] No thread given, defult as 4.\n";
	$threads = 4;
}
##Check dependencies.
$samtools = $samtools ? check_software("samtools", $samtools) : check_software("samtools");
die "[$task] Error! Samtools is not found.
[$task] Check your config file or set it in your environment.\n" if $samtools eq "-1";


##Check input files.
$dir = "gaep_${task}_$$" unless $dir;
if (! -e $dir){
	if (system "mkdir -p $dir"){
		die "[$task] Error! Can't make directory:\"$dir\"\n";
	}
}
$prefix   = "bkp_output_$$" unless $prefix;
my $prefix_out   = "$dir/$prefix" if $dir;
print "\n[$task] Output directory is $dir.\n";

## reads mapping
if (!$bam) {
	my $map_cmd = "$RealBin/run_minimap2.pl -r $assembly ";
	if (@sequences) {
		$map_cmd .= "-i $_ " for @sequences;
	}
	$map_cmd .= "-l $file_list " if $file_list;
	$map_cmd .= "-x $read_type " if $read_type;
	$map_cmd .= "-t $threads "   if $threads;
	$map_cmd .= "-d $dir/mapping/TGS_mapping ";
	$map_cmd .= "-o $prefix ";
	_system($map_cmd, $mode);
	$bam = "$dir/mapping/TGS_mapping/${prefix}_TGS_mapping.bam"
}


#---------------------------------variates-------------------------------------#

my $temp_dir = "$dir/temp";
if (! -d $temp_dir) {
	mkdir "$dir/temp" or die "[$task] Error! Can't make directory:\"$dir\"\n";
}
my $posi_prefix      = "${temp_dir}/${prefix}_temp_clip_site";
my $sort_prefix      = "${temp_dir}/${prefix}_temp_clip_site_sorted";
my $ab_prefix        = "${temp_dir}/${prefix}_temp_ab_region";
my $bkp_t_prefix     = "${temp_dir}/${prefix}_temp_breakpoints";

my $abr_prefix       = "${prefix_out}_ab_region";
my $result_prefix    = "${prefix_out}_breakpoints";
#if (-s $result){
#	unlink $result || die "[$task]Error! File already exists: $result.\n";
#}
my %seq_len   = ();
my @ctg_list  = ();
my $tolen     = 0;
my $win_queue = Thread::Queue->new();
my $bkp_queue = Thread::Queue->new();
my $flt_queue = Thread::Queue->new();
my $semaphore = Thread::Semaphore->new($threads);
#-----------------------------------main---------------------------------------#

my $pip = "$samtools view -H $bam ";
open my $BAM, '-|', $pip || die "[$task 1] Fail to read header of $bam.\n";
while (<$BAM>) {
	if (/^\@SQ\s+SN:(\S+)\s+LN:(\d+)/) {
		$tolen += $2;
		$seq_len{$1} = $2;
		$bkp_queue->enqueue($1);
		$flt_queue->enqueue($1);
		push @ctg_list, $1;
	}
}
close $BAM;
#if ($gzip ne "-1") {
#	read_header_gzip($bam, $tolen, \%seq_len, $bkp_queue, $flt_queue, \@ctg_list);
#}else {
#	read_header_samtools($bam, $tolen, \%seq_len, $bkp_queue, $flt_queue, \@ctg_list);
#}
$bkp_queue->end();
$flt_queue->end();
die "[$task]Fail to read header of $bam.\n" unless (%seq_len);

#If the bam has indexed?
unless (-s "${bam}.bai"){
	my $cmd = "$samtools index ";
	$cmd   .= "-@ $threads " if $threads;
	$cmd   .= "$bam";
	_system($cmd);
}

##calculate mean coverage and standard deviation
my @tocov :shared;
@tocov = ();

for my $ctg (keys %seq_len) {
	my $window_len = 0;
	while ($window_len < $seq_len{$ctg} * 0.1) {
		$win_queue->enqueue([ran_window($ctg, $seq_len{$ctg})]);
		$window_len += 1000;
	}
}
$win_queue->end();

for (my $i = 0; $i < $threads; $i++) {
	$semaphore->down(1);
	threads->create(\&cov_process, $win_queue) -> detach();
}

$semaphore->down($threads);
$semaphore->up($threads);

my ($mean, $sd) = cal_cov_stat(\@tocov);
print STDERR "Mean coverage: $mean\nStandard deviation: $sd\n";
$bam_depth = $mean if !$bam_depth;
#print $mean,"\t", $sd, "\n";

##detect breakpoints
print STDERR "[$task] Parsing bam to get clip site.\n";

for (my $i = 0; $i < $threads; $i++) {
	$semaphore->down(1);
	threads->create(\&bkp_process,$bkp_queue) -> detach();
}

$semaphore->down($threads);
$semaphore->up($threads);

##breakpoints filter
print STDERR "[$task] Filtering breakpoint regions.\n";

for (my $i = 0; $i < $threads; $i++) {
	$semaphore->down(1);
	threads->create(\&br_filter_process, $bam, $flt_queue, $mean, $sd) -> detach();
}

$semaphore->down($threads);
$semaphore->up($threads);



##Merge output
print "[$task] Merge output of each contig into ${result_prefix}.txt.\n";
open  ABR, '>', "${abr_prefix}.txt"    || die "Can't open such file: ${abr_prefix}.txt";
open  OUT, '>', "${result_prefix}.txt" || die "Can't open such file: ${result_prefix}.txt";
print OUT "#GAEP misassembly breakpoint detection.\n";
print OUT "#INFO: span: the count of reads that span the breakpoint region.\n";
print OUT "#INFO: no_span: the count of reads that do not span the breakpoint region.\n";
print OUT "#INFO: Lcov: the mean coverage of a 1k region to the left of the breakpoint region.\n";
print OUT "#INFO: Bcov: mean coverage of the breakpoint region, -1 if the region length less than 1K.\n";
print OUT "#INFO: Rcov: the mean coverage of a 1k region to the right of the breakpoint region.\n";
print OUT "#INFO: type: breakpoint type. 1: No reads span the breakpoint; 2: Only a low proportion of reads supports the contig; 3: Breakpoint caused by repeats.\n";
print OUT join("\t", "#contig", "start", "end", "span", "no_span", "Lcov", "Bcov", "Rcov", "type"),"\n";
for my $ctg (@ctg_list) {
	if (-s "${ab_prefix}_$ctg.txt") {
		open RES, '<', "${ab_prefix}_$ctg.txt" || die "Can't open such file: ${ab_prefix}_$ctg.txt";
		while (<RES>) {
			print ABR $_;
		}
		close RES;
	}
	if (-s "${bkp_t_prefix}_$ctg.txt") {
		open RES, '<', "${bkp_t_prefix}_$ctg.txt" || die "Can't open such file: ${bkp_t_prefix}_$ctg.txt";
		while (<RES>) {
			my $type = (split)[8];
			print OUT if $type > 0;
		}
		close RES;
	}
}
close ABR;
close OUT;

##-------------------------------subroutines----------------------------------##
                                          
sub read_header_samtools {
	my $bam       = shift;
	my $tolen     = shift;
	my $seq_len   = shift;
	my $bkp_queue = shift;
	my $flt_queue = shift;
	my $ctg_list  = shift;
	my $pip = "$samtools view -H $bam ";
	open my $BAM, '-|', $pip || die "[$task 1]Fail to read header of $bam.\n";
	while (<$BAM>) {
		if (/^\@SQ\s+SN:(\S+)\s+LN:(\d+)/) {
			$seq_len->{$1} = $2;
			$tolen += $2;
			$bkp_queue->enqueue($1);
			$flt_queue->enqueue($1);
			push @$ctg_list, [$1,$2];			
		}
	}
	return;
}

##sub read_header_gzip {
##	my $bam       = shift;
##	my $tolen     = shift;
##	my $seq_len   = shift;
##	my $bkp_queue = shift;
##	my $flt_queue = shift;
##	my $ctg_list  = shift;
##	my $text = '';
##	my $n_targets = 0;
##	my $bam_is_be = 0;
##	if (byteorder() eq 'b') {
##		$bam_is_be = 1;
##	}elsif (byteorder() eq 'l') {
##		$bam_is_be = 0;
##	}else {
##		die "Can't decide of the byteorder.\n";
##	}
##	
##	open my $pipe, '-|', "$gzip -dc $bam" || die $!;
##	read($pipe, my $buf, 4);
##	if ($buf ne "BAM\001") {
##		die "$bam is not a BAM file.\n";
##	}
##	read($pipe, my $l_text, 4);
##	$bam_is_be ? read($pipe, $text, unpack("N", $l_text)):
##                 read($pipe, $text, unpack("L", $l_text));
##	
##	read($pipe, $n_targets, 4);
##	$n_targets = $bam_is_be ? unpack "N", $n_targets:
##                              unpack "L", $n_targets;
##	close $pipe;
##	
##	while ($text =~ /\@SQ\s+SN:(\S+)\s+LN:(\d+)/g) {
##		$seq_len->{$1} = $2;
##		$tolen += $2;
##		$bkp_queue->enqueue($1);
##		$flt_queue->enqueue($1);
##		push @$ctg_list, [$1,$2];
##	}
##	return;
##}
##
##sub byteorder {
##	my $test = "\x00\x01";
##	open IN, '<', \$test;
##	read IN, my $endian, 2;
##	close IN;
##	my $unpack = unpack "S", $endian;
##	if ($unpack == 0x0001) {
##		return 'b';
##	}elsif ($unpack == 0x0100){
##		return 'l';
##	}else {
##		return 'n';
##	}
##}

#sub length_check_fai{
#	my %seq_len = ();
#	unless (-s "${assembly}.fai"){
#		my $cmd = "$samtools faidx ";
#		$cmd   .= "$assembly";
#		send_to_system($cmd);
#	}
#	open my $FAI, '<', "${assembly}.fai" || die "Can't open such file: ${assembly}.fai";
#	while (<$FAI>){
#		my ($chr, $length) = (split)[0,1];
#		if (!exists $seq_len{$chr}){
#			$seq_len{$chr} = $length;
#		}else {
#			die "[$task]Error!More than one contig have a same name. Make sure your fasta file is valid.\n";
#		}
#	}
#	close $FAI;
#	return %seq_len;
#}

sub ran_window {
	my $ctg = shift;
	my $ctg_len = shift;
	my $ran_start = 0;
    do {
        $ran_start = int(rand ($ctg_len)) + 1; 
        #print join("\t", $ran_ctg, $ran_start, $ran_start + 999),"\n";
    } while ($ctg_len - $ran_start < 1000);
    return ($ctg, $ran_start, $ran_start + 999);
}

sub cov_process {
	my $SAM = Bio::DB::HTS->new(-bam => $bam);
	my $queue = shift;
	while (defined(my $region = $queue->dequeue())) {
		my ($ctg, $start, $end) = @$region;
		cal_cov($ctg, $start, $end, \$SAM);
	}
	$semaphore->up(1);
}

sub cal_cov {
	my ($ref, $start, $end, $SAM) = @_;
	$start ++;
	my $cov = 0;
	my @alignments = $$SAM->get_features_by_location(-seq_id => $ref,
                                                     -start  => $start,
                                                     -end    => $end);
	for my $align (@alignments) {
		next if  $align->unmapped;
		my $rstart = $align->pos+1;
		my $rend   = $align->calend;
		#print join("\t",$ref, $start, $end, $rstart, $rend),"\n";
		my $left  = $start >= $rstart ? $start : $rstart;
		my $right = $end   <= $rend   ? $end   : $rend;
		$cov += $right - $left + 1;
	}
	push @tocov, $cov/($end-$start+1);
	return;
}

sub cal_cov_stat {
	my $cov_list = shift;
	my $ave = 0;
	my $sd  = 0;
	@$cov_list = sort{$a<=>$b} @$cov_list;
	@$cov_list = @{$cov_list}[(@$cov_list*0.1)..(@$cov_list-@$cov_list*0.1)];
	grep {$ave += $_}@$cov_list;
	$ave = $ave/@$cov_list;
	grep {$sd += ($_-$ave)**2}@$cov_list;
	$sd = sqrt($sd/@$cov_list);
	return ($ave, $sd);
}

sub bkp_process {
	my $queue        = shift;
	my $BAM          = Bio::DB::HTSfile->open($bam);
	my $header       = $BAM->header_read;
	my $target_count = $header->n_targets;
	my $target_names = $header->target_name;
	my $index        = Bio::DB::HTSfile->index_load($BAM);
	
	while (defined(my $contig = $queue->dequeue())) {
		my ($tid, $start, $end) = $header->parse_region($contig);
		open my $POSI , '>', "${posi_prefix}_$contig.txt" or die "[$task]Can't open such file: ${posi_prefix}_$contig.txt";
		$index->fetch($BAM,$tid,$start,$end,\&parse_alignment,[$contig,$POSI]);
		close $POSI;
		
		if (-s "${posi_prefix}_$contig.txt") {
			my $cmd = "sort -k3,3n ${posi_prefix}_$contig.txt > ${sort_prefix}_$contig.txt";
			_system($cmd);
			merge_sites("${sort_prefix}_$contig.txt", "${ab_prefix}_$contig.txt") if (-s "${sort_prefix}_$contig.txt");
			unlink "${posi_prefix}_$contig.txt";
		}elsif (-z "${posi_prefix}_$contig.txt") {
			unlink "${posi_prefix}_$contig.txt";
		}
		
	}
	$semaphore->up(1);
	return;
}

sub parse_alignment {
	my $align = shift;
	my $data = shift;
	my ($contig, $POSI) = @{$data};
	
	my $flag = $align->flag;
	return if $flag&260;
	## ref
	my $ref    = $contig;
	my $rstart = $align->start;
	my $cigar  = $align->cigar_str;

	#parse cigar
	my %atype        = ();
	$atype{$_}       = 0 for ("S","H","M","I","D","N","P","=","X");
	while ($cigar =~ /(\d+)([SHMIDNP=X])/g){
		if ( $1 >= 1000 ) {
			if ( $2 eq 'D' ) {
				my $lbreakpoint = $rstart + $atype{"M"} + $atype{"X"} + $atype{"="} + $atype{"D"} - 1;
				my $rbreakpoint = $lbreakpoint + $1 + 1;
				print $POSI join("\t", $ref, "Dr", $lbreakpoint),"\n";
				print $POSI join("\t", $ref, "Dl", $rbreakpoint),"\n";
			}elsif ( $2 eq 'I' ) {
				my $lbreakpoint = $rstart + $atype{"M"} + $atype{"X"} + $atype{"="} + $atype{"D"} - 1;
				my $rbreakpoint = $lbreakpoint + 1;
				print $POSI join("\t", $ref, "Ir", $lbreakpoint),"\n";
				print $POSI join("\t", $ref, "Il", $rbreakpoint),"\n";
			}
		}
		$atype{"$2"} += $1;
	}
	##cal_end
	$atype{"M"}     += $atype{"X"} + $atype{"="};
	$atype{"D"}     += $atype{"N"};
	#my $qlength      = $atype{"S"}+$atype{"H"}+$atype{"M"}+$atype{"I"};
	my $align_length = $atype{"M"} + $atype{"D"};
	return if $align_length <= 500;
	my $rend         = $rstart + $align_length - 1;
	my $lclip_len = $cigar =~ /^(\d+)[SH]/ ? $1 : 0;
	my $rclip_len = $cigar =~ /(\d+)[SH]$/ ? $1 : 0;

	if ( $rstart >= $ignore_end && ( $seq_len{$ref} - $rend ) >= $ignore_end ) {
		if ( $lclip_len >= $clip_threshold && $rclip_len >= $clip_threshold ) {
			print $POSI join("\t", $ref, "dl", $rstart),"\n";
			print $POSI join("\t", $ref, "dr", $rend),"\n";
		}elsif ( $lclip_len >= $clip_threshold ) {
			print $POSI join("\t", $ref, "l", $rstart),"\n";
		}elsif ( $rclip_len >= $clip_threshold ) {
			print $POSI join("\t", $ref, "r", $rend),"\n";
		}
	}
	return;
}

sub merge_sites{
	
	my $posi = shift;
	my $out  = shift;
	
	open my $POSI,'<', $posi || die "Can't open such file: $posi.\n";
	open my $AB_COOR, '>', $out || die "Can't open such file: $out.\n";
	
	my $rseq  = "";
	my $tcoor = 0;
	my @clip;	    #number of clips in a region
	my %btype = ("Dl" => 0,"Dr" => 0,"Il" => 0,"Ir" => 0,"l"  => 0,"r"  => 0,"dl" => 0,"dr" => 0);

	while (<$POSI>){
		chomp;
		my @line  = split;
		next if ( $seq_len{$line[0]} <= 5000 );
		my $rcoor = $line[2];
		
		if ($line[0] ne $rseq){
			if (@clip && @clip + 1 >= $bam_depth/10){
				print $AB_COOR join("\t", $rseq, $clip[0], $tcoor, $tcoor-$clip[0]+1, @clip+1,  
				                    $btype{"Dl"}, $btype{"Dr"}, $btype{"Il"}, $btype{"Ir"}, 
									$btype{"l" }, $btype{"r" }, $btype{"dl"}, $btype{"dr"},
				                   ) ,"\n";
				#depth_check($rseq, $clip[0]-150, $tcoor+150);
				#break_reads_info($rseq, $clip[0]-150, $tcoor+150);
			}
			$rseq = $line[0];
			@clip = ();
			$tcoor = 0;
			%btype = ("Dl" => 0,"Dr" => 0,"Il" => 0,"Ir" => 0,"l"  => 0,"r"  => 0,"dl" => 0,"dr" => 0);
			$btype{$line[1]} += 1;
		}else{
			if ($rcoor - $tcoor <= $site_distance){
				push @clip, $tcoor;
				$tcoor = $rcoor;
				$btype{$line[1]} += 1;
			}elsif (@clip){
				if (@clip && @clip + 1 >= $bam_depth/10){
					print $AB_COOR join("\t", $rseq, $clip[0], $tcoor, $tcoor-$clip[0]+1, @clip+1,  
										$btype{"Dl"}, $btype{"Dr"}, $btype{"Il"}, $btype{"Ir"}, 
										$btype{"l" }, $btype{"r" }, $btype{"dl"}, $btype{"dr"},
										) ,"\n";
					#depth_check($rseq, $clip[0]-150, $tcoor+150);
					#break_reads_info($rseq, $clip[0]-150, $tcoor+150);
				}
				@clip = ();
				%btype = ("Dl" => 0,"Dr" => 0,"Il" => 0,"Ir" => 0,"l"  => 0,"r"  => 0,"dl" => 0,"dr" => 0);
				$btype{$line[1]} += 1;
			}
		}
		$tcoor = $rcoor;
	}
	
	if (@clip && @clip + 1 >= $bam_depth/10){
		print $AB_COOR join("\t", $rseq, $clip[0], $tcoor, $tcoor-$clip[0]+1, @clip+1,  
							$btype{"Dl"}, $btype{"Dr"}, $btype{"Il"}, $btype{"Ir"}, 
							$btype{"l" }, $btype{"r" }, $btype{"dl"}, $btype{"dr"},
							) ,"\n";
		#depth_check($rseq, $clip[0]-150, $tcoor+150);
		#break_reads_info($rseq, $clip[0]-150, $tcoor+150);
	}
	
	close $POSI;
	close $AB_COOR;
	return;
}

sub br_filter_process {
	my $bam      = shift;
	my $queue    = shift;
	my $mean_cov = shift;
	my $sd_cov   = shift;
	my $sam = Bio::DB::HTS->new(-bam => $bam);
	while (defined(my $contig = $queue->dequeue())) {
		if (-s "${ab_prefix}_$contig.txt") {
			open BKP, '<', "${ab_prefix}_$contig.txt" or die "[$task]Can't open such file: ${posi_prefix}_$contig.txt";
			open my $RES, '>', "${bkp_t_prefix}_$contig.txt" or die "[$task]Can't open such file: ${bkp_t_prefix}_$contig.txt";
			while (<BKP>) {
				chomp;
				####my ($ctg, $start, $end) = split;
				my ($ctg, $start, $end, @count) = split;
				next if ( ($count[8]+$count[9])/$count[1] >= 0.5 ); # $count[1] != 0
				my @alignments = $sam->get_features_by_location(-seq_id => $ctg,
																-start  => $start-1000,
																-end    => $end+1000);
				br_filter(\@alignments, $ctg, $start, $end, $mean_cov, $sd_cov, $RES);
				
			}	
			close BKP;
			close $RES;
		}
	}
	$semaphore->up(1);
}

sub br_filter {
	my $alignment = shift;
	my $ctg       = shift;
	my $start     = shift;
	my $end       = shift;
	my $mean_cov  = shift;
	my $sd_cov    = shift;
	my $FILE      = shift;
	
	my %read_info = (
		span    => 0,
		nospan  => 0,
	);
	my $Lcoverage = 0;
	my $Bcoverage = 0;
	my $Rcoverage = 0;
	my %cov = (
		Rcov => 0,
		Bcov => 0,
		Lcov => 0,
	);
	
	LINE: for my $align (@{$alignment}) {
		# ref
		my $ref    = $align->seq_id;
		my $rstart = $align->start;
		my $rend   = $align->end;
		my $cigar  = $align->cigar_str;
		# query
		next if ($rend - $rstart <= 500);
		
		## span or no span
		if ($rstart < $start && $rend > $end) {
			$read_info{span} += 1;
			#insertion or deletion
			my %atype  = ();
			$atype{$_}       = 0 for ("S","H","M","I","D","N","P","=","X");
			while ($cigar =~ /(\d+)([SHMIDNP=X])/g){
				if ( $1 >= 1000 ) {
					if ( $2 eq 'D' ) {
						my $lbreakpoint = $rstart + $atype{"M"} + $atype{"X"} + $atype{"="} + $atype{"D"} - 1 ;
						my $rbreakpoint = $lbreakpoint + $1 + 1;
						next if ($rbreakpoint < $start);
						last if ($lbreakpoint > $end);
						#next LINE if ($lbreakpoint < $start && $rbreakpoint > $end);
						$read_info{nospan} += 1;
						$read_info{span} -= 1;
						last;
					}elsif ( $2 eq 'I' ) {
						my $breakpoint  = $rstart + $atype{"M"} + $atype{"X"} + $atype{"="} + $atype{"D"} - 1;
						next if ($breakpoint + 1 < $start);
						last if ($breakpoint > $end);
						$read_info{nospan} += 1;
						$read_info{span} -= 1;
						last;
					}
				}
				$atype{"$2"} += $1;
			}
		}else {
			$read_info{nospan} += 1 unless ($rend < $start || $rstart > $end);
		}
		
		##clip length
		my $llclip = $cigar =~ /^(\d+)[SH]/ ? $1 : 0;
		my $lrclip = $cigar =~ /(\d+)[SH]$/ ? $1 : 0;
				
		##  coverage   ($rstart, $rend, $start, $end)
		##                     L       B        R
	    ##  =============[==========](===)[==========]===============
		##                    1K       br      1K
		if ( $rend < $start ) {
			$Lcoverage += $rstart <= $start - 1000 ? ( 1000 - ( $start - $rend ) ) :
			                                         ( $rend - $rstart );
		}elsif ( $rend >= $rstart && $rend <= $end ) {
			$Lcoverage += $rstart <= $start - 1000 ? 1000                 :
			              $rstart <  $start        ? ( $start - $rstart ) :
						                             0                    ;
			$Bcoverage += $rstart <  $start ? ( $rend - $start + 1 ) :
			                                  ( $rend - $rstart + 1 );
		}elsif ( $rend > $end && $rend <= $end + 1000 ) {
			$Lcoverage += $rstart <= $start - 1000 ? 1000                 :
			              $rstart <  $start        ? ( $start - $rstart ) :
						                             0                    ;
			
			$Bcoverage += $rstart <  $start ? ( $end - $start + 1 )  :
			              $rstart <= $end   ? ( $end - $rstart + 1 ) :
											  0                      ;
			
            $Rcoverage += $rstart <= $end   ? ( $rend - $end + 1 )    :
			                                  ( $rend - $rstart + 1 ) ;
		}else {
			$Lcoverage += $rstart <= $start - 1000 ? 1000                 :
			              $rstart <  $start        ? ( $start - $rstart ) :
						                             0                    ;
			
			$Bcoverage += $rstart <  $start ? ( $end - $start + 1 )  :
			              $rstart <= $end   ? ( $end - $rstart + 1 ) :
											  0                      ;
			
            $Rcoverage += $rstart <= $end   ? 1000                       : 
			              $rstart >= $end   ? ( $end + 1000 - $rstart )  :
											  0                          ;								  
		}
	}
	$read_info{Lcov} = $Lcoverage/1000;
	$read_info{Bcov} = $end-$start+1 < 1000 ? -1 : sprintf("%.2f",$Bcoverage/($end-$start+1));
	$read_info{Rcov} = $Rcoverage/1000;
	
	##filter
	my ($max_cov, $min_cov) = max_min_three($read_info{Lcov}, $read_info{Bcov}, $read_info{Rcov});
	my $total = $read_info{nospan}+$read_info{span};
	my $bkp_type = $read_info{span} == 0                                              ?  1 :
			       $read_info{nospan}/$total < 0.25                                   ? -1 :
			       $read_info{nospan}/$total > 0.75                                   ?  2 :
			       ($max_cov > $mean_cov+3*$sd_cov || $max_cov < $mean_cov-3*$sd_cov) ?  3 :
		                                                                                -2 ;
	print $FILE join("\t",$ctg,$start,$end,$read_info{span},$read_info{nospan},$read_info{Lcov},$read_info{Bcov},$read_info{Rcov},$bkp_type), "\n";
}

sub max_min_three {
	my $max = shift;
	my $min = $max;
	shift if $_[0] == -1;
	for my $n (@_) {
		$max = $n > $max ? $n : $max;
		$min = $n < $min ? $n : $min;
	}
	return ($max, $min);
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