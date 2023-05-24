#! /usr/bin/perl -w

use strict;
use 5.010;
my $fa = shift;
my $fai = "$fa.fai";
if (! -e $fai) {
	die '[Error] Can\'t find the fai file. You can run "samtools faidx" to generate a fai file.'
}

my @base = ('A','T','C','G');

my %hash;
my $id = '';
open FA, '<', $fa or die $!;
while (<FA>) {
	chomp;
	if (/^>(\S+)/) {
		$id = $1;
	}else {
		$hash{$id} .= $_;
	}
}
my %hash_simu = %hash;


my %fai = ();
my @ctgList = ();
open FAI, '<', $fai || die $!;
while ( <FAI> ) {
	next if /^\s*$/;
	chomp;
	my @line = split;
	$fai{$line[0]} = [@line[1,2,3,4]];
	push @ctgList, $line[0];
}
close $fai;
my @chooseList = @ctgList; 
@ctgList = grep { $fai{$_}->[0] >= 500000 }@ctgList;


my $ref_adjust = 0;
my $simu_adjust = 0;
##repeats expanded
my @posi_list = ();
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(1000,2000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_3"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 15000);
	my $insLen = rand_r(2000,5000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_3"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 30000);
	my $insLen = rand_r(5000,10000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_4"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 60000);
	my $insLen = rand_r(10000,20000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_4"];
}
for ( 1..25 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 150000);
	my $insLen = rand_r(20000,50000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_4"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 8000);
	my $insLen = rand_r(1000,2000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_3"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 20000);
	my $insLen = rand_r(2000,5000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_3"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 40000);
	my $insLen = rand_r(5000,10000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_4"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 80000);
	my $insLen = rand_r(10000,20000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_4"];
}
for ( 1..25 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 150000);
	my $insLen = rand_r(20000,50000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"expanded_4"];
}

##repeats collasped
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(1000,2000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"collasped"];
}
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 15000);
	my $insLen = rand_r(2000,5000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"collasped"];
}
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 30000);
	my $insLen = rand_r(5000,10000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"collasped"];
}
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 60000);
	my $insLen = rand_r(10000,20000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"collasped"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 150000);
	my $insLen = rand_r(20000,50000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"collasped"];
}

##inversion
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(1000,2000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"inv"];
}
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 15000);
	my $insLen = rand_r(2000,5000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"inv"];
}
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 30000);
	my $insLen = rand_r(5000,10000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"inv"];
}
for ( 1..100 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 60000);
	my $insLen = rand_r(10000,20000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"inv"];
}
for ( 1..50 ) {
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 150000);
	my $insLen = rand_r(20000,50000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"inv"];
}

##Insertion
for ( 1..100 ) {               #1000,2000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(1000,2000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"ins"];
}
for ( 1..100 ) {               #2000,5000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(2000,5000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"ins"];
}
for ( 1..100 ) {               #5000,10000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(5000,10000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"ins"];
}
for ( 1..100 ) {               #10000,20000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(10000,20000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"ins"];
}
for ( 1..100 ) {               #20000,50000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(20000,50000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"ins"];
}
##Deletion
for ( 1..100 ) {               #1000,2000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(1000,2000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"del"];
}
for ( 1..100 ) {               #2000,5000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(2000,5000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"del"];
}
for ( 1..100 ) {               #5000,10000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(5000,10000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"del"];
}
for ( 1..100 ) {               #10000,20000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(10000,20000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"del"];
}
for ( 1..100 ) {               #20000,50000
	my $ranCtg = $ctgList[rand(scalar @ctgList)];
	my $ran_position = rand_r(5000,$fai{$ranCtg}->[0] - 6000);
	my $insLen = rand_r(20000,50000);
	push @posi_list, [$ranCtg,$ran_position,$insLen,"del"];
}


for (@posi_list) {
	print join("\t",@{$_}),"\n";
}
#----------------------------------subroutine----------------------------------#

sub rand_r {
	my $start = shift;
	my $end = shift;
	return if $end <= $start;
	my $int = @_ ? shift : 1;
	return $int ? int (rand $end-$start+1) + $start : rand ($end-$start) + $start;
}

sub delete_seq {                    #delete_seq(seq, start position, seq length)
	my $seq = shift;
	my $start = shift;
	my $len = shift;
	my $left = substr($seq, 0, $start);
	my $right = substr($seq, $start+$len, (length $seq)-($start+$len));
	return $seq = join('', $left, $right);
}

sub insert_seq {
	my $seq = shift;
	my $target_seq = shift;
	my $start = shift;
	my $left = substr($target_seq, 0, $start);
	my $right = substr($target_seq, $start, (length $target_seq) - $start);
	return $target_seq = join('', $left, $seq , $right);
}

sub random_seq {
	my $len = @_; 
	my $random_seq = '';
	for ( 1..$len ) {
		$random_seq .= $base[int(rand(4))];
	}
	return $random_seq;
}
