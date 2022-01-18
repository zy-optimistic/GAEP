#! /usr/bin/perl -w

use strict;

my $seq_name;
my $total_length = 0;
my $sum_GC = 0;
my $GC_cotent = 0;
my $count = 0;
my $num_of_N = 0;
my @name;
my @length;
my %seq = ();

die "No such file: $ARGV[0].\n" unless (-s $ARGV[0]);
open IN ,'<', $ARGV[0] || die ("Can't open such file: $ARGV[0].\n");
##Read file and write sequences in hash.
print STDERR "Reading file.........\n";
while (<IN>){
	chomp;
	if (/^>(.*)/){
		$seq_name = $1;
		push @name , $seq_name;
		$count ++;
	}else {
		$length[$count-1] += length $_;
		$sum_GC += $_ =~ tr/gcGC/gcGC/;
		$num_of_N += $_ =~ tr/Nn/Nn/;
		$total_length += length $_;
	}
}
print STDERR "#############################\n";



@length = sort {$a <=> $b} @length;

$GC_cotent = $sum_GC / $total_length; ##GC(%)

#printf ( "Total length:\t%d\nN's:\t%d\nTotal length(without N):\t%d\nNumber of GC:\t%d\nGC%%:\t%.2f%%\nLongest:\t%d\n" ,$total_length,$num_of_N,$total_length-$num_of_N,$sum_GC,$GC_cotent * 100,$length[-1]);
printf ("%-25s%d\n","Total sequences:",$count);
printf ("%-25s%d\n" ,"Total length:",$total_length);
printf ("%-25s%d\n","N's:",$num_of_N);
printf ("%-25s%d\n" ,"Total length(without N):",$total_length-$num_of_N);
printf ("%-25s%d\n" ,"Number of GC:",$sum_GC );
printf ("%-25s%.2f%%\n" ,"GC%:",$GC_cotent * 100);
printf ("%-25s%d\n" ,"Longest:",$length[-1]);


## N50 :
my @rank = (10,20,30,40,50,60,70,80,90,100);
my $Nx = 0;
my $j = 0;

for (my $i = @length-1; $i >= 0; $i --){
	$Nx += $length[$i];
	if ($Nx >= ($total_length*($rank[$j]/100))){
		printf ("%-25s%d\n", "N$rank[$j](bp):", $length[$i]);
		printf ("%-25s%d\n", "L$rank[$j]:"    , @length-$i);
		$j++;
	}
	last if $j == @rank;
}

close IN;
#close OUT;

