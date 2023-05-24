#! /usr/bin/perl -w

use strict;
use 5.010;


my $pre_line = [''];
while (<>) {
	next if /^\s*$/;
	chomp;
	my @line = split;
	if ($line[0] ne $pre_line->[0]) {
		$pre_line = \@line;
		next;
	}
	
	
	if (  $pre_line->[3] eq 'ins' && $line[1] <= ($pre_line->[1] + ($pre_line->[2]*2+50000) + 5000) ) {
		next;
	}elsif ( $line[1] <= ($pre_line->[1] + ($pre_line->[2]*3) + 5000) ) {
		next;
	}else {
		say join("\t",@{$pre_line});
		$pre_line = \@line;
	}
}

