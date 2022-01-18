#! /usr/bin/perl -w

use strict;

my $id  = '';
my $seq = '';
while (<>) {
	chomp;
	if ( /^>(\S+)/ ) {
		
		if ( $seq ) {
			while ( $seq =~ /(n+)/gi ) {
				my $Nlength = length $1;
				my $Nend = pos $seq;
				my $Nstart = $Nend - $Nlength + 1;
				print join("\t", $id, $Nstart, $Nend),"\n";
			}
		}
		
		$id  = $1;
		$seq = '';
	}else {
		$seq .= $_;
	}
}