#! /usr/bin/perl -w

use strict;
use Data::Dumper;

my $task = 'GAEP';

my $dir = shift;
my $prefix = shift;
my $prefix_out = "$dir/$prefix";
my $checkpoint = shift;

my %pipe = (
    0x0001 => "${prefix_out}_genome_stat.txt",                  #"basic",
    0x0002 => "$dir/mapping/tgs_mapping/flagstat_tgs.txt",      #"tmr", 
    0x0004 => "$dir/breakpoint_output/${prefix}_breakpoints.txt",                  #"bkp",
    0x0008 => "$dir/mapping/ngs_mapping/flagstat_ngs.txt",      #"nmr",
    0x0010 => "$dir/snv_depth_tgs/$prefix.png",                 #"tsd plot",
    0x0020 => "$dir/snv_depth_ngs/$prefix.png",                 #"nsd plot",
    0x0040 => "$dir/mapping/trans_mapping/flagstat_trans.txt",  #"trmr",
    0x0080 => "$dir/busco_out/$prefix/short_summary.*$prefix.txt",     #"busco",
	0x0100 => "$prefix_out.base.qv",
	0x0200 => "$dir/merqury_output/$prefix.qv",
);
$pipe{0x0080} = glob $pipe{0x0080};

my %result = (
	l1 => {},
	l2 => {},
	l3 => {},
);
my $sinfo  = 0;

if ($checkpoint & 0x0001 == 0x0001 && -e $pipe{0x0001}) {
	open STAT, '<', $pipe{0x0001} or die "[$task] Can't find $pipe{0x0001}.\n";
	while (<STAT>) {
	#while (<>) {
		chomp;
		my ($para, $value) = split /:\s+/, $_;
		if ($para eq 'Total length') {
			$result{l1}->{1} = [$para, 1];
			$result{l3}->{1} = [$value, 1];
		}
		if ($para eq 'Total sequences') {
			$result{l1}->{2} = [$para, 1];
			$result{l3}->{2} = [$value, 1];
		}
		if ($para eq 'Ungapped length') {
			$result{l1}->{3} = [$para, 1];
			$result{l3}->{3} = [$value, 1];
		}
		if ($para eq 'N50') {
			$result{l1}->{4} = [$para, 1];
			$result{l3}->{4} = [$value, 1];
		}
		if ($para eq 'GC%') {
			$result{l1}->{5} = [$para, 1];
			$result{l3}->{5} = [$value, 1];
		}

		if ($para eq 'N/Gb') {
			$result{l1}->{6} = [$para, 1];
			$result{l3}->{6} = [$value, 1];
		}
	}
	close STAT;
}

if ($checkpoint & 0x0002 == 0x0002 && -e $pipe{0x0002}) {
	my $para  = 'TGS';
	my $value = '';
	open STAT, '<', $pipe{0x0002} or warn "[$task] Can't find $pipe{0x0002}.\n";
	while (<STAT>) {
		$value = $1 if (/mapped\W*?\((\S+%)/);
	}
	close STAT;
	if ($value) {
		$result{l2}->{2} = [$para, 1];
		$result{l3}->{10} = [$value, 1];
	}
	
}

if ($checkpoint & 0x0004 == 0x0004 && -e $pipe{0x0004}) {
	my $para  = 'Breakpoints per Mb';
	my $value = 0;
	open STAT, '<', $pipe{0x0004} or warn "[$task] Can't find $pipe{0x0004}.\n";
	while (<STAT>) {
		next if (/^#/);
		my @line = split;
		$value += 1 if ($line[-1] > 0);
	}
	close STAT;
	if ($result{l3}->{1}->[0] != 0) {
		$value = $value / $result{l3}->{1}->[0] * 1000000;
		$value = sprintf("%0.2f", $value);
		$result{l1}->{8} = [$para, 1];
		$result{l3}->{8} = [$value, 1];
	}else {
		warn "Can't get total length.\n";
	}
}

if ($checkpoint & 0x0008 == 0x0008 && -e $pipe{0x0008}) {
	my $para  = 'NGS';
	my $value = '';
	open STAT, '<', $pipe{0x0008} or warn "[$task] Can't find $pipe{0x0008}.\n";
	while (<STAT>) {
		$value = $1 if (/mapped\W*?\((\S+%)/);
	}
	close STAT;
	if ($value) {
		$result{l2}->{1} = [$para, 1];
		$result{l3}->{9} = [$value, 1];
	}
}

if ($checkpoint & 0x0040 == 0x0040 && -e $pipe{0x0040}) {
	my $para  = 'Trans';
	my $value = '';
	open STAT, '<', $pipe{0x0040} or warn "[$task] Can't find $pipe{0x0040}.\n";
	while (<STAT>) {
		$value = $1 if (/mapped\W*?\((\S+%)/);
	}
	close STAT; 
	if ($value) {
		$result{l2}->{3}  = [$para, 1];
		$result{l3}->{11} = [$value, 1];	
	}
}

if ($checkpoint & 0x0080 == 0x0080 && -e $pipe{0x0080}) {
	my ($C, $D, $F) = (0,0,0);
	open STAT, '<', $pipe{0x0080} or warn "[$task] Can't find $pipe{0x0080}.\n";
	while (<STAT>) {
		next if (/^#/);
		($C, $D, $F) = ($1, $2, $3) if (/C:(\S+?%).*?D:(\S+?%).*?F:(\S+?%)/);
	}
	close STAT;
	$result{l2}->{4}  = ["Complete", 1];
	$result{l3}->{12} = [$C, 1];
	$result{l2}->{5}  = ["Fragmented", 1];
	$result{l3}->{13} = [$F, 1];
	$result{l2}->{6}  = ["Duplicated", 1];
	$result{l3}->{14} = [$D, 1];
}

if ($checkpoint & 0x0100 == 0x0100 && -e $pipe{0x0100}) {
	my $para  = 'QV (reads mapping)';
	my $value = 0;
	
	open STAT, '<', $pipe{0x0100} or warn "[$task] Can't find $pipe{0x0100}.\n";
	while (<STAT>) {
		if (/^QV/) {
			$value = (split)[1];
		}
	}
	close STAT; 
	
	if ($value) {
		$result{l1}->{6} = [$para, 1];
		$result{l3}->{6} = [$value, 1];	
	}
}

if ($checkpoint & 0x0200 == 0x0200 && -e $pipe{0x0200}) {
	my $para  = 'QV (K-mer)';
	my $value = 0;
	
	open STAT, '<', $pipe{0x0200} or warn "[$task] Can't find $pipe{0x0200}.\n";
	while (<STAT>) {
		if (/\S+/) {
			$value = (split)[3];
		}
	}
	close STAT; 
	if ($value) {
		$result{l1}->{7} = [$para, 1];
		$result{l3}->{7} = [$value, 1];	
	}
}

my $cl = 0;
for (1..3) {
	$cl += 1 if $result{l2}->{$_};
}
$result{l1}->{9} = ["Mapping ratio", $cl] if $cl > 0;

if ($result{l2}->{4} && $result{l2}->{5} && $result{l2}->{6}) {
	$result{l1}->{10} = ["Busco", 3];
}


#print Dumper(\%result);		
create_table(\%result);

#add_a_element(1, 1, "Total length");

##---------------------------------subroutine---------------------------------##

sub head {
	my $date = `date`;
	$date =~ s/\n$//;
	print "<!DOCTYPE html>\n\n";
	print "<html>\n";
	print "  <head>\n";
	print "    <meta charset=\"UTF-8\">\n";
#	print "    <title>GAEP: Genome Assembly Assessment Pipeline</title>\n";
	print "  </head>\n\n";
	
	print "  <body>\n";
	print "    <table width=\"100%\" border=\"1\" cellspacing=\"0\" cellpadding=\"0\" align=\"left\">\n";
	print "      <tr>\n";
	print "        <td align=\"center\" class=\"title\" height=\"60\">GAEP: Genome Assembly Evaluating Pipeline</td>\n";
	print "      </tr>\n";
	print "      <tr>\n";
	print "        <td align=\"right\" height=\"25\">",$date,"</td>\n";
	print "      </tr>\n";
	print "    </table>\n";
	return;
}

sub table_start {
	print "    <table width=\"100%\" border=\"1\" cellspacing=\"0\" cellpadding=\"4\"\ align=\"center\">\n";
	return;
}

sub line_start {
	print "      <tr>\n";
	return;
}

sub add_a_element {
	my $colspan = shift;
	my $rowspan = shift;
	my $element = shift;
	my $align   = "center";
	print "        <td ";
	print "align=", "\"$align\"";
	print " colspan=", "\"$colspan\"";
	print " rowspan=", "\"$rowspan\"";
	print ">";
	print "$element";
	print "</td>\n";
	return;
}

sub line_end {
	print "      </tr>\n";
	return;
}

sub table_end {
	print"    </table>\n";
	return;
}

sub tail {
	print "  </body>\n";
	print "</html>\n";
	return;
}

sub create_table {
	my $table = shift;
	head();
	table_start();

	for my $layer (sort keys %$table) {
		line_start();
		
		for (sort {$a <=> $b} keys %{$table->{$layer}}) {
			my $ele = $table->{$layer}->{$_}->[0];
			my $col = $table->{$layer}->{$_}->[1];
			my $row = 1;
			if ($layer eq 'l1') {
				if ($table->{$layer}->{$_}->[0] ne 'Mapping ratio' && $table->{$layer}->{$_}->[0] ne 'Busco' ){
					$row = 2;
				}
			}
			add_a_element($col, $row, $ele);
		}
		
		line_end();
	}

	table_end();
	tail();
}
