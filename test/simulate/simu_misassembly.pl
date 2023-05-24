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
	$fai{$line[0]} = [@line[1,2]];
	push @ctgList, $line[0];
}
close $fai;
@ctgList = grep { $fai{$_}->[0] >= 500000 }@ctgList;


my $ref_adjust = 0;
my $simu_adjust = 0;
my $prectg = '';
my $id_seq = 1;
open POSI, '<', "position_redun.txt" || die $!;
while (<POSI>) {
	chomp;
	my ($ctg,$posi,$length,$type) = split;
	if ( $prectg ne $ctg ) {
	$ref_adjust = 0; 
	$simu_adjust = 0;
	}
	my $posi_ref = $posi + $ref_adjust;
	my $posi_simu = $posi + $simu_adjust;
	say $ref_adjust; 
	say $simu_adjust;
	
	if ( $type eq 'expanded_3' ) {
		my $repeat = substr($hash{$ctg},$posi_ref,$length);
		my $trans_len = rand_r(0,1000);
		my $transition = substr($hash{$ctg},$posi_ref+$length,$trans_len);
		my $repeat_var = induce_variants($repeat,$length);
		#ref
		insert_seq($repeat_var, $hash{$ctg}, $posi_ref+$length+(length $transition));
		#simu
		insert_seq($repeat_var.$transition.$repeat_var, $hash_simu{$ctg}, $posi_simu+$length+(length $transition));
		
		say join("\t",$id_seq,1,$ctg,$type,$length,$posi_ref,$posi_ref+$length+(length $transition)+(length $repeat_var));
		say join("\t",$id_seq,2,$ctg,$type,$length,$posi_simu,$posi_simu+$length+(length $transition)*2+(length $repeat_var)*2,
		              $posi_simu+$length+(length $transition)+(length $repeat_var),$posi_simu+$length+(length $transition)*2+(length $repeat_var));
		
		$ref_adjust += length $repeat_var;
		$simu_adjust += length ($repeat_var.$transition.$repeat_var);
	}
	
	if ( $type eq 'expanded_4' ) {
		my $repeat = substr($hash{$ctg},$posi_ref,$length);
		my $trans_len = rand_r(0,1000);
		my $transition = substr($hash{$ctg},$posi_ref+$length,$trans_len);
		my $repeat_var = induce_variants($repeat,$length);
		#ref
		insert_seq($repeat_var, $hash{$ctg}, $posi_ref+$length+(length $transition));
		#simu
		my $ins_seq = $repeat_var.$transition.$repeat.$transition.$repeat_var;
		insert_seq($ins_seq, $hash_simu{$ctg}, $posi_simu+$length+(length $transition));
		
		say join("\t",$id_seq,1,$ctg,$type,$length,$posi_ref,$posi_ref+$length+(length $transition)+(length $repeat_var));
		say join("\t",$id_seq,2,$ctg,$type,$length,$posi_simu,$posi_simu+$length*2+(length $transition)*3+(length $repeat_var)*2,
		              $posi_simu+$length+(length $transition)+(length $repeat_var),$posi_simu+$length+(length $transition)*2+(length $repeat_var));
		
		$ref_adjust += length $repeat_var;
		$simu_adjust += length $ins_seq;
	}
	
	if ( $type eq 'collasped' ) {
		my $repeat = substr($hash{$ctg},$posi_ref,$length);
		my $trans_len = rand_r(0,1000);
		my $transition = substr($hash{$ctg},$posi_ref+$length,$trans_len);
		my $repeat_var = induce_variants($repeat,$length);
		#ref
		insert_seq($repeat_var, $hash{$ctg}, $posi_ref+$length+(length $transition));
		#simu
		delete_seq($hash_simu{$ctg}, $posi_simu+$length,(length $transition));
		
		say join("\t",$id_seq,1,$ctg,$type,$length,$posi_ref,$posi_ref+$length+(length $transition)+(length $repeat_var));
		say join("\t",$id_seq,2,$ctg,$type,$length,$posi_simu,$posi_simu+$length,$posi_simu,$posi_simu+$length);
		
		$ref_adjust += (length $repeat_var);
		$simu_adjust -= (length $transition);
	}
	
	if ( $type eq 'inv' ) {
		my $repeat = substr($hash{$ctg},$posi_ref,$length);
		my $trans_len = rand_r(1000,10000);
		my $repeat_var = induce_variants($repeat,$length);
		$repeat_var = reverse $repeat_var;
		#$repeat_var =~ tr/atcgATCG/tagcTAGC/;
		#ref
		insert_seq($repeat_var, $hash{$ctg}, $posi_ref+$length+$trans_len);
		##simu
		insert_seq($repeat_var, $hash_simu{$ctg}, $posi_simu+$length+$trans_len);
		my $ctg_len = length $hash_simu{$ctg};
		seq_comp($hash_simu{$ctg},$ctg_len,$posi_simu+$length,$trans_len);
		random_INV($hash_simu{$ctg},$ctg_len,$posi_simu,$length+$trans_len+(length $repeat_var));
		#
		say join("\t",$id_seq,1,$ctg,$type,$length,$trans_len,$posi_ref,$posi_ref+$length+$trans_len+(length $repeat_var));
		say join("\t",$id_seq,2,$ctg,$type,$length,$trans_len,$posi_simu,$posi_simu+$length+$trans_len+(length $repeat_var), $posi_simu+(length $repeat_var),$posi_simu+(length $repeat_var)+$trans_len);
		#
		$ref_adjust += length $repeat_var;
		$simu_adjust += length $repeat_var;
	}
	
	if ( $type eq 'ins' ) {
		my $insert_seq = random_seq($length);
		
		#ref
		
		#simu
		insert_seq($insert_seq, $hash_simu{$ctg}, $posi_simu);
		
		say join("\t",$id_seq,1,$ctg,$type,$length,$posi_ref,$posi_ref+1);
		say join("\t",$id_seq,2,$ctg,$type,$length,$posi_simu,$posi_simu+$length);
		
		$simu_adjust += length $insert_seq;
	}
	
	if ( $type eq 'del' ) {

		#ref
		
		#simu
		delete_seq($hash_simu{$ctg}, $posi_simu,$length);
		
		say join("\t",$id_seq,1,$ctg,$type,$length,$posi_ref,$posi_ref+$length);
		say join("\t",$id_seq,2,$ctg,$type,$length,$posi_simu,$posi_simu+1);
		
		$simu_adjust -= $length;
	}
	
	$id_seq += 1;
	$prectg = $ctg;
}


open REF, '>', $fa.$$.".ref.fasta" or die $!;
for my $ctg (@ctgList) {
	print REF ">$ctg\n";
	print REF $hash{$ctg},"\n";
}
close REF;

open SIMU, '>', $fa.$$.".simu.fasta" or die $!;
for my $ctg (@ctgList) {
	print SIMU ">$ctg\n";
	print SIMU $hash_simu{$ctg},"\n";
}
close SIMU;


##repeats expanded


##repeats collasped


##inversion


#----------------------------------subroutine----------------------------------#

sub rand_r {
	my $start = shift;
	my $end = shift;
	return if $end <= $start;
	my $int = @_ ? shift : 1;
	return $int ? int (rand $end-$start+1) + $start : rand ($end-$start) + $start;
}

sub random_INV_rev {
	my $seq = \shift;
	my $seqLen = @_ ? shift @_ : length $$seq;
	my $invStart = shift ;
	my $invLen = shift;
	
	my $inv = substr($$seq, $invStart, $invLen);
	$inv = reverse $inv;
	my $left = substr($$seq, 0, $invStart);
	my $right = substr($$seq, $invStart+$invLen, $seqLen-($invStart+$invLen));
	$$seq = join('', $left, $inv, $right);
}

sub random_INV {
	my $seq = \shift;
	my $seqLen = @_ ? shift @_ : length $$seq;
	my $invStart = shift ;
	my $invLen = shift;
	
	my $inv = substr($$seq, $invStart, $invLen);
	$inv = reverse $inv;
	my $left = substr($$seq, 0, $invStart);
	my $right = substr($$seq, $invStart+$invLen, $seqLen-($invStart+$invLen));
	$$seq = join('', $left, $inv, $right);
}

sub seq_comp {
	my $seq = \shift;
	my $seqLen = @_ ? shift @_ : length $$seq;
	my $start = shift ;
	my $len = shift;
	
	my $comp = substr($$seq, $start, $len);
	$comp =~ tr/atcgATCG/tagcTAGC/;
	my $left = substr($$seq, 0, $start);
	my $right = substr($$seq, $start+$len, $seqLen-($start+$len));
	$$seq = join('', $left, $comp, $right);
}

sub delete_seq {                     #delete_seq(seq, start position, seq length)
	my $seq = \(shift @_);
	my $start = shift;
	my $len = shift;
	my $left = substr($$seq, 0, $start);
	my $right = substr($$seq, $start+$len, (length $$seq)-($start+$len));
	$$seq = join('', $left, $right);
}

sub insert_seq {                    #insert_seq(ins seq, target seq, start position)
	my $ins_seq = \shift;
	my $target_seq = \shift;
	my $start = shift;
	my $left = substr($$target_seq, 0, $start);
	my $right = substr($$target_seq, $start, (length $$target_seq) - $start);
	$$target_seq = join('', $left, $$ins_seq , $right);
}

sub random_seq {
	my $len = shift; 
	my $random_seq = '';
	for ( 1..$len ) {
		$random_seq .= $base[int(rand(4))];
	}
	return $random_seq;
}

sub random_base_subtitude {
	my $seq = \(shift);
	my $seqLen = length $$seq;
	my $ran_position = int(rand($seqLen)) + 1;
	
	my $left = substr($$seq, 0, $ran_position-1);
	my $right = substr($$seq, $ran_position, $seqLen - $ran_position);
	
	my $oriBase = substr($$seq,$ran_position-1,1);
	my $subBase = $base[int(rand(4))];
	$subBase = $oriBase eq $subBase ? '' : $subBase;

	$$seq = join('', $left, $subBase, $right);
}

sub induce_variants {                    #induce_variants(seq, seq length)
	my $seq = shift;
	my $seqLen = @_ ? shift @_ : length $seq;
	my $identity = rand(0.25) + 0.75;
	my $varThreshold = int($seqLen*(1-$identity));
	while (1) {
		last if $varThreshold <= 0;
		my $varType = int rand(7);
		if ( $varType < 5 ) { ##snv 0,1,2,3,4
			random_base_subtitude($seq);
			$varThreshold -= 1;
		}elsif ( $varType == 5 ) {  ##ins 5
			$seqLen = length $seq;
			my $ran_position = int(rand($seqLen)) + 1;
			my $insLen = int (rand(100)) + 1;
			insert_seq(random_seq($insLen),$seq,$ran_position);
			$varThreshold -= $insLen;
		}else {  ##del 6
			$seqLen = length $seq;
			my $ran_position = int (rand($seqLen)) + 1;
			my $delLen = int (rand(100)) + 1;
			$delLen = $seqLen - $ran_position if  $ran_position+$delLen > $seqLen;
			delete_seq($seq,$ran_position,$delLen);
			$varThreshold -= $delLen;
		}
	}
	return $seq;
}
