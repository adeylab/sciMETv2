#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("P:A:m:O:", \%opt);

#defaults
$min = 50;

$die = "
Usage:

sciMET_annotPairwise.pl -P [pairwise file] -A [annot file] -O [output prefix]

Produces all-by-all cell overlap and shared cytosine methylation.

Options:

-P   [STR]   Pairwise methylation call file (req)
-A   [STR]   Annot file (req)
-m   [INT]   Min overlapping sites of cell pair (def = $min)

";

# PARSE OPTS

if (!defined $opt{'P'} ||
	!defined $opt{'A'} ||
	!defined $opt{'O'}) {die $die};
if (defined $opt{'m'}) {$min = $opt{'m'}};
	
open IN, "$opt{'A'}";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot{$cellID} = $annot;
} close IN;

open IN, "$opt{'P'}";
open OUT, ">$opt{'O'}.fracs.txt";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $CELLID_annot{$P[0]} && defined $CELLID_annot{$P[1]} && $P[3] >= $min) {
		print OUT "$CELLID_annot{$P[0]}\t$CELLID_annot{$P[1]}\t$P[4]\n";
		$ANNOT1_ANNOT2_ct{$CELLID_annot{$P[0]}}{$CELLID_annot{$P[1]}}++;
		$ANNOT1_ANNOT2_sum{$CELLID_annot{$P[0]}}{$CELLID_annot{$P[1]}}+=$P[4];
	}
} close IN; close OUT;

foreach $annot1 (keys %ANNOT1_ANNOT2_ct) {
	foreach $annot2 (keys %{$ANNOT1_ANNOT2_ct{$annot1}}) {
		$ANNOT1_ANNOT2_mean{$annot1}{$annot2} = $ANNOT1_ANNOT2_sum{$annot1}{$annot2}/$ANNOT1_ANNOT2_ct{$annot1}{$annot2};
	}
}

open IN, "$opt{'P'}";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if (defined $CELLID_annot{$P[0]} && defined $CELLID_annot{$P[1]} && $P[3] >= $min) {
		$ANNOT1_ANNOT2_stdevSum{$CELLID_annot{$P[0]}}{$CELLID_annot{$P[1]}}+=(($P[3]-$ANNOT1_ANNOT2_mean{$CELLID_annot{$P[0]}}{$CELLID_annot{$P[1]}})**2);
	}
} close IN;

open OUT, ">$opt{'O'}.summary.txt";
foreach $annot1 (keys %ANNOT1_ANNOT2_ct) {
	foreach $annot2 (keys %{$ANNOT1_ANNOT2_ct{$annot1}}) {
		$stdev = sprintf("%.3f", sqrt($ANNOT1_ANNOT2_stdevSum{$annot1}{$annot2}/($ANNOT1_ANNOT2_ct{$annot1}{$annot2}-1)));
		$mean = sprintf("%.3f", $ANNOT1_ANNOT2_mean{$annot1}{$annot2});
		print OUT "$annot1\t$annot2\t$mean\t$stdev\n";
	}
} close OUT;
