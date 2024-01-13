#!/usr/bin/perl

#chr10   101288606       66.67   1       2

$die = "

ARGV0 = cellcall folder
ARGV1 = annot of clusters
ARGV2 = out folder
ARGV3 = CG or CH (def = CG)

will generate a bed of sites covered by cluster

";

if (!defined $ARGV[2]) {die $die};

if (defined $ARGV[3]) {$context = $ARGV[3]} else {$context = "CG"}; 

open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$ANNOT_CELLID{$annot}{$cellID} = 1;
} close IN;

system("mkdir $ARGV[2]");

open STATS, ">$ARGV[2]/$ARGV[2].stats";
foreach $annot (keys %ANNOT_CELLID) {
	%SITE_cov = (); $annot_cov = 0;
	foreach $cellID (keys %{$ANNOT_CELLID{$annot}}) {
		open IN, "$ARGV[0]/$cellID.$context.cov";
		while ($l = <IN>) {
			chomp $l;
			($chr,$pos,$meth,$C,$mC) = split(/\t/, $l);
			$site = "$chr.$pos";
			if (!defined $SITE_cov{$site}) {$annot_cov++};
			$SITE_cov{$site} += ($C+$mC);
			$SITE_meth{$site} = $meth;
		} close IN;
	}
	print STATS "$annot\r$annot_cov\n";
	open OUT, "| sort -k1,1 -k2,2n > $ARGV[2]/$annot.meth.bed";
	foreach $site (keys %SITE_cov) {
		($chr,$pos) = split(/\./, $site);
		print OUT "$chr\t$pos\t$pos\t$SITE_cov{$site}\t$SITE_meth{$site}\n";
	} close OUT;
} close STATS;
