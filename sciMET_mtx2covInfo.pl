#!/usr/bin/perl

$die = "

sciMET_mtx2covInfo.pl [cov mtx file] > [output stats on cells]

";

if (!defined $ARGV[0]) {die $die};

open IN, "$ARGV[0]";
$h = <IN>; @H = split(/\t/, $h);
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$siteID = shift(@P);
	for ($i = 0; $i < @H; $i++) {
		if ($P[$i] > 0) {
			$CELLID_cov{$H[$i]}++;
		}
	}
	$siteCT++;
} close IN;

foreach $cellID (keys %CELLID_cov) {
	if ($CELLID_cov{$cellID}>0) {
		$frac = sprintf("%.4f", $CELLID_cov{$cellID}/$siteCT);
	}
	print "$cellID\t$CELLID_cov{$cellID}\t$frac\n";
}
