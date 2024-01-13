#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("M:C:O:A:", \%opt);

$die = "

sciMET_mtx2methylKitAnnot.pl -M [methylation percent matrix] -C [coverage matrix] -O [output prefix] -A [cluster annot file]

Gets methylation calls into a format that methylKit can use for DMR detection on a
cluster level.

Options:
   -O   [STR]   Output prefix (req)
   -M   [STR]   Methylation matrix (req)
   -C   [STR]   Coverage matrix
   -A   [STR]   Cluster annotation file (req)
   
";


if (!defined $opt{'M'} ||
	!defined $opt{'C'} ||
	!defined $opt{'A'} ||
	!defined $opt{'O'}) {die $die};
	
open IN, "$opt{'A'}";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot{$cellID} = $annot;
	$ANNOTS{$annot}++;
} close IN;

$header = "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT";
foreach $annot (keys %ANNOTS) {
	open OUT, ">$opt{'O'}.$annot.methylKit.txt"; print OUT "$header\n"; close OUT;
	$handle = $annot;
#	open $handle, "| sort -k 2,2 -k 3,3n >> $opt{'O'}.$annot.methylKit.txt";
	open $handle, ">>$opt{'O'}.$annot.methylKit.txt";
}

open M, "$opt{'M'}";
$cellIDs = <M>; chomp $cellIDs; @CELLIDS = split(/\t/, $cellIDs);

open C, "$opt{'C'}";
$null = <C>;

while ($l = <M>) {
	chomp $l; @M = split(/\t/, $l); $site = shift(@M);
	$c = <C>; chomp $c; @C = split(/\t/, $c); $null = shift(@C);
	
	($chr,$start,$end) = split(/_/, $site);
	$pos = int(($end+$start)/2);
	
	foreach $annot (keys %ANNOTS) {
		$ANNOT_cov{$annot} = 0;
		$ANNOT_C{$annot} = 0;
	}
	
	for ($i = 0; $i < @CELLIDS; $i++) {
		if (defined $CELLID_annot{$CELLIDS[$i]}) {
			$annot = $CELLID_annot{$CELLIDS[$i]};
			
			$ANNOT_cov{$annot} += $C[$i];
			$ANNOT_C{$annot} += ($C[$i]*$M[$i]);
			
		}
	}
	
	foreach $annot (keys %ANNOTS) {
		
		if ($ANNOT_cov{$annot}>0) {
			$freqC = sprintf("%.2f", (($ANNOT_C{$annot}/$ANNOT_cov{$annot})*100));
			$freqT = sprintf("%.2f", (100-$freqC));
			$handle = $annot;
			print $handle "$chr.$pos\t$chr\t$pos\tF\t$ANNOT_cov{$annot}\t$freqC\t$freqT\n";
		}
	}
}

foreach $annot (keys %ANNOTS) {
	$handle = $annot;
	close $handle;
}

