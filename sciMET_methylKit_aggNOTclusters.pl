#!/usr/bin/perl

$die = "

ARGV0-N = methylKit files for each cluster
Outputs a file with the same name, but adds NOT before methylKit.txt
this is for each cluster - to be used as a comparison in
methylKit, loading the cluster and the NOT aggregate cells
then performing the DMR analysis

Min cov is set to 10 (hardcoded)

";


#chrBase chr     base    strand  coverage        freqC   freqT
#chr1.25000      chr1    25000   F       1139    1.58    98.42
#chr1.75000      chr1    75000   F       370     2.70    97.30

if (!defined $ARGV[2]) {die $die};

for ($i = 0; $i < @ARGV; $i++) {
	%POS_cov = (); %POS_C = ();
	for ($j = 0; $j < @ARGV; $j++) {
		if ($i != $j) {
			open IN, "$ARGV[$j]";
			$null = <IN>;
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				$POS_cov{$P[0]} += $P[4];
				$POS_C{$P[0]} += ($P[4]*($P[5]/100));
			} close IN;
		}
	}
	open IN, "$ARGV[$i]";
	$out = $ARGV[$i]; $out =~ s/methylKit.txt/NOT.methylKit.txt/;
	open OUT, ">$out";
	$header = <IN>; print OUT "$header";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if (defined $POS_cov{$P[0]} && $POS_cov{$P[0]} >= 10 && $P[4] >= 10) {
			$freqC = sprintf("%.2f", (($POS_C{$P[0]}/$POS_cov{$P[0]})*100));
			$freqT = sprintf("%.2f", (100-$freqC));
			print OUT "$P[0]\t$P[1]\t$P[2]\t$P[3]\t$POS_cov{$P[0]}\t$freqC\t$freqT\n";
		}
	} close IN; close OUT;
}
