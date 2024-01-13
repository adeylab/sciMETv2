#!/usr/bin/perl

$die = "

sciMET_methylKit2bgph.pl [input methylKit window file] [window size] [min cov] > [out bedgraph]

";

if (!defined $ARGV[2]) {die $die};

open IN, "$ARGV[0]";

$null = <IN>;

#chrBase 		 chr     base    strand  coverage  freqC   freqT
#chr1.11000      chr1    11000   F       90        35.55   64.45

$half = int($ARGV[1]/2);

while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[4]>=$ARGV[2]) {
		$start = $P[2]-$half;
		$end = $P[2]+$half;
		print "$P[1]\t$start\t$end\t$P[5]\n";
	}
} close IN;
