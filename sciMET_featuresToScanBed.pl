#!/usr/bin/perl

#defaults
$featureWin = 100;
$minSize = 1000;
$maxSize = 50000;
$upstreamSize = 0;
$upstreamWin = 0;
$downstreamSize = 0;
$downstreamWin = 0;

use Getopt::Std; %opt = ();

getopts("N:m:M:S:s:E:e:", \%opt);

$die = "

sciMET_featuresToScanBed.pl (options) [Bed of peaks or features] [output bed file]

Defaults set for peaks. Use -S -E -s -e for genes.
It is recommended to sort the bed file after generation.

Options:

-N   [INT]   Number of windows within feature (def = $featureWin)

-m   [INT]   Minimum size of feature (def = $minSize)
-M   [INT]   Max size of feature (def = $maxSize)

-S   [INT]   Bases upstream of feature (def = $upstreamSize)
-s   [INT]   Number of windows in upstream region (def = $upstreamWin)

-E   [INT]   Bases downstream of feature (def = $downstreamSize)
-e   [INT]   Number of windows in downstream region (def = $downstreamWin)

";

if (!defined $ARGV[1]) {die $die};

if (defined $opt{'N'}) {$featureWin = $opt{'N'}};
if (defined $opt{'M'}) {$maxSize = $opt{'M'}};
if (defined $opt{'m'}) {$minSize = $opt{'m'}};
if (defined $opt{'S'}) {$upstreamSize = $opt{'S'}};
if (defined $opt{'s'}) {$upstreamWin = $opt{'s'}};
if (defined $opt{'E'}) {$downstreamSize = $opt{'E'}};
if (defined $opt{'e'}) {$downstreamWin = $opt{'e'}};

if ($upstreamSize > 1) {
	$upstreamIncrement = int($upstreamSize/$upstreamWin);
}

if ($downstreamSize > 1) {
	$downstreamIncrement = int($downstreamSize/$downstreamWin);
}

open OUT, ">$ARGV[1]";
open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winID = 0;
	if (defined $P[5] && $P[5] eq "-") { #strand info and is reverse
		if ($upstreamSize > 1) {
			$start = $P[2]+$upstreamSize;
			for ($i = 0; $i < $upstreamWin; $i++) {
				$end = $start-$upstreamIncrement+1;
				print OUT "$P[0]\t$end\t$start\t$winID\n";
				$winID++;
				$start = $end-1;
			}
		}
		$featureSize = $P[2]-$P[1];
		if ($featureSize < $minSize) {
			$start = int((($P[1]+$P[2])/2)+($minSize/2));
			$P[1] = int((($P[1]+$P[2])/2)-($minSize/2));
			$featureSize = $minSize;
		} elsif ($featureSize > $maxSize) {
			$start = int((($P[1]+$P[2])/2)+($maxSize/2));
			$P[1] = int((($P[1]+$P[2])/2)-($maxSize/2));
			$featureSize = $maxSize;
		} else {
			$start = $P[2];
		}
		$featureIncrement = int($featureSize/$featureWin);
		for ($i = 0; $i < $featureWin; $i++) {
			$end = $start-$featureIncrement+1;
			print OUT "$P[0]\t$end\t$start\t$winID\n";
			$winID++;
			$start = $end-1;
		}
		if ($downstreamSize > 1) {
			$start = $P[1];
			for ($i = 0; $i < $downstreamWin; $i++) {
				$end = $start-$downstreamIncrement+1;
				if ($start>0&&$end>0) {
					print OUT "$P[0]\t$end\t$start\t$winID\n";
				}
				$winID++;
				$start = $end-1;
			}
		}
	} else { # no strand info or forward
		if ($upstreamSize > 1) {
			$start = $P[1]-$upstreamSize;
			for ($i = 0; $i < $upstreamWin; $i++) {
				$end = $start+$upstreamIncrement-1;
				if ($start>0&&$end>0) {
					print OUT "$P[0]\t$start\t$end\t$winID\n";
				}
				$winID++;
				$start = $end+1;
			}
		}
		$featureSize = $P[2]-$P[1];
		if ($featureSize < $minSize) {
			$start = int((($P[1]+$P[2])/2)-($minSize/2));
			$P[2] = int((($P[1]+$P[2])/2)+($minSize/2));
			$featureSize = $minSize;
		} elsif ($featureSize > $maxSize) {
			$start = int((($P[1]+$P[2])/2)-($maxSize/2));
			$P[2] = int((($P[1]+$P[2])/2)+($maxSize/2));
			$featureSize = $maxSize;
		} else {
			$start = $P[1];
		}
		$featureIncrement = int($featureSize/$featureWin);
		for ($i = 0; $i < $featureWin; $i++) {
			$end = $start+$featureIncrement-1;
			print OUT "$P[0]\t$start\t$end\t$winID\n";
			$winID++;
			$start = $end+1;
		}
		if ($downstreamSize > 1) {
			$start = $P[2];
			for ($i = 0; $i < $downstreamWin; $i++) {
				$end = $start+$downstreamIncrement-1;
				print OUT "$P[0]\t$start\t$end\t$winID\n";
				$winID++;
				$start = $end+1;
			}
		}
	}
} close IN; close OUT;
exit;