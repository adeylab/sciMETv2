#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:t:q:", \%opt);

$threads = 1;
$minq = 10;

$die = "

sciMET_rmdup.pl (options) [namesort bam file]

Performs barcode-aware duplicate read removal.
Duplicate reads will be chosen based on the longest
aligned read.

Options:
   -O   [STR]   Output prefix (def = input bam prefix)
   -t   [INT]   Threads for the sorting process. (def = $threads)
   -q   [INT]   Min read alignment quality (def = $minq)
   
";

if (!defined $ARGV[0]) {die $die};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.bam$//;
	$opt{'O'} =~ s/\.nsrt$//;
}

if (defined $opt{'q'}) {$minq = $opt{'q'}};
if (defined $opt{'t'}) {$threads = $opt{'t'}};

open OUT, "| samtools view -bSu - | samtools sort -@ $threads -T $opt{'O'}.TMP - > $opt{'O'}.bbrd.q10.bam";

open HEAD, "samtools view -H $ARGV[0] |";
while ($l = <HEAD>){print OUT "$l"};
close HEAD;

open IN, "samtools view -q $minq $ARGV[0] |";

$currentBarc = "null";
$total_aligned = 0;
$total_kept = 0;

while ($l = <IN>) {
	$total_aligned++;
	
	chomp $l;
	@P = split(/\t/, $l);
	$barc = $P[0]; $barc =~ s/:.+$//;
	
	if ($currentBarc ne $barc) {
		# print stored output
		if ($currentBarc ne "null") {
			foreach $coord (keys %COORD_read) {
				print OUT "$COORD_read{$coord}\n";
				$total_kept++;
			}
		}
		# reset hashes
		%COORD_read = (); %COORD_length = ();
		$currentBarc = $barc;
	}
	
	# process read
	$BARC_total{$barc}++;
	
	$length = length($P[9]);
	$XR = chop $P[14]; $XG = chop $P[15];
	
	# flags different w/ bismark; reads aligning to - strand of reference
	# have the start of the sequence read as the end of the alignment.
	# this matters if read lengths sequenced are different, where the
	# alignment coordinate will be the end of the read.
	
	if ($XR eq $XG) {	# read is reference-strand and thus +
		$start = $P[3];
	} else {			# read is not reference strand and thus - (tag/prime site is end of alignment)
		$start = $P[3]+$length;
	}
	
	$coord = "$P[2].$start.$XR.$XG";
	
	
	if (!defined $COORD_read{$coord}) {
		$COORD_read{$coord} = $l;
		$COORD_length{$coord} = $length;
		$BARC_kept{$barc}++;
	} else {
		if ($COORD_length{$coord} < $length) { # retain longer alignment
			$COORD_length{$coord} = $length;
			$COORD_read{$coord} = $l;
		}
	}
} close IN; close OUT;

open OUT, ">$opt{'O'}.complexity.txt";
$rank = 1;
foreach $barc (sort {$BARC_kept{$b}<=>$BARC_kept{$a}} keys %BARC_kept) {
	$pct = sprintf("%.2f", ($BARC_kept{$barc}/$BARC_total{$barc})*100);
	print OUT "$rank\t$barc\t$BARC_total{$barc}\t$BARC_kept{$barc}\t$pct\n";
	$rank++;
} close OUT;