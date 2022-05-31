#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:t:q:", \%opt);

$threads = 1;
$minq = 10;

$die = "

sciMET_rmdup_pe.pl (options) [namesort bam file]

Performs barcode-aware duplicate read removal.
Paired end alignment mode only.

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

open OUT, "| samtools view -bSu - | samtools sort -@ $threads -T $opt{'O'}.TMP -m 4G - > $opt{'O'}.bbrd.q10.bam";

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
	$tag = $P[0]; $tag =~ s/\#.+$//;
	$barc = $P[0]; $barc =~ s/:.+$//;
	
	if ($currentBarc ne $barc) {
		# print stored output
		if ($currentBarc ne "null") {
			foreach $tag (keys %READ_keep) {
				print OUT "$READ_keep{$tag}\n";
				$total_kept++;
			}
		}
		# reset hashes
		%COORD_read = (); %READ_keep = ();
		$currentBarc = $barc;
	}
	
	# process read
	$BARC_total{$barc}++;
	
	$start = $P[3];
	$coord = "$P[2].$start";
	
	if (defined $READ_keep{$tag}) { # mate found, keep it
		$tag .= ".m"; # make tag unique
		$READ_keep{$tag} = $l;
		$COORD_read{$coord} = 1;
		$BARC_kept{$barc}++;
	}
	
	if (!defined $COORD_read{$coord}) { # coord not observed, keep it
		$READ_keep{$tag} = $l;
		$COORD_read{$coord} = 1;
		$BARC_kept{$barc}++;
	}
} close IN; close OUT;

open OUT, ">$opt{'O'}.complexity.txt";
$rank = 1;
foreach $barc (sort {$BARC_kept{$b}<=>$BARC_kept{$a}} keys %BARC_kept) {
	$pct = sprintf("%.2f", ($BARC_kept{$barc}/$BARC_total{$barc})*100);
	print OUT "$rank\t$barc\t$BARC_total{$barc}\t$BARC_kept{$barc}\t$pct\n";
	$rank++;
} close OUT;