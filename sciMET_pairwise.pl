#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("F:O:", \%opt);

#defaults


$die = "
Usage:

sciMET_pairwise.pl -F [chroms] -O [output file] 2> [logfile of progress]

Produces all-by-all cell overlap and shared cytosine methylation.

Options:

-F   [DIR]   Directory with methylation by cell for each chromosome (req)
             Can be a comma separated list to process multiple experiments.
             Or can be a file ending in 'txt' with aline for each folder
-O   [STR]   Output file (req)

";

# PARSE OPTS

if (!defined $opt{'F'} ||
	!defined $opt{'O'}) {die $die};
	
if ($opt{'F'} =~ /txt$/) {
	@FOLDERS = ();
	open F, "$opt{'F'}";
	while ($l = <F>) {chomp $l; push @FOLDERS, $l};
	close F;
} else {
	@FOLDERS = split(/,/, $opt{'F'});
}

%CELLID1_CELLID2_overlap = ();
%CELLID1_CELLID2_shared = ();

for ($chr = 1; $chr <= 22; $chr++) {
	read_chr($chr);
} read_chr("X");

sub read_chr {
	$in_chr = $_[0];
	$ts = localtime(time);
	print STDERR "$ts Chr $in_chr being read...\n";
	%POS_CELLID = ();
	%POS_STATUS = ();
	$pos_ct = 0;
	$pos_read = 0;
	foreach $chr_dir (@FOLDERS) {
		if (-e "$chr_dir/chr$in_chr.bed.gz") {
			open IN, "zcat $chr_dir/chr$in_chr.bed.gz |";
			$ts = localtime(time);
			print STDERR "$ts      File: $chr_dir/chr$in_chr.bed.gz\n";
			while ($l = <IN>) {
				chomp $l;
				($chrom,$start,$end,$cellID,$call) = split(/\t/, $l);
				if (!defined @{$POS_CELLID{$start}}[0]) {
					$pos_ct++;
				}
				push @{$POS_CELLID{$start}}, $cellID;
				push @{$POS_STATUS{$start}}, $call;
			} close IN;
		} else {print STDERR "$chr_dir/chr$in_chr.bed.gz does not exist. Skipping...\n"};
	}
	$ts = localtime(time);
	print STDERR "$ts Chr $in_chr read in, $pos_ct sites. Computing cell-cell overlaps and adding to totals.\n";
	$increment = 0.05; $report = $increment;
	foreach $pos (keys %POS_CELLID) {
		$pos_read++;
		if (($pos_read/$pos_ct) >= $report) {
			$ts = localtime(time);
			print STDERR "$ts      $report complete ($pos_read / $pos_ct)\n";
			$report += $increment;
		}
		for ($i = 0; $i < @{$POS_CELLID{$pos}}; $i++) { 
			for ($j = 0; $j < @{$POS_CELLID{$pos}}; $j++) {
				if ($i != $j) {
					$CELLID1_CELLID2_overlap{$POS_CELLID{$pos}[$i]}{$POS_CELLID{$pos}[$j]}++;
					if ($POS_STATUS{$pos}[$i] eq $POS_STATUS{$pos}[$j]) {
						$CELLID1_CELLID2_shared{$POS_CELLID{$pos}[$i]}{$POS_CELLID{$pos}[$j]}++;
					}
				}
			}
		}
	}
}

%CELLID_POS_call = ();
%CELLID1_CELLID2_check = ();

open OUT, ">$opt{'O'}";
foreach $cellID1 (keys %CELLID1_CELLID2_overlap) {
	foreach $cellID2 (keys %CELLID1_CELLID2_overlap) {
		if ($cellID1 ne $cellID2 && $CELLID1_CELLID2_overlap{$cellID1}{$cellID2} > 0) {
			if (defined $CELLID1_CELLID2_shared{$cellID1}{$cellID2}) {
				$shared = $CELLID1_CELLID2_shared{$cellID1}{$cellID2};
			} else {$shared = 0};
			$frac = sprintf("%.3f", $CELLID1_CELLID2_shared{$cellID1}{$cellID2}/$CELLID1_CELLID2_overlap{$cellID1}{$cellID2});
			print OUT "$cellID1\t$cellID2\t$CELLID1_CELLID2_overlap{$cellID1}{$cellID2}\t$shared\t$frac\n";
		}
	}
} close OUT;

exit;
