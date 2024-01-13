#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("1:O:q:m:w:n:", \%opt);

$qcut = "0.0001";
$mcut = "25";
$wsize = 1500;

$die = "

sciMET_collateDMRclusters.pl [options] -1 [DMR...pg_1].txt

Options:

-1   [STR]   DMR output for cluster 'pg_1' (req)
                Must have 'pg_1' as cluster and all other
                cluster files must have the same nomenclature
                Expects a header line
-O   [STR]   Output prefix (req)
                Adds: .qNN.m##.bed and .qNN_m##.merged.bed
-n   [STR]   Name of sample (included in output bed, def = pg_1 file prefix)
-q   [FLT]   Q-value cutoff (def = $qcut)
-m   [INT]   Percent methylation cut (def = $mcut)
-w   [INT]   Window size (def = $wsize)

";

if (!defined $opt{'1'} || !defined $opt{'O'}) {die $die};

if (defined $opt{'q'}) {$qcut = $opt{'q'}};
if (defined $opt{'m'}) {$mcut = $opt{'m'}};
if (defined $opt{'w'}) {$wsize = $opt{'w'}};

open OUT, "| sort -k1,1 -k 2,2n - > $opt{'O'}.q$qcut.m$mcut.bed";

($in_pfx,$in_sfx) = split(/pg_1/, $opt{'1'});

if (defined $opt{'n'}) {$name = $opt{'n'}} else {$name = $in_pfx};

$pad = int($wsize/2);

for ($i = 1; $i <= 100; $i++) {
	$file = "$in_pfx"."pg_$i"."$in_sfx";
	if (-e "$file") {
		print STDERR "\nLoading $file ... \n";
		open IN, "$file";
		$null = <IN>;
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			if ($P[6]<=$qcut&&abs($P[7])>=$mcut) {
				print OUT "$P[1]\t".($P[2]-$pad)."\t".($P[2]+$pad)."\t$name\tpg_$i\t$P[6]\t$P[7]\n";
			}
		} close IN;
	} else {
		print STDERR "\nFILE: $file not found. Ending cluster loading.\n";
		$i+=100;
	}
} close OUT;

system("cat $opt{'O'}.q$qcut.m$mcut.bed | mergeBed -c 5 -o distinct > $opt{'O'}.q$qcut.m$mcut.merged.bed");

exit;
