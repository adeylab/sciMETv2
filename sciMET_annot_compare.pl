#!/usr/bin/perl

$die = "

sciMET_annot_compare.pl [truth annot] [test annot] [outfile]

";

if (!defined $ARGV[2]) {die $die};

open IN, "$ARGV[0]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$trueAnnot{$annot} = $annot;
	$CELLID_trueAnnot{$cellID} = $annot;
} close IN;

open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	if (defined $CELLID_trueAnnot{$cellID}) {
		$TRUE_TEST_ct{$CELLID_trueAnnot{$cellID}}{$annot}++;
	}
	$testAnnot_ct{$annot}++;
} close IN;

open OUT, ">$ARGV[2]";
foreach $trueAnnot (keys %trueAnnot) {
	foreach $testAnnot (keys %testAnnot_ct) {
		if (defined $TRUE_TEST_ct{$trueAnnot}{$testAnnot}) {
			$count = $TRUE_TEST_ct{$trueAnnot}{$testAnnot};
			$pct = sprintf("%.2f", (($count/$testAnnot_ct{$testAnnot})*100));
		} else {
			$count = 0; $pct = "0.00";
		}
		print OUT "$trueAnnot\t$testAnnot\t$count\t$pct\n";
	}
} close OUT;
