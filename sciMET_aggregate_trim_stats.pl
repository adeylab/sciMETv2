#!/usr/bin/perl

$die = "

sciMET_aggregate_trim_stats.pl [list of trim stats txt files, min 2] > [aggregated stats by cellID]

";

if (!defined $ARGV[1]) {die $die};

foreach $in (@ARGV) {
	open IN, "$in";
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$CELLID_total{$P[0]} += $P[1];
		$CELLID_retained{$P[0]} += $P[2];
	} close IN;
}

foreach $cellID (sort {$CELLID_total{$b}<=>$CELLID_total{$a}} keys %CELLID_total) {
	$pct = sprintf("%.2f", ($CELLID_retained{$cellID}/$CELLID_total{$cellID})*100);
	print "$cellID\t$CELLID_total{$cellID}\t$CELLID_retained{$cellID}\t$pct\n";
}

exit;
