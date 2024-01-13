#!/usr/bin/perl

$die = "

sciMET_BSBolt_StrandBias.pl [myBam] [out prefix] (optional: CellID list)

Generates:

[out].cellStrandBias.txt (per cell)
[out].strandBias.txt (global)

If a text file of cellIDs is provided it will only include information
from those cell IDs in both global and cell-level outputs.

";

if (!defined $ARGV[1]) {die $die};
%GLOBAL = ( 'W_C2T' => 0,
			'W' => 0,
			'C' => 0,
			'W_G2A' => 0,
			'C_C2T' => 0,
			'C_G2A' => 0,
			'C2T' => 0,
			'G2A' => 0);
			
if (defined $ARGV[2]) {
	open IN, "$ARGV[2]";
	while ($l = <IN>) {
		chomp $l;
		$l =~ s/\s.+$//;
		$CELLID_FILT{$l} = 1;
	} close IN;
}

open IN, "samtools view $ARGV[0] |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$cellID = $P[0]; $cellID =~ s/:.+$//;
	if (!defined $ARGV[2] || defined $CELLID_FILT{$cellID}) {
		for ($i = 11; $i < @P; $i++) {
			if ($P[$i] =~ /^YS/) {
				$strand = $P[$i];
				$strand =~ s/YS:Z://;
			}
		}
		if (defined $GLOBAL{$strand}) {
			$GLOBAL{$strand}++;
			$total++;
			$CELLID_STRAND{$cellID}{$strand}++;
			$CELLID_sum{$cellID}++;
			if ($strand =~ /^W/) {$GLOBAL{'W'}++; $CELLID_STRAND{$cellID}{'W'}++};
			if ($strand =~ /^C/) {$GLOBAL{'C'}++; $CELLID_STRAND{$cellID}{'C'}++};
			if ($strand =~ /C2T/) {$GLOBAL{'C2T'}++; $CELLID_STRAND{$cellID}{'C2T'}++};
			if ($strand =~ /G2A/) {$GLOBAL{'G2A'}++; $CELLID_STRAND{$cellID}{'G2A'}++};
		}
	}
} close IN;

open OUT, ">$ARGV[1].strandBias.txt";
open OUT2, ">$ARGV[1].cellStrandBias.txt";
print OUT2 "#CellID";
foreach $strand (sort keys %GLOBAL) {
	$pct = sprintf("%.2f", ($GLOBAL{$strand}/$total)*100);
	print OUT "$strand\t$GLOBAL{$strand}\t$pct\n";
	print OUT2 "\t$strand\_count\t$strand\_pct";
}
close OUT;
print OUT2 "\n";

foreach $cellID (keys %CELLID_STRAND) {
	print OUT2 "$cellID";
	if ($CELLID_sum{$cellID} > 0) {
		foreach $strand (keys %GLOBAL) {
			$pct = sprintf("%.2f", ($CELLID_STRAND{$cellID}{$strand}/$CELLID_sum{$cellID})*100);
			print OUT2 "\t$CELLID_STRAND{$cellID}{$strand}\t$pct";
		}
	}
	print OUT2 "\n";
}
close OUT2;
