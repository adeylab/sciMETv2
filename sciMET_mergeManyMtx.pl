#!/usr/bin/perl

$die = "

sciMET_mergeManyMtx.pl [output name] [matrix 1] [matrix 2] (matrix 3) ...

Must have no overlapping cells. Must be same rows.

Designed for merging the cells split into chroms
files during the extract process.

";

if (!defined $ARGV[2]) {die $die};

for ($i = 1; $i < @ARGV; $i++) {
	
}

open OUT, ">$ARGV[0]";
close OUT;