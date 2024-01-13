#!/usr/bin/perl

$die = "

sciMET_mergeSVD.pl [output matrix] [input matrix ] [input matrix 2] (input matrix N) ...

Merges matrixes, retaining only the union of cells and all rows from
all matrixes. Can be used to merge SVD matrixes from CH and CG contexts.
Appends fileID to the end of the windowID.

";

if (!defined $ARGV[2]) {die $die};
$fileCT = 0;
for ($i = 1; $i < @ARGV; $i++) {
	$fileCT++;
	open IN, "$ARGV[$i]";
	$head = <IN>;
	chomp $head;
	@CELLS = split(/\t/, $head);
	foreach $cellID (@CELLS) {$CELLID_fileCT{$cellID}++};
	while ($l = <IN>) {
		chomp $l;
		@P = split(/\t/, $l);
		$rowID = shift(@P);
		$rowID = $rowID."_$i";
		for ($j = 0; $j < @P; $j++) {
			$ROWID_CELLID_val{$rowID}{$CELLS[$j]} = $P[$j];
		}
	} close IN;
}

@CELLS = (); $out_head = "";
foreach $cellID (keys %CELLID_fileCT) {
	$cellCT++;
	if ($CELLID_fileCT{$cellID} == $fileCT) {
		$includedCells++;
		push @CELLS, $cellID;
		$out_head .= "$cellID\t";
	}
} $out_head =~ s/\t$//;

open OUT, ">$ARGV[0]";
print OUT "$out_head\n";
foreach $rowID (keys %ROWID_CELLID_val) {
	print OUT "$rowID";
	for ($i = 0; $i < @CELLS; $i++) {
		$cellID = $CELLS[$i];
		print OUT "\t$ROWID_CELLID_val{$rowID}{$cellID}";
	} print OUT "\n";
} close OUT;