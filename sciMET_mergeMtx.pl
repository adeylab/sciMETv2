#!/usr/bin/perl

$die = "

sciMET_mergeMtx.pl [output name] [matrix 1] [matrix 2] (opt: [annot file out],[matrix1 annot],[matrix2 annot])

Must have no overlapping cells.
Only shared rows will be included in output.
Optional to output an annot file for the
   two input matrixes and cell IDs.

";

if (!defined $ARGV[2]) {die $die};

open M, "$ARGV[1]";
$head = <M>; chomp $head;
@M1_CELLS = split(/\t/, $head);
foreach $cellID (@M1_CELLS) {
	$M1_CELLIDS{$cellID} = 1;
}
while ($l = <M>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winName = shift(@P);
	@{$M1_WINDOW_vals{$winName}} = @P;
} close M;

$passing_rows = 0;
open M, "$ARGV[2]";
$head = <M>; chomp $head;
@M2_CELLS = split(/\t/, $head);
foreach $cellID (@M2_CELLS) {
	if (defined $M1_CELLIDS{$cellID}) {
		die "\nERROR: $cellID is present in BOTH matrix files!\n$die";
	}
}

open OUT, ">$ARGV[0]";
$out_header = "";
for ($i = 0; $i < @M1_CELLS; $i++) {$out_header .= "$M1_CELLS[$i]\t"};
for ($i = 0; $i < @M2_CELLS; $i++) {$out_header .= "$M2_CELLS[$i]\t"};
$out_header =~ s/\t$//;
print OUT "$out_header\n";

while ($l = <M>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winName = shift(@P);
	if (defined $M1_WINDOW_vals{$winName}) {
		print OUT "$winName";
		for ($i = 0; $i < @M1_CELLS; $i++) {
			print OUT "\t$M1_WINDOW_vals{$winName}[$i]";
		}
		for ($i = 0; $i < @M2_CELLS; $i++) {
			print OUT "\t$P[$i]";
		}
		print OUT "\n";
		$passing_rows++;
	}
} close M; close OUT;

print STDERR "INFO: $passing_rows rows present in both matrixes.\n";

if (defined $ARGV[3]) {
	($file,$annot1,$annot2) = split(/,/, $ARGV[3]);
	open OUT, ">$file";
	foreach $cellID (@M1_CELLS) {print OUT "$cellID\t$annot1\n"};
	foreach $cellID (@M2_CELLS) {print OUT "$cellID\t$annot2\n"};
	close OUT;
}

exit;