#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("C:W:M:O:w:L:F:", \%opt);

#defaults
$minCells = 0.75;
$minCov = 20;
$minFrac = 0.75;

$die = "

sciMET_filtMtx.pl (options) -C [coverage matrix] -M [methylation or ratio matrix] -O [output prefix]

Options:

-C   [STR]   Coverage matrix (req)
-M   [STR]   Methylation or ratio matrix (to be filtered, req)
-O   [STR]   Output prefix (adds a stats and mtx file, req)

-w   [INT]   Min coverage of cell in window to pass (def = $minCov)
-W   [FLT]   Min fraction of passing windows for cell to pass (def = $minFrac)
-F   [FLT]   Min fraction of passing cells meeting min coverage in window to include (def = $minCells)
-L   [STR]   List of CellIDs to include (does not factor into other calc)

";

if (!defined $opt{'O'} ||
	!defined $opt{'C'} ||
	!defined $opt{'M'}) {
		die $die;
}

if (defined $opt{'w'}) {$minCov = $opt{'w'}};
if (defined $opt{'W'}) {$minFrac = $opt{'W'}};
if (defined $opt{'F'}) {$minCells = $opt{'F'}};

if (defined $opt{'L'}) {
	open CELLS, "$opt{'L'}";
	while ($l = <CELLS>) {
		chomp $l;
		$INCLUDE_LIST{$l} = 1;
	} close CELLS;
}

%CELLID_passwin = ();
$totalWin = 0;
open COV, "$opt{'C'}";
$head = <COV>; chomp $head;
@CELLS = split(/\t/, $head);
while ($l = <COV>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winID = shift(@P);
	$totalWin++;
	for ($i = 0; $i < @P; $i++) {
		if (!defined $opt{'L'} || defined $INCLUDE_LIST{$CELLS[$i]}) {
			if ($P[$i]>=$minCov) {$CELLID_passwin{$CELLS[$i]}++};
		}
	}
} close COV;

%CELLID_pass = ();
$out_header = ""; $pass_cell_ct = 0;
foreach ($i = 0; $i < @CELLS; $i++) {
	$cellID = $CELLS[$i];
	if (!defined $opt{'L'} || defined $INCLUDE_LIST{$CELLS[$i]}) {
		if (($CELLID_passwin{$cellID}/$totalWin)>=$minFrac) {
			$CELLID_pass{$cellID} = 1;
			$out_header .= "$cellID\t";
			$pass_cell_ct++;
		}
	}
} $out_header =~ s/\t$//;

open COV, "$opt{'C'}";
$pass_win_ct = 0;
$null = <COV>;
while ($l = <COV>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winID = shift(@P);
	$pass = 0;
	for ($i = 0; $i < @P; $i++) {
		if (!defined $opt{'L'} || defined $INCLUDE_LIST{$CELLS[$i]}) {
			if (defined $CELLID_pass{$CELLS[$i]}) {
				if ($P[$i] >= $minCov) {
					$pass++;
				}
			}
		}
	}
	if (($pass/$pass_cell_ct)>=$minCells) {
		$WINID_pass{$winID} = 1;
		$pass_win_ct++;
	}
} close COV;

open OUT, ">$opt{'O'}.mtx";
open IN, "$opt{'M'}";
print OUT "$out_header\n";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$winID = shift(@P);
	if (defined $WINID_pass{$winID}) {
		print OUT "$winID";
		for ($i = 0; $i < @P; $i++) {
			if (defined $CELLID_pass{$CELLS[$i]}) {
				print OUT "\t$P[$i]";
			}
		}
		print OUT "\n";
	}
} close OUT; close IN;

open OUT, ">$opt{'O'}.cellIDs.txt";
foreach ($i = 0; $i < @CELLS; $i++) {
	if (defined $CELLID_pass{$CELLS[$i]}) {
		print OUT "$CELLS[$i]\n";
	}
} close OUT;

open OUT, ">$opt{'O'}.stats.txt";
$totalCells = @CELLS;
$pass_cell_pct = sprintf("%.2f", ($pass_cell_ct/$totalCells)*100);
$pass_win_pct = sprintf("%.2f", ($pass_win_ct/$totalWin)*100);
print OUT "sciMET_filtMtx.pl Stats:

-O $opt{'O'}
-C $opt{'C'}
-M $opt{'M'}

Min cov per cell epr win (-w) = $minCov
Min fraction of win per cell (-W) = $minFrac
Min frac of pass cells with pass cov (-F) = $minCells

Total cells = $totalCells
Passing cells = $pass_cell_ct ($pass_cell_pct)
Total windows = $totalWin
Passing windows = $pass_win_ct ($pass_win_pct)
";
close OUT;