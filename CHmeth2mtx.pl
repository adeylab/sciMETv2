#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("f:c:F:L:O:B:", \%opt);

#defaults
$minCalled = 100;
$qualFrac = 0.9;

$die = "
Usage:

CHmeth2mtx.pl [options] -B [window bed file] -O [output prefix] -F [methylation chrom folder] ...

Options:

-F   [DIR]   Directory with CH methylation chromosome windows (req)
-B   [STR]   Bed file of windows for computing methylaiton (req)
-O   [STR]   Output prefix (req)

-c   [INT]   Min bases called per cell per window to qualify (def = $minCalled)
-f   [STR]   Fraction of cells that qualify to include window  (def = $qualFrac)

-L   [STR]   File of list of passing cellIDs (will filter to them)

";

# PARSE OPTS

if (!defined $opt{'F'} ||
	!defined $opt{'O'} ||
	!defined $opt{'B'}) {die $die};


if (defined $opt{'L'}) {
	open IN, "$opt{'L'}";
	while ($l = <IN>) {
		chomp $l;
		$l =~ s/\s.+$//;
		$PASSING_CELLS{$l} = 1;
	} close IN;
}

if (defined $opt{'c'}) {$minCalled = $opt{'c'}};
if (defined $opt{'f'}) {$qualFrac = $opt{'f'}};

$winID = 0;
open BED, "$opt{'B'}";
while ($l = <BED>) {
	chomp $l;
	@P = split(/\t/, $l);
	$win_name = "$P[0]_$P[1]_$P[2]";
	$winID++;
	$WINID_name{$winID} = $win_name;
	$WINNAME_winID{$win_name} = $winID;
	$CHR_START_win{$P[0]}{$P[1]} = $winID;
	$WINID_end{$winID} = $P[2];

	$CHR_HANDLES{$P[0]} = 1;
	
}
$lastWin = $winID;

$CONTEXT_total{'H'}=0;
$CONTEXT_total{'h'}=0;
$CONTEXT_total{'X'}=0;
$CONTEXT_total{'x'}=0;
$not_in_window=0;
$win_with_cov=0;
@CELLIDs = ();

foreach $chr (keys %CHR_HANDLES) {
	$chr_dir = $opt{'F'};
	
	if (-e "$chr_dir/$chr.CH.bed.gz") {
		$ts = localtime(time);
		print STDERR "$ts\tFound $chr_dir/$chr.CH.bed.gz, processing...\n";
		
		open INT, "bedtools intersect -a $chr_dir/$chr.CH.bed.gz -b $ARGV[0] -wa -wb |";
		# 0   1    2     3      4   5    6    7
		#chr pos1 pos2 cellID meth chr start end ...
		
		while ($l = <INT>) {
			chomp $l;
			@P = split(/\t/, $l);
			$cellID = $P[3];
			if (!defined $opt{'L'} || defined $PASSING_CELLS{$cellID}) {
				$winName = "$P[5]_$P[6]_$P[7]";
				$winID = $WINNAME_winID{$winName};
				if (!defined $WINID_totalCov{$winID}) {
					$win_with_cov++;
				}
				$WINID_totalCov{$winID}++;
				
				if (!defined $WINID_CELLID_covstr{$winID}{$cellID}) {
					$WINID_cellCT{$winID}++;
				}
				$WINID_totalCov{$winID}++;
				$WINID_CELLID_covstr{$winID}{$cellID} .= $P[4];
				
				if (defined $opt{'F'}) {
					if (!defined $CELLID_totalCov{$cellID}) {
						$CELLID_totalCov{$cellID} = 1;
						$cellCT++;
						push @CELLIDs, $cellID;
						$CELLID_CONTEXT_cov{$cellID}{'H'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'h'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'X'} = 0;
						$CELLID_CONTEXT_cov{$cellID}{'x'} = 0;
					} else {
						$CELLID_totalCov{$cellID}++;
					}
					$CELLID_CONTEXT_cov{$cellID}{$P[4]}++;
					$CONTEXT_total{$P[4]}++;
				}
			}
			
		} close INT;
		
	} else {
		print STDERR "WARNING: Cannot find $chr_dir/$chr.CH.bed.gz!\n";
	}
}

open CHM, ">$ARGV[1].CH.mtx";
open CHMn, ">$ARGV[1].CH.ratio.mtx";

$header = "";
for ($i = 0; $i < @CELLIDs; $i++) {
	$cellID = $CELLIDs[$i];
	$header .= "$cellID\t";
	$CELLID_methTotal{$cellID} = $CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'};
	if (($CELLID_CONTEXT_cov{$cellID}{'h'}+$CELLID_CONTEXT_cov{$cellID}{'x'}+$CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})<1) {$CELLID_CONTEXT_cov{$cellID}{'h'}=100};
		# pseudocount if 0; shouldn't ever happen, but will cancel cell out in final matrix (everything will be mean, but cell will get filtered later)
	if (($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})<1) {$CELLID_CONTEXT_cov{$cellID}{'H'}=1}; # pseudocount if 0
	$CELLID_mCH{$cellID} = sprintf("%.3f", ($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})/($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'}+$CELLID_CONTEXT_cov{$cellID}{'h'}+$CELLID_CONTEXT_cov{$cellID}{'x'}));
	$CELLID_qualWin{$cellID} = 0;
}
$header =~ s/\t$//;

print CHM "$header\n";
print CHMn "$header\n";

for ($winID = 1; $winID <= $lastWin; $winID++) {
	$qualCells = 0;
	%CELLID_winMeth = ();
	foreach $cellID (@CELLIDs) {
		$cov = length($WINID_CELLID_covstr{$winID}{$cellID});
		if ($cov>=$minCalled) {
			$qualCells++;
			$CELLID_qualWin{$cellID}++;
		}
		@M = split(//, $WINID_CELLID_covstr{$winID}{$cellID});
		$meth = 0;
		foreach $base (@M) {
			if ($base eq "H" || $base eq "X") {$meth++};
		}
		if ($cov > 0) {$frac = sprintf("%.3f", $meth/$cov)} else {$frac = $CELLID_mCH{$cellID}};
		$CELLID_winMeth{$cellID} = $frac;
	}
	if (($qualCells/$cellCT)>=$qualFrac) {
		print CHM "$WINID_name{$winID}";
		print CHMn "$WINID_name{$winID}";
		for ($i = 0; $i < @CELLIDs; $i++) {
			print CHM "\t$CELLID_winMeth{$CELLIDs[$i]}";
			if ($CELLID_mCH{$CELLIDs[$i]} == 0) {
				$ratio = 1;
			} else {
				$ratio = sprintf("%.3f", $CELLID_winMeth{$CELLIDs[$i]}/$CELLID_mCH{$CELLIDs[$i]});
			}
			print CHMn "\t$ratio";
		}
		print CHM "\n";
		print CHMn "\n";
	} else {
		$failWinCT++;
	}
} close CHM; close CHMn;

open CELLS, ">$ARGV[1].cellStats.txt";
foreach $cellID (@CELLIDs) {
	print CELLS "$cellID\t$CELLID_totalCov{$cellID}\t$CELLID_methTotal{$cellID}\t$CELLID_mCH{$cellID}\t$CELLID_qualWin{$cellID}\n";
} close CELLS;

open LOG, ">$ARGV[1].otherStats.txt";
print LOG "Bed file: $ARGV[0]
Meth Files:";
for ($i = 2; $i < @ARGV; $i++) {print LOG " $ARGV[$i]"};
$passWin = sprintf("%.2f", (($lastWin-$failWinCT)/$lastWin)*100);
$passWinCT = ($lastWin-$failWinCT);
print LOG "
Min called per cell per window: $minCalled
Min qual cell frac for window: $qualFrac

Total Cells: $cellCT
Total Windows: $lastWin
Failing Windows: $failWinCT
Passing Windows: $passWinCT
Passing Window Percent: $passWin\n";
close LOG;