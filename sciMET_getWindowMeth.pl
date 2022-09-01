#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("F:B:A:a:O:m:", \%opt);

# defaults
$minCov=10;

$die = "

Usage:

sciMET_getWindowMeth.pl -F [chromosome meth call folder(s)] -O [output prefix / folder] [gene1] -B [Bed of targets]

Required options:

-F   [STR]   Folder with chromosome-split methylation calls for cells
             generated using sciMET_extract.pl
			 Can be multiple folders comma-separated.
             Or can be a single text file (end in .txt) with list of folders
-O   [STR]   Output prefix, creates or uses directory.
-B   [STR]   Bed file with targets to profile.
             Generate using: sciMET_featuresToScanBed.pl
             Each feature needs 4x columns - chrom, start, end, windowPosition
-m   [INT]   Min coverage to report call, annotation-based if -A (def = 10)

Other Options:
-A   [STR]   Annot file
-a   [STR]   Comma separated list of annotations to include (subset to)


";

if (!defined $opt{'F'} ||
	!defined $opt{'O'} ||
	!defined $opt{'B'}) {die $die};

if (defined $opt{'m'}) {$minCov = $opt{'m'}};

$opt{'F'} =~ s/\/$//;

# read annot file
if (defined $opt{'A'}) {
	open IN, "$opt{'A'}";
	while ($l = <IN>) {
		chomp $l;
		($cellID,$annot) = split(/\t/, $l);
		$CELLID_annot{$cellID} = $annot;
		$ANNOT_ct{$annot}++;
	} close IN;
}

if (defined $opt{'a'}) {
	if (!defined $opt{'A'}) {die "\nERROR: -a annotaiton list requires an annotation file!\n$die"};
	@ANNOT_LIST = split(/,/, $opt{'a'});
	foreach $annot (@ANNOT_LIST) {
		$ANNOT_include{$annot} = 1;
	}
}

for ($i = 1; $i <= 22; $i++) {
	$CHRS{"chr$i"} = 1;
} $CHRS{"chrX"} = 1; $CHRS{"chrY"} = 1;

if ($opt{'F'} =~ /txt$/) {
	@FOLDERS = ();
	open F, "$opt{'F'}";
	while ($l = <F>) {chomp $l; push @FOLDERS, $l};
	close F;
} else {
	@FOLDERS = split(/,/, $opt{'F'});
}
foreach $chr_dir (@FOLDERS) {
	foreach $chr (keys %CHRS) {
		if (-e "$chr_dir/$chr.bed.gz") {
			$chr_file = "$chr_dir/$chr.bed.gz";
			print STDERR "Running intersect: bedtools intersect -a $chr_file -b $opt{'B'} -wa -wb\n";
			open INT, "bedtools intersect -a $chr_file -b $opt{'B'} -wa -wb |";
			
			# chr pos pos cellID methID chr start end winID
			#  0   1   2     3      4    5    6    7    8
			
			while ($l = <INT>) {
				chomp $l;
				@P = split(/\t/, $l);
				if (defined $opt{'A'}) {
					$annot = $CELLID_annot{$P[3]};
					if (!defined $opt{'a'} || defined $ANNOT_include{$annot}) {
						$WINID_ANNOT_BASE_count{$P[8]}{$annot}{$P[4]}++;
					}
				}
				$WINID_BASE_count{$P[8]}{$P[4]}++;
				
				
			} close INT;
		} else {
			print STDERR "WARNING: Cannot find $chr.bed.gz in target directory $chr_dir! Skipping!\n";
		}
	}
}


open OUT, ">$opt{'O'}.CG.txt";
foreach $winID (sort {$a<=>$b} keys %WINID_BASE_count) {
	if (($WINID_BASE_count{$winID}{'Z'}+$WINID_BASE_count{$winID}{'z'})>=$minCov) {
		$meth = sprintf("%.3f", ($WINID_BASE_count{$winID}{'Z'}/($WINID_BASE_count{$winID}{'Z'}+$WINID_BASE_count{$winID}{'z'}))*100);
		$cov = ($WINID_BASE_count{$winID}{'Z'}+$WINID_BASE_count{$winID}{'z'});
		print OUT "$winID\tALL\t$meth\t$cov\n";
		if (defined $opt{'A'}) {
			foreach $annot (keys %{$WINID_ANNOT_BASE_count{$winID}}) {
				if ($WINID_ANNOT_BASE_count{$winID}{$annot}{'z'}>=$minCov) {
					$meth = sprintf("%.3f", ($WINID_ANNOT_BASE_count{$winID}{$annot}{'Z'}/($WINID_ANNOT_BASE_count{$winID}{$annot}{'Z'}+$WINID_ANNOT_BASE_count{$winID}{$annot}{'z'}))*100);
					$cov = ($WINID_ANNOT_BASE_count{$winID}{$annot}{'Z'}+$WINID_ANNOT_BASE_count{$winID}{$annot}{'z'});
					print OUT "$winID\t$annot\t$meth\t$cov\n";
				}
			}
		}
	}
} close OUT;



open OUT, ">$opt{'O'}.CH.txt";
foreach $winID (sort {$a<=>$b} keys %WINID_BASE_count) {
	$cov = ($WINID_BASE_count{$winID}{'H'}+$WINID_BASE_count{$winID}{'X'}+$WINID_BASE_count{$winID}{'h'}+$WINID_BASE_count{$winID}{'x'});
	if ($cov >= $minCov) {
		$meth = sprintf("%.3f", (($WINID_BASE_count{$winID}{'H'}+$WINID_BASE_count{$winID}{'X'})/$cov)*100);
		print OUT "$winID\tALL\t$meth\t$cov\n";
		if (defined $opt{'A'}) {
			foreach $annot (keys %{$WINID_ANNOT_BASE_count{$winID}}) {
				if ($WINID_ANNOT_BASE_count{$winID}{$annot}{'x'}>0 || $WINID_ANNOT_BASE_count{$winID}{$annot}{'h'}>0) {
					$cov = ($WINID_ANNOT_BASE_count{$winID}{$annot}{'X'}+$WINID_ANNOT_BASE_count{$winID}{$annot}{'H'}+$WINID_ANNOT_BASE_count{$winID}{$annot}{'x'}+$WINID_ANNOT_BASE_count{$winID}{$annot}{'h'});
					if ($cov >= $minCov) {
						$meth = sprintf("%.3f", (($WINID_ANNOT_BASE_count{$winID}{$annot}{'X'}+$WINID_ANNOT_BASE_count{$winID}{$annot}{'H'})/$cov)*100);
						print OUT "$winID\t$annot\t$meth\t$cov\n";
					}
				}
			}
		}
	}
} close OUT;
