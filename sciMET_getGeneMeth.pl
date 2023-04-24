#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("A:F:O:G:S:s:E:Vm:rL:xpC:a:", \%opt);

# defaults
$refGene = "/home/groups/ravnica/refs/hg38/hg38.refGene.txt";
$startPad = 20000;
$endPad = 20000;
$minCov = 50;

$die = "

Usage:

sciMET_getGeneMeth.pl -F [chromosome meth call folder] -O [output prefix / folder] [gene1] ... [chrN:#-#] ...

regions will be explicitly that window with no padding.
regions must be: chr:start-end

Required options:

-F   [STR]   Folder with chromosome-split methylation calls for cells
             generated using sciMET_extract.pl
			 Can be multiple folders comma-separated.
             Or can be a file ending in 'txt' with aline for each folder
-O   [STR]   Output prefix, creates or uses directory. Sugegsted to include CH or CG in name.

Default Options:

-A   [STR]   Annot file - will aggregate over clusters as well as cell-based
-a   [STR]   List of annotations to include, comma-separated. Excludes all other cells.
-G   [STR]   RefGene file (def = $refGene)
-S   [STR]   Bases from the start of TSS to begin window (def = $startPad, if -s suggested: 5000)
-s   [STR]   Bases from the end of TSS to end window (for CG; def = null, overrides -E, suggested: 2500)
-E   [STR]   Bases from TES to end window (for CH; def = $endPad)
-m   [INT]   Min number of bases in region to make a methylaiton call (def = $minCov)
-r           Replace genes in annot methylaiton file (def = keep previous)
-L   [STR]   Text file listing targets (can be supplemented with those listed in ARGS)
             (e.g. /home/groups/ravnica/refs/hg38/hg38.refGene.uniqNames.txt)
-C   [STR]   File listing CellIDs to include (def = all)
-x           Do not print out values for individual cells (requires -A)
-p           Plot after (def = np)

-V           Verbose

";

if (!defined $opt{'F'} ||
	!defined $opt{'O'}) {die $die};

if (defined $opt{'G'}) {$refGene = $opt{'G'}};

$opt{'F'} =~ s/\/$//;

if (defined $opt{'E'} && defined $opt{'s'}) {
	die "\nOptions -s and -E are mutually-exclusive!\n$die";
}

if (defined $opt{'E'}) {$endPad = $opt{'E'}};
if (defined $opt{'S'}) {$startPad = $opt{'S'}};
if (defined $opt{'s'}) {$startPad2 = $opt{'s'}};
if (defined $opt{'m'}) {$minCov = $opt{'m'}};
if (defined $opt{'x'} && !defined $opt{'A'}) {
	die "ERROR: To skip individual cell plotting, an annotaiton file MUST be provided!\n$die";
}

if (defined $opt{'C'}) {
	open IN, "$opt{'C'}";
	while ($l = <IN>) {
		chomp $l; $l =~ s/\s.+$//;
		$INCLUDE{$l} = 1;
	} close IN;
}

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

# pull targets to a bed:
open TAR, "| sort -k1,1 -k2,2n > $opt{'O'}.targets.bed";
$genes = 0; $regions = 0;
if (defined $opt{'L'}) {
	open LIST, "$opt{'L'}";
	while ($l = <LIST>) {
		chomp $l;
		$l =~ s/\s.+$//;
		push @ARGV, $l;
	}
	close LIST;
}
foreach $target (@ARGV) {
	if ($target =~ /:/) {
		($chr,$start,$end) = split(/[:-]/, $target);
		print TAR "$chr\t$start\t$end\t$chr\_$start\_$end\n";
		print STDERR "Region: $chr\t$start\t$end\t$chr\_$start\_$end\n";
		$regions++;
		$CHR_READS{$chr}++;
	} else {
		$genes++;
		$GENES{uc($target)} = 1;
		print STDERR "Gene: $target\n";
	}
}

print STDERR "Found $regions region targets and $genes gene targets.\n";

if ($genes > 0) {
	print STDERR "Opening refgene: $refGene\n";
	open GEN, "$refGene";
	# 585     NR_046018       chr1    +       11873   14409   14409   14409   3       11873,12612,13220,      12227,12721,14409,      0       DDX11L1 unk     unk     -1,-1,-1,
	#  0          1             2     3          4       5      6       7     8               9                       10              11        12     13     14         15
	while ($l = <GEN>) {
		chomp $l;
		@P = split(/\t/, $l);
		if (defined $GENES{uc($P[12])}) {
			$gene = uc($P[12]);
			if (!defined $FOUND{$gene}) {
				
				if (defined $opt{'s'}) {
					$start = $P[4]-$startPad;
					$end = $P[5]+$endPad;
				} else {
					if ($P[3] eq "+") {
						$start = $P[4]-$startPad;
						$end = $P[4]+$startPad2;
					} else {
						$start = $P[5]-$startPad2;
						$end = $P[5]+$startPad;
					}
				}
				
				if ($start<=0) {$start=1}; if ($end<=0) {$end=1};
				print TAR "$P[2]\t$start\t$end\t$gene\n";
				print STDERR "Gene region: $P[2]\t$start\t$end\t$gene\n";
				$CHR_READS{$P[2]}++;
				$FOUND{$gene} = 1;
			}
		}
	} close GEN;
}

close TAR;

#chr11   120821  120821  TCCAGAAGCCGTAGTAGTCCGTTCCAAT    X     chr start end target 
# 0        1       2                     3               4      5    6    7     8

if ($opt{'F'} =~ /txt$/) {
	@FOLDERS = ();
	open F, "$opt{'F'}";
	while ($l = <F>) {chomp $l; push @FOLDERS, $l};
	close F;
} else {
	@FOLDERS = split(/,/, $opt{'F'});
}
foreach $chr_dir (@FOLDERS) {
	foreach $chr (keys %CHR_READS) {
		if (-e "$chr_dir/$chr.bed.gz") {
			$chr_file = "$chr_dir/$chr.bed.gz";
			print STDERR "Running intersect: bedtools intersect -a $chr_file -b $opt{'O'}.targets.bed -wa -wb\n";
			open INT, "bedtools intersect -a $chr_file -b $opt{'O'}.targets.bed -wa -wb |";
			while ($l = <INT>) {
				chomp $l;
				@P = split(/\t/, $l);
				if (defined $opt{'V'}) {print STDERR "Intersect line: $l\n"};

				if (!defined $opt{'C'} || defined $INCLUDE{$P[8]}) {
					if (!defined $opt{'a'} || defined $ANNOT_include{$CELLID_annot{$P[3]}}) {
						$TARGET_CELLID_cov{$P[8]}{$P[3]}++;
						if ($P[4] =~ /[XHZ]/) {
							$TARGET_CELLID_meth{$P[8]}{$P[3]}++;
						}
						
						if (defined $CELLID_annot{$P[3]}) {
							$annot = $CELLID_annot{$P[3]};
							if (!defined $opt{'a'} || defined $ANNOT_include{$annot}) {
								$TARGET_ANNOT_cov{$P[8]}{$annot}++;
								if ($P[4] =~ /[XHZ]/) {
									$TARGET_ANNOT_meth{$P[8]}{$annot}++;
								}
							}
					}
					}
				}
			} close INT;
		} else {
			print STDERR "WARNING: Cannot find $chr.bed.gz in target directory $chr_dir! Skipping!\n";
		}
	}
}

# calc and print
if (!defined $opt{'x'}) {
	if (-d "$opt{'O'}") {
		print STDERR "$opt{'O'} directory exists! Will add to directory.\n";
	} else {
		print STDERR "Creating $opt{'O'} directory.\n";
		system("mkdir $opt{'O'}");
	}

	foreach $target (keys %TARGET_CELLID_cov) {
		open OUT, ">$opt{'O'}/$target.vals";
		foreach $cellID (keys %{$TARGET_CELLID_cov{$target}}) {
			if ($TARGET_CELLID_cov{$target}{$cellID} >= $minCov) {
				if (!defined $TARGET_CELLID_meth{$target}{$cellID}) {$TARGET_CELLID_meth{$target}{$cellID} = 0};
				$meth = sprintf("%.2f", ($TARGET_CELLID_meth{$target}{$cellID}/$TARGET_CELLID_cov{$target}{$cellID})*100);
				print OUT "$cellID\t$meth\n";
			}
		}
		close OUT;
	}
}

if (defined $opt{'A'}) {
	if (-e "$opt{'O'}.annot.meth.txt") {
		print STDERR "$opt{'O'}.annot.meth.txt exists! Looking for existing targets.\n";
		open CURR, "$opt{'O'}.annot.meth.txt";
		$CURR_OUT = ""; $new_file = 0;
		while ($l = <CURR>) {
			chomp $l;
			($target,$annot,$meth,$zscore) = split(/\t/, $l);
			if (!defined $ANNOT_ct{$annot}) {
				print STDERR "WARNING: Annot $annot is not in provided annot file. WIll create a new output file.\n";
				$new_file = 1;
			}
			if (defined $TARGET_ANNOT_cov{$target}) {
				if (!defined $opt{'r'}) {
					print STDERR "\tTarget $target exists in current file - original will be kept (unless annot mismatch detected)\n";
					$CURR_OUT .= "$l\n";
				} else {
					print STDERR "\tTarget $target exists in current file and will be rewritten with new results.\n";
				}
			}
		} close CURR;
		
		if ($new_file>0) {
			open OUT, ">$opt{'O'}.annot.meth.txt-ANNOT_CONFLICT\n";
		} else {
			open OUT, ">$opt{'O'}.annot.meth.txt";
			print OUT "$CURR_OUT";
		}
		
	} else {
		open OUT, ">$opt{'O'}.annot.meth.txt";
	}
	
	
	foreach $target (keys %TARGET_ANNOT_cov) {
		%ANNOT_meth = (); @METH = (); $mean = 0; $stdev = 0; $sum = 0; $stdev_sum = 0; $passing = 0;
		foreach $annot (keys %{$TARGET_ANNOT_cov{$target}}) {
			if ($TARGET_ANNOT_cov{$target}{$annot} >= $minCov) {
				if (!defined $TARGET_ANNOT_meth{$target}{$annot}) {$TARGET_ANNOT_meth{$target}{$annot} = 0};
				$ANNOT_meth{$annot} = sprintf("%.2f", ($TARGET_ANNOT_meth{$target}{$annot}/$TARGET_ANNOT_cov{$target}{$annot})*100);
				push @METH, $ANNOT_meth{$annot};
				$sum += $ANNOT_meth{$annot};
				$passing++;
			}
		}
		if ($passing >= 3) {
			$mean = $sum/$passing;
			foreach $val (@METH) {$stdev_sum+=abs($val-$mean)**2};
			$stdev = sqrt($stdev_sum/($passing-1));
			if ($stdev>0) {
				foreach $annot (keys %{$TARGET_ANNOT_cov{$target}}) {
					if ($TARGET_ANNOT_cov{$target}{$annot} >= $minCov) {
						$zscore = ($ANNOT_meth{$annot}-$mean)/$stdev;
						print OUT "$target\t$annot\t$ANNOT_meth{$annot}\t$zscore\n";
					}
				}
			}
		}
	} close OUT;
	if (defined $opt{'p'}) {system("plotGeneByAnnot.pl $opt{'O'}.annot.meth.txt")};
}
