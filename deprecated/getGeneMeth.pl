#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("A:F:O:G:S:s:E:Vm:", \%opt);

# defaults
$refGene = "/home/groups/oroaklab/refs/hg38/hg38.refGene.txt";
$startPad = 20000;
$endPad = 20000;
$minCov = 50;

$die = "

Usage:

getGeneMeth.pl -F [chromosome meth call folder] -O [output prefix] [gene1] ... [chrN:#-#] ...

regions will be explicitly that window with no padding.
regions must be: chr:start-end

Required options:

-F   [STR]   Folder with chromosome-split methylation calls for cells
             generated using CHmeth2mtx/CGmeth2mtx
-O   [STR]   Output prefix

Default Options:

-A   [STR]   Annot file - will aggregate over clusters as well as cell-based
-G   [STR]   RefGene file (def = $refGene)
-S   [STR]   Bases from the start of TSS to begin window (def = $startPad)
-s   [STR]   Bases from the end of TSS to end window (for CG; def = null, overrides -E)
-E   [STR]   Bases from TES to end window (for CH; def = $endPad)
-m   [INT]   Min number of bases in region to make a methylaiton call (def = $minCov)

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

# pull targets to a bed:
open TAR, ">$opt{'O'}.targets.bed";
$genes = 0; $regions = 0;
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

$chr_dir = $opt{'F'};
foreach $chr (keys %CHR_READS) {
	if (-e "$chr_dir/$chr.CH.bed.gz") {
		$chr_file = "$chr_dir/$chr.CH.bed.gz";
	} elsif (-e "$chr_dir/$chr.CG.bed.gz") {
		$chr_file = "$chr_dir/$chr.CG.bed.gz";
	} else {
		print STDERR "FATAL ERROR: Cannot find $chr.CH.bed.gz OR $chr.CG.bed.gz in target directory $chr_dir!\n";
		die;
	}
	print STDERR "Running intersect: bedtools intersect -a $chr_file -b $opt{'O'}.targets.bed -wa -wb\n";
	open INT, "bedtools intersect -a $chr_file -b $opt{'O'}.targets.bed -wa -wb |";
	while ($l = <INT>) {
		chomp $l;
		@P = split(/\t/, $l);
		if (defined $opt{'V'}) {print STDERR "Intersect line: $l\n"};

		$TARGET_CELLID_cov{$P[8]}{$P[3]}++;
		if ($P[4] =~ /[XHZ]/) {
			$TARGET_CELLID_meth{$P[8]}{$P[3]}++;
		}
		
		if (defined $CELLID_annot{$P[3]}) {
			$annot = $CELLID_annot{$P[3]};
			$TARGET_ANNOT_cov{$P[8]}{$annot}++;
			if ($P[4] =~ /[XHZ]/) {
				$TARGET_ANNOT_meth{$P[8]}{$annot}++;
			}
		}
	} close INT;
}

# calc and print

system("mkdir $opt{'O'}.target_cell_values");
foreach $target (keys %TARGET_CELLID_cov) {
	open OUT, ">$opt{'O'}.target_cell_values/$target.vals";
	foreach $cellID (keys %{$TARGET_CELLID_cov{$target}}) {
		if ($TARGET_CELLID_cov{$target}{$cellID} >= $minCov) {
			if (!defined $TARGET_CELLID_meth{$target}{$cellID}) {$TARGET_CELLID_meth{$target}{$cellID} = 0};
			$meth = sprintf("%.2f", ($TARGET_CELLID_meth{$target}{$cellID}/$TARGET_CELLID_cov{$target}{$cellID})*100);
			print OUT "$cellID\t$meth\n";
		}
	}
	close OUT;
}

if (defined $opt{'A'}) {
	open OUT, ">$opt{'O'}.annot.meth.txt";
	foreach $target (keys %TARGET_ANNOT_cov) {
		foreach $annot (keys %{$TARGET_ANNOT_cov{$target}}) {
			if ($TARGET_ANNOT_cov{$target}{$annot} >= $minCov) {
				if (!defined $TARGET_ANNOT_meth{$target}{$annot}) {$TARGET_ANNOT_meth{$target}{$annot} = 0};
				$meth = sprintf("%.2f", ($TARGET_ANNOT_meth{$target}{$annot}/$TARGET_ANNOT_cov{$target}{$annot})*100);
				print OUT "$target\t$annot\t$meth\n";
			}
		}
	} close OUT;
}