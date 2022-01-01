#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("L:O:m:Xt:C:H:bs", \%opt);

#defaults
$minSize = 10000000;
$threads = 1;

$die = "
Usage:

sciMET_extract.pl (options) [rmdup & filtered bam file]

Extracts methylation calls and produces a 'chroms' folder for
subsequent analyses. Will generate both a CG and CH folder.

Options:

-O   [STR]   Output prefix (def = bam prefix)
-L   [STR]   File of list of passing cellIDs (will filter to them)
-m   [INT]   Minimum chromosome size to retain (def = $minSize)
              Used to exclude random and other small contigs.
-t   [INT]   Threads for extraction (def = $threads)
-s           Sort the chrom files. This will take substantially
              longer to run. (def = no sort)

-b           Only run the bismark extract portion.

-C   [STR]   Do not run bismark extract for CG and instead
              provide the CG extracted file here.*
-H   [STR]   Do not run bismark extract for CH and instead
              provide the CH extracted file here.*

     *If providing -H and/or -C, a bam file is still
      required but only needs to have the correct header.
      If EITHER -H or -C is provided, then bismark extract
      will be skipped and only chrom extract and sorting
      will be performed with the provided file.
              
-X           Delete intermediate files (def = keep)

";

# PARSE OPTS

if (!defined $ARGV[0]) {die $die};
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.bam//;
}

if (defined $opt{'L'}) {
	open IN, "$opt{'L'}";
	while ($l = <IN>) {
		chomp $l;
		$l =~ s/\s.+$//;
		$PASSING_CELLS{$l} = 1;
	} close IN;
}

if (defined $opt{'m'}) {$minSize = $opt{'m'}};
if (defined $opt{'t'}) {$threads = $opt{'t'}};

open LOG, ">$opt{'O'}.extract.log";


# run bismark extract
if (!defined $opt{'C'} && !defined $opt{'H'}) {
	system("mkdir $opt{'O'}.bismark_extract");
	$extract_call = "bismark_methylation_extractor --multicore $threads --comprehensive --merge_non_CpG --single-end --no_header --gzip -o $opt{'O'}.bismark_extract $ARGV[0] >> $opt{'O'}.bismark_extract.log 2>> $opt{'O'}.bismark_extract.log";
	$ts = localtime(time);
	print LOG "$ts\tRunning bismark extraction:\n\t$extract_call\n";
	system($extract_call);
	if (defined $opt{'X'}) {system("rm -f $opt{'O'}.bismark_extract.log")};
}

if (defined $opt{'b'}) {
	exit;
}

# get chromosome names to prep for making files
$chromCT = 0;
open HEADER, "samtools view -H $ARGV[0] |";
while ($l = <HEADER>) {
	chomp $l;
	@P = split(/\t/, $l);
	if ($P[0] =~ /^\@SQ/) {
		$P[2] =~ s/^LN://;
		if ($P[2] >= $minSize) {
			$P[1] =~ s/^SN://;
			$CHROMS{$P[1]} = "$P[1]";
			$chromCT++;
		}
	}
} close HEADER;
$ts = localtime(time);
print LOG "$ts\tFound $chromCT chromosomes that met min size requirement.\n";


if (defined $opt{'H'} || (!defined $opt{'H'} && !defined $opt{'C'})) {
	$ts = localtime(time);
	print LOG "$ts\tBeginning CH chrom file generation.\n";

	$bases_processed = 0;
	$increment = 10000000; $report = $increment;

	$CONTEXT_total{'H'}=0;
	$CONTEXT_total{'h'}=0;
	$CONTEXT_total{'X'}=0;
	$CONTEXT_total{'x'}=0;

	system("mkdir $opt{'O'}.CH.chroms");
	foreach $chrom (keys %CHROMS) {
		$handle = $CHROMS{$chrom};
		if (!defined $opt{'s'}) {
			open $handle, "| gzip > $opt{'O'}.CH.chroms/$chrom.bed.gz";
		} else {
			open $handle, "| sort -k 1,1 -k2,2n | gzip > $opt{'O'}.CH.chroms/$chrom.bed.gz";
		}
	}

	if (defined $opt{'H'}) {
		$CH_file = $opt{'H'};
	} else {
		$extract_prefix = $ARGV[0];
		$extract_prefix =~ s/\.bam$//;
		$CH_file = "$opt{'O'}.bismark_extract/Non_CpG_context_$extract_prefix.txt.gz"
	}
	
	open IN, "zcat $CH_file |";
	while ($l = <IN>) {
		chomp $l;
		$bases_processed++;
		if ($bases_processed >= $report) {
			$report += $increment;
			$ts = localtime(time);
			print LOG "\t$ts\t$bases_processed bases processed. H=$CONTEXT_total{'H'},X=$CONTEXT_total{'X'},h=$CONTEXT_total{'h'},x=$CONTEXT_total{'x'}\n";
		}
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.+$//;
		if (!defined $opt{'L'} || defined $PASSING_CELLS{$cellID}) {
			if (defined $CHROMS{$P[2]}) {
				$handle = $CHROMS{$P[2]};
				print $handle "$P[2]\t$P[3]\t$P[3]\t$cellID\t$P[4]\n";
				
				if (!defined $CELLID_totalCHcov{$cellID}) {
					$CELLID_totalCHcov{$cellID} = 1;
					$cellCT++;
					$CELLID_CONTEXT_cov{$cellID}{'H'} = 0;
					$CELLID_CONTEXT_cov{$cellID}{'h'} = 0;
					$CELLID_CONTEXT_cov{$cellID}{'X'} = 0;
					$CELLID_CONTEXT_cov{$cellID}{'x'} = 0;
				} else {
					$CELLID_totalCHcov{$cellID}++;
				}
				$CELLID_CONTEXT_cov{$cellID}{$P[4]}++;
				$CONTEXT_total{$P[4]}++;
			}
		}
	} close IN;

	foreach $chrom (keys %CHROMS) {
		$handle = $CHROMS{$chrom};
		close $handle;
	}

	open VALS, ">$opt{'O'}.mCH.vals";
	open COV, ">$opt{'O'}.CH_cov.vals";
	foreach $cellID (keys %CELLID_totalCHcov) {
		$meth = sprintf("%.3f", (($CELLID_CONTEXT_cov{$cellID}{'H'}+$CELLID_CONTEXT_cov{$cellID}{'X'})/$CELLID_totalCHcov{$cellID})*100);
		print VALS "$cellID\t$meth\n";
		print COV "$cellID\t$CELLID_totalCHcov{$cellID}\n";
	} close VALS; close COV;

}


if (defined $opt{'C'} || (!defined $opt{'H'} && !defined $opt{'C'})) {

	$ts = localtime(time);
	print LOG "$ts\tBeginning CG chrom file generation.\n";

	$bases_processed = 0;
	$increment = 10000000; $report = $increment;

	$CONTEXT_total{'Z'}=0;
	$CONTEXT_total{'z'}=0;

	system("mkdir $opt{'O'}.CG.chroms");
	foreach $chrom (keys %CHROMS) {
		$handle = $CHROMS{$chrom};
		if (!defined $opt{'s'}) {
			open $handle, "| gzip > $opt{'O'}.CG.chroms/$chrom.bed.gz";
		} else {
			open $handle, "| sort -k 1,1 -k2,2n | gzip > $opt{'O'}.CG.chroms/$chrom.bed.gz";
		}
	}

	if (defined $opt{'C'}) {
		$CG_file = $opt{'C'};
	} else {
		$extract_prefix = $ARGV[0];
		$extract_prefix =~ s/\.bam$//;
		$CG_file = "$opt{'O'}.bismark_extract/CpG_context_$extract_prefix.txt.gz"
	}
	
	open IN, "zcat $CG_file |";
	while ($l = <IN>) {
		chomp $l;
		$bases_processed++;
		if ($bases_processed >= $report) {
			$report += $increment;
			$ts = localtime(time);
			print LOG "\t$ts\t$bases_processed bases processed. Z=$CONTEXT_total{'Z'},z=$CONTEXT_total{'z'}\n";
		}
		@P = split(/\t/, $l);
		$cellID = $P[0]; $cellID =~ s/:.+$//;
		if (!defined $opt{'L'} || defined $PASSING_CELLS{$cellID}) {
			if (defined $CHROMS{$P[2]}) {
				$handle = $CHROMS{$P[2]};
				print $handle "$P[2]\t$P[3]\t$P[3]\t$cellID\t$P[4]\n";
				
				if (!defined $CELLID_totalCGcov{$cellID}) {
					$CELLID_totalCGcov{$cellID} = 1;
					$cellCT++;
					$CELLID_CONTEXT_cov{$cellID}{'Z'} = 0;
					$CELLID_CONTEXT_cov{$cellID}{'z'} = 0;
				} else {
					$CELLID_totalCGcov{$cellID}++;
				}
				$CELLID_CONTEXT_cov{$cellID}{$P[4]}++;
				$CONTEXT_total{$P[4]}++;
			}
		}
	} close IN;

	foreach $chrom (keys %CHROMS) {
		$handle = $CHROMS{$chrom};
		close $handle;
	}

	open VALS, ">$opt{'O'}.mCG.vals";
	open COV, ">$opt{'O'}.CG_cov.vals";
	foreach $cellID (keys %CELLID_totalCGcov) {
		$meth = sprintf("%.3f", ($CELLID_CONTEXT_cov{$cellID}{'Z'}/$CELLID_totalCGcov{$cellID})*100);
		print VALS "$cellID\t$meth\n";
		print COV "$cellID\t$CELLID_totalCGcov{$cellID}\n";
	} close VALS; close COV;

}

if (defined $opt{'X'}) {system("rm -f -R $opt{'O'}.bismark_extract")};

exit;

