#!/usr/bin/perl

use Getopt::Std; %opt = ();
use List::Util qw(shuffle);

getopts("L:O:m:Xt:C:H:bN:G:", \%opt);

#defaults
$minSize = 10000000;
$threads = 1;
$mCH_max = 10;
$cells_per_folder = 100;

$die = "
Usage:

sciMET_extract.pl (options) [rmdup & filtered bam file]

Extracts methylation calls and produces a 'chroms' folder for
subsequent analyses. Will generate both a CG and CH folder.

Options:

-O   [STR]   Output prefix (def = bam prefix)
-L   [STR]   File of list of passing cellIDs (will filter to them)
              Also used to pre-set cells to folders which will make
              them more uniform in size. (recommended)
                       !!!! NOT YET FULLY IMPLEMENTED
-m   [INT]   Minimum chromosome size to retain (def = $minSize)
              Used to exclude random and other small contigs.
-t   [INT]   Threads for extraction (def = $threads)

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
	  
-G   [BED]   Bed file of GpCs - will create additional folders
             for GpC or HpC contexts.

-N   [INT]   Number of cells per chroms folder. (def = $cells_per_folder)
              If there are more cells than this value, it will
              produce multiple folders, each with the specified
              number of cells. This keeps files smaller and allows
              for better parallelization for downstream functions.
              Cells in folders will not match between CG and CH.
              
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
		push @CELLS, $l;
	} close IN;
	%PASSING_CELLS = ();
	@CELLS = shuffle(@CELLS);
	$folderID = sprintf("%03d", 1); $folderCT = 1;
	open OUT, ">$opt{'O'}.cellID_folderID.annot";
	foreach $cellID (@CELLS) {
		$folderCT++;
		if ($folderCT > $cells_per_folder) {
			$folderCT = 1; $folderID++;
		}
		$PASSING_CELLS{$cellID} = $folderID;
		print OUT "$cellID\t$folderID\n";
	}
	$last_folder = $folderID;
	close OUT;
}

if (defined $opt{'N'}) {$cells_per_folder = $opt{'N'}};
if (defined $opt{'M'}) {$mCH_max = $opt{'M'}};
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

if (defined $opt{'G'}) {
	if ($opt{'G'} =~ /gz$/) {
		open GC, "zcat $opt{'G'} |";
	} else {
		open GC, "$opt{'G'}";
	}
	while ($l = <GC>) {
		chomp $l;
		($chrom,$start,$end,$strand) = split(/\t/, $l);
		$GPC_CHR_POS{$chrom}{$start} = $strand;
		$inGpC++;
	} close GC;
}

if (defined $opt{'H'} || (!defined $opt{'H'} && !defined $opt{'C'})) {
	$ts = localtime(time);
	print LOG "$ts\tBeginning CH chrom file generation.\n";

	$bases_processed = 0;
	$increment = 10000000; $report = $increment;

	$CONTEXT_total{'H'}=0;
	$CONTEXT_total{'h'}=0;
	$CONTEXT_total{'X'}=0;
	$CONTEXT_total{'x'}=0;
	
	$folderID = sprintf("%03d", 1); $folderCT = 0; %CELLID_folderID = (); %HANDLES = ();
	%FOLDERID_cellIDs = ();
	@{$FOLDERID_cellIDs{$folderID}} = ();
	
	if (defined $opt{'L'}) {
		foreach $cellID (keys %PASSING_CELLS) {
			$CELLID_folderID{$cellID} = $PASSING_CELLS{$cellID};
		}
		for ($id = 1; $id <= $last_folder; $id++) {
			if (defined $opt{'G'}) {
				system("mkdir $opt{'O'}.HCH.$folderID.chroms");
				system("mkdir $opt{'O'}.GCH.$folderID.chroms");
				foreach $chrom (keys %CHROMS) {
					$handle = $folderID.$CHROMS{$chrom}."H";
					$HANDLES{$handle} = 1;
					open $handle, ">$opt{'O'}.HCH.$folderID.chroms/$chrom.bed";
					$handle = $folderID.$CHROMS{$chrom}."G";
					$HANDLES{$handle} = 1;
					open $handle, ">$opt{'O'}.GCH.$folderID.chroms/$chrom.bed";
				}
			} else {
				system("mkdir $opt{'O'}.CH.$folderID.chroms");
				foreach $chrom (keys %CHROMS) {
					$handle = $folderID.$CHROMS{$chrom};
					$HANDLES{$handle} = 1;
					open $handle, ">$opt{'O'}.CH.$folderID.chroms/$chrom.bed";
				}
			}
			$folderID++;
		}
	} else {
		if (defined $opt{'G'}) {
			system("mkdir $opt{'O'}.HCH.$folderID.chroms");
			system("mkdir $opt{'O'}.GCH.$folderID.chroms");
			foreach $chrom (keys %CHROMS) {
				$handle = $folderID.$CHROMS{$chrom}."H";
				$HANDLES{$handle} = 1;
				open $handle, ">$opt{'O'}.HCH.$folderID.chroms/$chrom.bed";
				$handle = $folderID.$CHROMS{$chrom}."G";
				$HANDLES{$handle} = 1;
				open $handle, ">$opt{'O'}.GCH.$folderID.chroms/$chrom.bed";
			}
		} else {
			system("mkdir $opt{'O'}.CH.$folderID.chroms");
			foreach $chrom (keys %CHROMS) {
				$handle = $folderID.$CHROMS{$chrom};
				$HANDLES{$handle} = 1;
				open $handle, ">$opt{'O'}.CH.$folderID.chroms/$chrom.bed";
			}
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
		
		if (!defined $CELLID_folderID{$cellID} && !defined $opt{'L'}) {
			$folderCT++;
			if ($folderCT > $cells_per_folder) { # new folder
				open FOLDER_LIST, ">$opt{'O'}.CH.$folderID.cellIDs.txt";
				foreach $printID (@{$FOLDERID_cellIDs{$folderID}}) {
					print FOLDER_LIST "$printID\n";
				} close FOLDER_LIST;
				$folderCT = 1;
				$folderID++;
				$CELLID_folderID{$cellID} = $folderID;
				push @{$FOLDERID_cellIDs{$folderID}}, $cellID;
				
				if (defined $opt{'G'}) {
					system("mkdir $opt{'O'}.HCH.$folderID.chroms");
					system("mkdir $opt{'O'}.GCH.$folderID.chroms");
					foreach $chrom (keys %CHROMS) {
						$handle = $folderID.$CHROMS{$chrom}."H";
						$HANDLES{$handle} = 1;
						open $handle, ">$opt{'O'}.HCH.$folderID.chroms/$chrom.bed";
						$handle = $folderID.$CHROMS{$chrom}."G";
						$HANDLES{$handle} = 1;
						open $handle, ">$opt{'O'}.GCH.$folderID.chroms/$chrom.bed";
					}
				} else {
					system("mkdir $opt{'O'}.CH.$folderID.chroms");
					foreach $chrom (keys %CHROMS) {
						$handle = $folderID.$CHROMS{$chrom};
						$HANDLES{$handle} = 1;
						open $handle, ">$opt{'O'}.CH.$folderID.chroms/$chrom.bed";
					}
				}
				
			} else {
				$CELLID_folderID{$cellID} = $folderID;
				push @{$FOLDERID_cellIDs{$folderID}}, $cellID;
			}
		}
		
		if (!defined $opt{'L'} || defined $PASSING_CELLS{$cellID}) {
			if (defined $CHROMS{$P[2]}) {
				if (defined $opt{'G'}) {
					if (defined $GPC_CHR_POS{$P[2]}{$P[3]}) {
						$handle = $CELLID_folderID{$cellID}.$CHROMS{$P[2]}."G";
					} else {
						$handle = $CELLID_folderID{$cellID}.$CHROMS{$P[2]}."H";
					}
				} else {
					$handle = $CELLID_folderID{$cellID}.$CHROMS{$P[2]};
				}
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
				
				if (defined $GPC_CHR_POS{$P[2]}{$P[3]}) {
					$CELLID_GCH_cov{$cellID}++;
					$CELLID_GCH_context_cov{$cellID}{$P[4]}++;
					$GCH_CONTEXT_total{$P[4]}++;
				}
			}
		}
	} close IN;

	
	foreach $handle (keys %HANDLES) {
		close $handle;
	}
	
	# print last folder cell list
	if (!defined $opt{'L'}) {
		open FOLDER_LIST, ">$opt{'O'}.CH.$folderID.cellIDs.txt";
		foreach $printID (@{$FOLDERID_cellIDs{$folderID}}) {
			print FOLDER_LIST "$printID\n";
		} close FOLDER_LIST;
	}
	
	if (defined $opt{'L'}) {
		open GC_VALS, ">$opt{'O'}.GmCH.vals";
		open HC_VALS, ">$opt{'O'}.HmCH.vals";
		open GC_COV, ">$opt{'O'}.GCH_cov.vals";
		open HC_COV, ">$opt{'O'}.HCH_cov.vals";
		foreach $cellID (keys %CELLID_totalCGcov) {
			# GCH first
			$meth = sprintf("%.3f", (($CELLID_GCH_context_cov{$cellID}{'H'}+$CELLID_GCH_context_cov{$cellID}{'X'})/$CELLID_GCH_cov{$cellID})*100);
			print GC_VALS "$cellID\t$meth\n";
			print GC_COV "$cellID\t$CELLID_GCH_cov{$cellID}\n";
			# HCG next
			$HC_cov = ($CELLID_totalCGcov{$cellID}-$CELLID_GCH_cov{$cellID});
			$meth = sprintf("%.3f", ((($CELLID_CONTEXT_cov{$cellID}{'H'}-$CELLID_GCH_context_cov{$cellID}{'H'})+($CELLID_CONTEXT_cov{$cellID}{'X'}-$CELLID_GCH_context_cov{$cellID}{'X'}))/$HC_cov)*100);
			print HC_VALS "$cellID\t$meth\n";
			print HC_COV "$cellID\t$HC_cov\n";
		} close GC_VALS; close GC_COV; close HC_VALS; close HC_COV;
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

	$folderID = sprintf("%03d", 1); $folderCT = 0; %CELLID_folderID = (); %HANDLES = ();
	%FOLDERID_cellIDs = ();
	@{$FOLDERID_cellIDs{$folderID}} = ();
	
	if (defined $opt{'L'}) {
		foreach $cellID (keys %PASSING_CELLS) {
			$CELLID_folderID{$cellID} = $PASSING_CELLS{$cellID};
		}
		for ($id = 1; $id <= $last_folder; $id++) {
			if (defined $opt{'G'}) {
				system("mkdir $opt{'O'}.HCG.$folderID.chroms");
				system("mkdir $opt{'O'}.GCG.$folderID.chroms");
				foreach $chrom (keys %CHROMS) {
					$handle = $folderID.$CHROMS{$chrom}."H";
					$HANDLES{$handle} = 1;
					open $handle, ">$opt{'O'}.HCG.$folderID.chroms/$chrom.bed";
					$handle = $folderID.$CHROMS{$chrom}."G";
					$HANDLES{$handle} = 1;
					open $handle, ">$opt{'O'}.GCG.$folderID.chroms/$chrom.bed";
				}
			} else {
				system("mkdir $opt{'O'}.CG.$folderID.chroms");
				foreach $chrom (keys %CHROMS) {
					$handle = $folderID.$CHROMS{$chrom};
					$HANDLES{$handle} = 1;
					open $handle, ">$opt{'O'}.CG.$folderID.chroms/$chrom.bed";
				}
			}
			$folderID++;
		}
	} else {
		if (defined $opt{'G'}) {
			system("mkdir $opt{'O'}.HCG.$folderID.chroms");
			system("mkdir $opt{'O'}.GCG.$folderID.chroms");
			foreach $chrom (keys %CHROMS) {
				$handle = $folderID.$CHROMS{$chrom}."H";
				$HANDLES{$handle} = 1;
				open $handle, ">$opt{'O'}.HCG.$folderID.chroms/$chrom.bed";
				$handle = $folderID.$CHROMS{$chrom}."G";
				$HANDLES{$handle} = 1;
				open $handle, ">$opt{'O'}.GCG.$folderID.chroms/$chrom.bed";
			}
		} else {
			system("mkdir $opt{'O'}.CH.$folderID.chroms");
			foreach $chrom (keys %CHROMS) {
				$handle = $folderID.$CHROMS{$chrom};
				$HANDLES{$handle} = 1;
				open $handle, ">$opt{'O'}.CG.$folderID.chroms/$chrom.bed";
			}
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
		
		if (!defined $CELLID_folderID{$cellID} && !defined $opt{'L'}) {
			$folderCT++;
			if ($folderCT > $cells_per_folder) { # new folder
				open FOLDER_LIST, ">$opt{'O'}.CG.$folderID.cellIDs.txt";
				foreach $printID (@{$FOLDERID_cellIDs{$folderID}}) {
					print FOLDER_LIST "$printID\n";
				} close FOLDER_LIST;
				$folderCT = 1;
				$folderID++;
				$CELLID_folderID{$cellID} = $folderID;
				push @{$FOLDERID_cellIDs{$folderID}}, $cellID;
				if (defined $opt{'G'}) {
					system("mkdir $opt{'O'}.HCG.$folderID.chroms");
					system("mkdir $opt{'O'}.GCG.$folderID.chroms");
					foreach $chrom (keys %CHROMS) {
						$handle = $folderID.$CHROMS{$chrom}."H";
						$HANDLES{$handle} = 1;
						open $handle, ">$opt{'O'}.HCG.$folderID.chroms/$chrom.bed";
						$handle = $folderID.$CHROMS{$chrom}."G";
						$HANDLES{$handle} = 1;
						open $handle, ">$opt{'O'}.GCG.$folderID.chroms/$chrom.bed";
					}
				} else {
					system("mkdir $opt{'O'}.CG.$folderID.chroms");
					foreach $chrom (keys %CHROMS) {
						$handle = $folderID.$CHROMS{$chrom};
						$HANDLES{$handle} = 1;
						open $handle, ">$opt{'O'}.CG.$folderID.chroms/$chrom.bed";
					}
				}
			} else {
				$CELLID_folderID{$cellID} = $folderID;
				push @{$FOLDERID_cellIDs{$folderID}}, $cellID;
			}
		}
		
		if (!defined $opt{'L'} || defined $PASSING_CELLS{$cellID}) {
			if (defined $CHROMS{$P[2]}) {
				if (defined $opt{'G'}) {
					if (defined $GPC_CHR_POS{$P[2]}{$P[3]}) {
						$handle = $CELLID_folderID{$cellID}.$CHROMS{$P[2]}."G";
					} else {
						$handle = $CELLID_folderID{$cellID}.$CHROMS{$P[2]}."H";
					}
				} else {
					$handle = $CELLID_folderID{$cellID}.$CHROMS{$P[2]};
				}
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
				
				if (defined $GPC_CHR_POS{$P[2]}{$P[3]}) {
					$CELLID_GCG_cov{$cellID}++;
					$CELLID_GCG_context_cov{$cellID}{$P[4]}++;
					$GCG_CONTEXT_total{$P[4]}++;
				}
			}
		}
	} close IN;

	foreach $handle (keys %HANDLES) {
		close $handle;
	}
	
	# print last folder cell list
	if (!defined $opt{'L'}) {
		open FOLDER_LIST, ">$opt{'O'}.CG.$folderID.cellIDs.txt";
		foreach $printID (@{$FOLDERID_cellIDs{$folderID}}) {
			print FOLDER_LIST "$printID\n";
		} close FOLDER_LIST;
	}

	if (defined $opt{'L'}) {
		open GC_VALS, ">$opt{'O'}.GmCG.vals";
		open HC_VALS, ">$opt{'O'}.HmCG.vals";
		open GC_COV, ">$opt{'O'}.GCG_cov.vals";
		open HC_COV, ">$opt{'O'}.HCG_cov.vals";
		foreach $cellID (keys %CELLID_totalCGcov) {
			# GCG first
			$meth = sprintf("%.3f", ($CELLID_GCG_context_cov{$cellID}{'Z'}/$CELLID_GCG_cov{$cellID})*100);
			print GC_VALS "$cellID\t$meth\n";
			print GC_COV "$cellID\t$CELLID_GCG_cov{$cellID}\n";
			# HCG next
			$HC_cov = ($CELLID_totalCGcov{$cellID}-$CELLID_GCG_cov{$cellID});
			$meth = sprintf("%.3f", (($CELLID_CONTEXT_cov{$cellID}{'Z'}-$CELLID_GCG_context_cov{$cellID}{'Z'})/$HC_cov)*100);
			print HC_VALS "$cellID\t$meth\n";
			print HC_COV "$cellID\t$HC_cov\n";
		} close GC_VALS; close GC_COV; close HC_VALS; close HC_COV;
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

