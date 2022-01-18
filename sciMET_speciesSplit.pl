#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:B:1:2:", \%opt);

$die = "

sciMET_speciesSplit.pl -B [barnyard_compare cell stats file] -1 [trimmed read 1] -2 [trimmed read 2] -O [out prefix]

All options required - will split fastq files by species.

";

if (!defined $opt{'O'} ||
	!defined $opt{'1'} ||
	!defined $opt{'2'} ||
	!defined $opt{'B'}) {
		die $die;
}

open IN, "$opt{'B'}";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$CELLID_species{$P[0]} = lc($P[5]);
} close IN;

$assigned1 = 0; $total1 = 0;
open H, "| gzip > $opt{'O'}.human.R1.fq.gz"; $human1 = 0;
open M, "| gzip > $opt{'O'}.mouse.R1.fq.gz"; $mouse1 = 0;
open X, "| gzip > $opt{'O'}.mixed.R1.fq.gz"; $mixed1 = 0;
open U, "| gzip > $opt{'O'}.unassigned.R1.fq.gz";

open IN, "zcat $opt{'1'} |";
while ($tag = <IN>) {
	$total1++;
	chomp $tag;
	$cellID = $tag; $cellID =~ s/:.+$//; $cellID =~ s/^\@//;
	$read = <IN>; chomp $read;
	$null = <IN>;
	$qual = <IN>; chomp $qual;
	if (defined $CELLID_species{$cellID}) {
		$assigned1++;
		if ($CELLID_species{$cellID} =~ /human/) {
			$human1++;
			print H "$tag\t$read\n\+\n$qual\n";
		} elsif ($CELLID_species{$cellID} =~ /mouse/) {
			$mouse1++;
			print M "$tag\t$read\n\+\n$qual\n";
		} elsif ($CELLID_species{$cellID} =~ /mixed/) {
			$mixed1++;
			print X "$tag\t$read\n\+\n$qual\n";
		} else {
			print U "$tag\t$read\n\+\n$qual\n";
		}
	} else {
		print U "$tag\t$read\n\+\n$qual\n";
	}
} close IN;
close H; close M; close X; close U;

if ($assigned1 < 1) {die "ERROR: ZERO reads were assigned for read 1! Check input files\n"};

$assigned2 = 0; $total2 = 0;
open H, "| gzip > $opt{'O'}.human.R2.fq.gz"; $human2 = 0;
open M, "| gzip > $opt{'O'}.mouse.R2.fq.gz"; $mouse2 = 0;
open X, "| gzip > $opt{'O'}.mixed.R2.fq.gz"; $mixed2 = 0;
open U, "| gzip > $opt{'O'}.unassigned.R2.fq.gz";

open IN, "zcat $opt{'2'} |";
while ($tag = <IN>) {
	$total1++;
	chomp $tag;
	$cellID = $tag; $cellID =~ s/:.+$//; $cellID =~ s/^\@//;
	$read = <IN>; chomp $read;
	$null = <IN>;
	$qual = <IN>; chomp $qual;
	if (defined $CELLID_species{$cellID}) {
		$assigned1++;
		if ($CELLID_species{$cellID} =~ /human/) {
			$human2++;
			print H "$tag\t$read\n\+\n$qual\n";
		} elsif ($CELLID_species{$cellID} =~ /mouse/) {
			$mouse2++;
			print M "$tag\t$read\n\+\n$qual\n";
		} elsif ($CELLID_species{$cellID} =~ /mixed/) {
			$mixed2++;
			print X "$tag\t$read\n\+\n$qual\n";
		} else {
			print U "$tag\t$read\n\+\n$qual\n";
		}
	} else {
		print U "$tag\t$read\n\+\n$qual\n";
	}
} close IN;
close H; close M; close X; close U;

$assigned1_pct = sprintf("%.2f", ($assigned1/$total1)*100);
$human1_pct = sprintf("%.2f", ($human1/$assigned1)*100);
$mouse1_pct = sprintf("%.2f", ($mouse1/$assigned1)*100);
$mixed1_pct = sprintf("%.2f", ($mixed1/$assigned1)*100);
$assigned2_pct = sprintf("%.2f", ($assigned2/$total2)*100);
$human2_pct = sprintf("%.2f", ($human2/$assigned2)*100);
$mouse2_pct = sprintf("%.2f", ($mouse2/$assigned2)*100);
$mixed2_pct = sprintf("%.2f", ($mixed2/$assigned2)*100);

open LOG, ">$opt{'O'}.species_split.txt";
print LOG "Total R1 reads: $total1
Assigned R1: $assigned1 ($assigned1_pct)
Human assigned: $human1 ($human1_pct)
Mouse assigned: $mouse1 ($mouse1_pct)
Mixed assigned: $mixed1 ($mixed1_pct)

Total R2 reads: $total2
Assigned R2: $assigned2 ($assigned2_pct)
Human assigned: $human2 ($human2_pct)
Mouse assigned: $mouse2 ($mouse2_pct)
Mixed assigned: $mixed2 ($mixed2_pct)
"; close LOG;

exit;