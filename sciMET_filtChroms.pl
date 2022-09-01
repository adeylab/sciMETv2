#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("F:B:O:", \%opt);

$die = "

Usage:

sciMET_filtChroms.pl -F [chromosome meth call folder(s)] -O [output folder] -B [Bed of targets]

Filters chrom folders to only include calls within regions and outputs a single
chroms folder for those calls. Needs to be sorted.

Required options:

-F   [STR]   Folder with chromosome-split methylation calls for cells
             generated using sciMET_extract.pl
			 Can be multiple folders comma-separated.
             Or can be a single text file (end in .txt) with list of folders
-O   [STR]   Output prefix, creates or uses directory.
-B   [STR]   Bed file with targets to restrict chrom meth calls to.

";

if (!defined $opt{'F'} ||
	!defined $opt{'O'} ||
	!defined $opt{'B'}) {die $die};
	
if ($opt{'F'} =~ /txt$/) {
	@FOLDERS = ();
	open F, "$opt{'F'}";
	while ($l = <F>) {chomp $l; push @FOLDERS, $l};
	close F;
} else {
	@FOLDERS = split(/,/, $opt{'F'});
}

for ($i = 1; $i <= 22; $i++) {$CHRS{"chr$i"} = 1}; $CHRS{"chrX"} = 1;
system("mkdir $opt{'O'}");

foreach $chr_dir (@FOLDERS) {
	foreach $chr (keys %CHRS) {
		if (-e "$chr_dir/$chr.bed.gz") {
			$chr_file = "$chr_dir/$chr.bed.gz";
			open OUT, ">> $opt{'O'}/$chr.bed";
			print STDERR "Running intersect: bedtools intersect -a $chr_file -b $opt{'B'} -wa \n";
			open INT, "bedtools intersect -a $chr_file -b $opt{'B'} -wa |";
			while ($l = <INT>) {
				print OUT $l;
			} close INT;
			close OUT;
		} else {
			print STDERR "WARNING: Cannot find $chr.bed.gz in target directory $chr_dir! Skipping!\n";
		}
	}
}