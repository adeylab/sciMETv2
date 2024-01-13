#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:A:F:O:", \%opt);

$die = "

sciMET_Annot2methylKit.pl -F [CG methyl call folder] -O [output prefix] -A [cluster annot file]

Gets methylation calls into a format that methylKit can use for DMR detection on a
cluster level.

Options:
   -O   [STR]   Output prefix (req)
   -F   [STR]   Folder with chromosome-split methylation calls for cells (req)
                  generated using sciMET_extract.pl - use CG one only
   -A   [STR]   Cluster annotation file (req)
   
";

if (!defined $opt{'F'} ||
	!defined $opt{'A'} ||
	!defined $opt{'O'}) {die $die};
	
open IN, "$opt{'A'}";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$CELLID_annot{$cellID} = $annot;
	$ANNOTS{$annot}++;
} close IN;

$header = "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT";
foreach $annot (keys %ANNOTS) {
	open OUT, ">$opt{'O'}.$annot.methylKit.txt"; print OUT "$header\n"; close OUT;
}
	
@chr_files = <$opt{'F'}/*.gz>;

#chr1    10486   10486   TACCGATAGAGTAGTAGTCCGCGCAAGC    Z
# out =
# chrBase	    chr	    base	strand	coverage	freqC	freqT
# chr21.9764539	chr21	9764539	   R	12	        25.00	75.00

foreach $chr_file (@chr_files) {
	%ANNOT_POS_cov = (); %ANNOT_POS_C = ();
	open IN, "zcat $chr_file |";
	while ($l = <IN>) {
		chomp $l;
		($chr,$pos,$null,$cellID,$call) = split(/\t/, $l);
		$annot = $CELLID_annot{$cellID};
		$chrBase = "$chr.$pos";
		$ANNOT_POS_cov{$annot}{$chrBase}++;
		if ($call =~ /Z/) {
			$ANNOT_POS_C{$annot}{$chrBase}++;
		}
	}
	close IN;

	foreach $annot (keys %ANNOTS) {
		open OUT, "| sort -k 2,2 -k 3,3n >> $opt{'O'}.$annot.methylKit.txt";
		foreach $chrBase (keys %{$ANNOT_POS_cov{$annot}}) {
			($chr,$base) = split(/\./, $chrBase);
			if (!defined $ANNOT_POS_C{$annot}{$chrBase}) {$ANNOT_POS_C{$annot}{$chrBase} = 0};
			$freqC = sprintf("%.2f", ($ANNOT_POS_C{$annot}{$chrBase}/$ANNOT_POS_cov{$annot}{$chrBase})*100);
			$freqT = sprintf("%.2f", (($ANNOT_POS_cov{$annot}{$chrBase}-$ANNOT_POS_C{$annot}{$chrBase})/$ANNOT_POS_cov{$annot}{$chrBase})*100);
			print OUT "$chrBase\t$chr\t$base\tF\t$ANNOT_POS_cov{$annot}{$chrBase}\t$freqC\t$freqT\n";
		}
		close OUT;
	}
}
