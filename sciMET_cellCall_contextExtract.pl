#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:X:N:B:C:o:d", \%opt);

$die = "

sciMET_cellCall_contextExtract.pl (options) [folder with cellcall files]

Options:

   -O  [DIR]  Output directory for new cellcall files (def = input directory)
   -C  [CG/CH] Context that will be searched. Can also search custom contexts. (req)
   -N  [STR]  Name of new context (e.g. CAC), (req)
   -B  [BED]  Bed file of sites, can be gzipped (req)
   -X  [STR]  Name of context for reciprocal (will only output the reciprocal
                if this option is specified)

Eg for GC methylaiton for accessibility saving to a new folder, run for CH and CG:
>sciMET_cellCall_contextExtract.pl -O MyOut -C CH -N GCH -B myGCsites.bed -X HCH MyIn
   and
>sciMET_cellCall_contextExtract.pl -O MyOut -C CG -N GCG -B myGCsites.bed -X HCG MyIn

Or for CAC methylation keeping regular CH files intact and adding to current folder:
>sciMET_cellCall_contextExtract.pl -C CH -N CAC -B myCACsites.bed MyIn

";

# parse input options
if (!defined $opt{'C'} ||
	!defined $opt{'N'} ||
	!defined $opt{'B'} ||
	!defined $ARGV[0]) {die $die};

# if issues with 0 or 1 based coordinates
if (!defined $opt{'o'}) {
	$offset = 0;
} else {
	$offset = $opt{'o'};
}

$ARGV[0] =~ s/\/$//;
if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
} else {
	$opt{'O'} =~ s/\/$//;
	system("mkdir $opt{'O'}");
}

# load bed sites
if ($opt{'B'} =~ /gz$/) {
	open IN, "zcat $opt{'B'} |";
} else {
	open IN, "$opt{'B'}";
}
%POS = ();
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	$P[1] += $offset;
	$pos = "$P[0].$P[1]";
	$POS{$pos} = 0;
	$sitesLoaded++;
} close IN;

if (defined $opt{'d'}) {
	print STDERR "$sitesLoaded sites loaded.\n";
}

# pull files from cellCalls folder
opendir(CCD, "$ARGV[0]");
@CELLCALLS = readdir(CCD);
closedir(CCD);

%BARC_sitesIncluded = (); %BARC_cov = (); %BARC_meth = ();
open INFO, ">$opt{'O'}.$opt{'N'}.cellInfo.txt";
if (defined $opt{'X'}) {
	%BARC_sitesExcluded = (); %BARC_covExcluded = (); %BARC_methExcluded = ();
	open XINFO, ">$opt{'O'}.$opt{'X'}.cellInfo.txt";
}
foreach $cellCall (@CELLCALLS) {
	($barc,$context,$cov) = split(/\./, $cellCall);
	if ($context eq $opt{'C'}) {
		$cellsSearched++;
		if (defined $opt{'d'}) {
			print STDERR "Readin cov file: $cellCall, cell barcode $barc\nRand 1% subset:\n";
		}
		if ($cellCall =~ /gz$/) {
			open IN, "zcat $ARGV[0]/$cellCall |" || die "ERROR! Cannot open $cellCall!\n";
		} else {
			open IN, "$ARGV[0]/$cellCall" || die "ERROR! Cannot open $cellCall!\n";
		}
		if (-e "$opt{'O'}/$barc.$opt{'N'}.cov") {
			die "\n\nERROR: $opt{'O'}/$barc.$opt{'N'}.cov file exists! Exiting so as not to overwrite!\n\n";
		}
		open OUT, "| sort -k1 -k 2,2n > $opt{'O'}/$barc.$opt{'N'}.cov";
		if (defined $opt{'X'}) {
			if (-e "$opt{'O'}/$barc.$opt{'X'}.cov") {
				die "\n\nERROR: $opt{'O'}/$barc.$opt{'X'}.cov file exists! Exiting so as not to overwrite!\n\n";
			}
			open XOUT, "| sort -k1 -k 2,2n > $opt{'O'}/$barc.$opt{'X'}.cov";
		}
		while ($l = <IN>) {
			chomp $l;
			($chr,$base,$pct,$C,$mC) = split(/\t/, $l);
			$pos = "$chr.$base";
			if (defined $opt{'d'}) {
				$rand = rand();
				if ($rand <0.01) {
					print STDERR "\t$barc\t$context\t$chr\t$base\t$pos\n";
				}
			}
			if (defined $POS{$pos}) {
				$BARC_sitesIncluded{$barc}++;
				$BARC_cov{$barc}+=($C+$mC);
				$BARC_meth{$barc}+=$mC;				
				print OUT "$l\n";
			} elsif (defined $opt{'X'}) {
				$BARC_sitesExcluded{$barc}++;
				$BARC_covExcluded{$barc}+=($C+$mC);
				$BARC_methExcluded{$barc}+=$mC;
				print XOUT "$l\n";
			}
		}
		close IN; close OUT;
		if ($BARC_cov{$barc}>0) {
			$mCpct = sprintf("%.2f", ($BARC_meth{$barc}/$BARC_cov{$barc})*100);
			print INFO "$barc\t$BARC_sitesIncluded{$barc}\t$BARC_cov{$barc}\t$mCpct\n";
		}
		if (defined $opt{'X'}) {
			close XOUT;
			if ($BARC_covExcluded{$barc}>0) {
				$mCpct = sprintf("%.2f", ($BARC_methExcluded{$barc}/$BARC_covExcluded{$barc})*100);
				print XINFO "$barc\t$BARC_sitesExcluded{$barc}\t$BARC_covExcluded{$barc}\t$mCpct\n";
			}
		}
	}
}
close INFO;
if (defined $opt{'X'}) {close XINFO};

exit;
