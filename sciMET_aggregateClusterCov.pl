#!/usr/bin/perl

#chr10   101288606       66.67   1       2

use Getopt::Std; %opt = ();

getopts("O:F:A:C:R:f:v", \%opt);

$die = "

sciMET_aggregateClusterCov.pl -O [output folder] -A [annot for cells to clusters] -F [cov file folder]
   Optional: -C [context] - def = CG
             -R [chrom] - only run this chromosome (adds chrom name to output prefix)
             -f [cov/bed/both] output format (def = bed)

will generate a bed of sites covered by cluster

";

if (!defined $opt{'O'} ||
	!defined $opt{'F'} ||
	!defined $opt{'A'}) {die $die};

if (defined $opt{'C'}) {$context = $opt{'C'}} else {$context = "CG"};
if (defined $opt{'f'}) {$format = $opt{'f'}} else {$format = "bed"};

if ($format !~ /bed/i && $format !~ /both/i && $format !~ /cov/i) {die "\nFormat must be 'bed', 'cov' or 'both'\n"};
	
open IN, "$opt{'A'}";
while ($l = <IN>) {
	chomp $l;
	($cellID,$annot) = split(/\t/, $l);
	$ANNOT_CELLID{$annot}{$cellID} = 1;
} close IN;

if (defined $opt{'R'}) {
	system("mkdir $opt{'O'}.$opt{'R'}");
	open STATS, ">$opt{'O'}.$opt{'R'}/$opt{'O'}.$opt{'R'}.stats";
} else {
	system("mkdir $opt{'O'}");
	open STATS, ">$opt{'O'}/$opt{'O'}.stats";
}
if (defined $opt{'v'}) {$ts = localtime(time); print STDERR "$ts READING COV FILES\n"}; 
foreach $annot (keys %ANNOT_CELLID) {
	if (defined $opt{'v'}) {$ts = localtime(time); print STDERR "$ts READING CELLS for $annot\n"};
	%SITE_cov = (); $annot_cov = 0;
	foreach $cellID (keys %{$ANNOT_CELLID{$annot}}) {
		if (defined $opt{'v'}) {$ts = localtime(time); print STDERR "\t$ts READING: $cellID ...\n"};
		open IN, "$opt{'F'}/$cellID.$context.cov";
		while ($l = <IN>) {
			chomp $l;
			($chr,$pos,$meth,$C,$mC) = split(/\t/, $l);
			if (!defined $opt{'R'} || $chr eq $opt{'R'}) {
				$site = "$chr.$pos";
				if (!defined $SITE_cov{$site}) {$annot_cov++};
				$SITE_cov{$site} += ($C+$mC);
				$SITE_meth{$site} = $meth;
			}
		} close IN;
	}
	print STATS "$annot\r$annot_cov\n";
	if ($format =~ /both/i || $format =~ /bed/i) {
		if (defined $opt{'v'}) {$ts = localtime(time); print STDERR "$ts PRINTING BED FILES for $annot ...\n"};
		if (!defined $opt{'R'}) {
			open OUT, "| sort -k1,1 -k2,2n > $opt{'O'}/$annot.meth.bed";
		} else {
			open OUT, "| sort -k2,2n > $opt{'O'}.$opt{'R'}/$annot.$opt{'R'}.meth.bed";
		}
		foreach $site (keys %SITE_cov) {
			($chr,$pos) = split(/\./, $site);
			print OUT "$chr\t$pos\t$pos\t$SITE_cov{$site}\t$SITE_meth{$site}\n";
		} close OUT;
	}
	if ($format =~ /both/i || $format =~ /cov/i) {
		if (defined $opt{'v'}) {$ts = localtime(time); print STDERR "$ts PRINTING COV FILES for $annot ...\n"};
		if (!defined $opt{'R'}) {
			open OUT, "| sort -k1,1 -k2,2n > $opt{'O'}/$annot.$context.cov";
		} else {
			open OUT, "| sort -k2,2n > $opt{'O'}.$opt{'R'}/$annot.$opt{'R'}.$context.cov";
		}
		foreach $site (keys %SITE_cov) {
			($chr,$pos) = split(/\./, $site);
			$ctC = int($SITE_meth{$site}*$SITE_cov{$site});
			$ctT = $SITE_cov{$site}-$ctC;
			print OUT "$chr\t$pos\t$SITE_meth{$site}\t$ctT\t$ctC\n";
		} close OUT;
	}
} close STATS;
