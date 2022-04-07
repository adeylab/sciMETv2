#!/usr/bin/perl

$die = "

sciMET_split_GpC_calls.pl [GpC bed file] [chroms folder in] (chroms folder in 2) ...

Will output two new chroms folders using the input as the
output prefix for each chroms folder.

";

if (!defined $ARGV[1]) {die $die};

open GC, "$ARGV[0]";
while ($l = <GC>) {
	chomp $l;
	($chrom,$start,$end,$strand) = split(/\t/, $l);
	$CHR_POS{$chrom}{$start} = $strand;
	$inGpC++;
} close GC;

for ($i = 1; $i < @ARGV; $i++) {

	$C_count = 0; $inC_meth = 0; $GpC_count = 0; $GpC_meth = 0; $HpC_count = 0; $HpC_meth = 0;
	
	opendir(CHROMS, $ARGV[$i]) or die "Cannot open $ARGV[$i] chroms folder!\n";

	$pfx = $ARGV[$i];
	$pfx =~ s/\.chroms$//;
	system("mkdir $pfx.GC.chroms");
	system("mkdir $pfx.HC.chroms");

	while ($chrom_file = readdir(CHROMS)) {
		if ($chrom_file =~ /\.bed.gz$/) {
			open GC, "| gzip > $pfx.GC.chroms/$chrom_file";
			open HC, "| gzip > $pfx.HC.chroms/$chrom_file";
			open IN, "zcat $ARGV[$i]/$chrom_file |" or die "Cannot open $ARGV[$i]/$chrom_file!!\n";
			while ($l = <IN>) {
				chomp $l;
				@P = split(/\t/, $l);
				$C_count++;
				if ($P[4] =~ /[HXZ]/) {
					$inC_meth++;
				}
				if (defined $CHR_POS{$P[0]}{$P[1]}) {
					$GpC_count++;
					print GC "$l\n";
					if ($P[4] =~ /[HXZ]/) {
						$GpC_meth++;
					}
				} else {
					$HpC_count++;
					print HC "$l\n";
					if ($P[4] =~ /[HXZ]/) {
						$HpC_meth++;
					}
				}
			} close IN; close GC; close HC;
		}
	}
	closedir(CHROMS);

$inC_meth_pct = sprintf("%.2f", ($inC_meth/$C_count)*100);
$HpC_meth_pct = sprintf("%.2f", ($HpC_meth/$HpC_count)*100);
$GpC_meth_pct = sprintf("%.2f", ($GpC_meth/$GpC_count)*100);
$GpC_pct = sprintf("%.2f", ($GpC_count/$C_count)*100);

open OUT, ">$pfx.GpC_stats.txt";
print OUT "Chroms folder: $ARGV[$i]
Bed of GpC sites: $ARGV[0]
GpC sites loaded: $inGpC
Total coverage of input: $C_count
Methylation of input: $inC_meth_pct%
GpC site coverage: $GpC_count ($GpC_pct%)
GpC methylation level: $GpC_meth_pct%
HpC methylation level: $HpC_meth_pct%\n";
close OUT;

}

exit;