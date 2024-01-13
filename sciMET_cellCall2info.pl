#!/usr/bin/perl

$die = "

sciMET_cellCall2info.pl [cellCall folder] > [cellInfo output]

Will print CG then CH contexts then any additional contexts

";

if (!defined $ARGV[0]) {die $die};

opendir(CCD, "$ARGV[0]");
@CELLCALLS = readdir(CCD);
closedir(CCD);
$alt = 0;
foreach $cellCall (@CELLCALLS) { if ($cellCall =~ /cov/) {
	($cellID,$context,$cov) = split(/\./, $cellCall);
	print STDERR "DEBUG: file=$ARGV[0]/$cellCall cellID=$cellID\n";
	if ($context ne "CG" || $context ne "CH") {$alt++; push @ALTS, $context};
	if ($cellCall =~ /gz$/) {
		open IN, "zcat $ARGV[0]/$cellCall |" || die "ERROR! Cannot open $cellCall!\n";
	} else {
		open IN, "$ARGV[0]/$cellCall" || die "ERROR! Cannot open $cellCall!\n";
	}
	while ($l = <IN>) {
		chomp $l;
		($chr,$pos,$pct,$C,$mC) = split(/\t/, $l);
		$CELLID_ctx_cov{$cellID}{$context}+=($C+$mC);
		$CELLID_ctx_mC{$cellID}{$context}+=$mC;
		$CELLID_ctx_sites{$cellID}{$context}++;
		$CELLID_totalCov{$cellID}+=($C+$mC);
	} close IN;
}}

$header = "#CellID\tCoverage\tCG_Cov\tCG_mC_Pct\tCH_Cov\tCH_mC_Pct";
if ($alt>0) {
	for ($i = 0; $i < @ALTS; $i++) {
		$header .= "\t$ALTS[$i]_Cov\t$ALTS[$i]_mC_Pct";
	}
}

foreach $cellID (keys %CELLID_ctx_cov) {
	if (defined $CELLID_totalCov{$cellID}) {
		print "$cellID\t$CELLID_totalCov{$cellID}";
	} else {
		print "$cellID\t0";
	}
	if (defined $CELLID_ctx_cov{$cellID}{'CG'}) {
		$CGsites = $CELLID_ctx_sites{$cellID}{'CG'};
		$CGpct = sprintf("%.2f", ($CELLID_ctx_mC{$cellID}{'CG'}/$CELLID_ctx_cov{$cellID}{'CG'})*100);
	} else {
		$CGsites = 0; $CGpct = "0.00";
	}
	if (defined $CELLID_ctx_cov{$cellID}{'CH'}) {
		$CHsites = $CELLID_ctx_sites{$cellID}{'CH'};
		$CHpct = sprintf("%.2f", ($CELLID_ctx_mC{$cellID}{'CH'}/$CELLID_ctx_cov{$cellID}{'CH'})*100);
	} else {
		$CHsites = 0; $CHpct = "0.00";
	}
	print "\t$CGsites\t$CGpct\t$CHsites\t$CHpct";
	if ($alt > 0) {
		for ($i = 0; $i < @ALTS; $i++) {
			if (defined $CELLID_ctx_cov{$cellID}{$ALTS[$i]}) {
				$ALTsites = $CELLID_ctx_sites{$cellID}{$ALTS[$i]};
				$ALTpct = sprintf("%.2f", ($CELLID_ctx_mC{$cellID}{$ALTS[$i]}/$CELLID_ctx_cov{$cellID}{$ALTS[$i]})*100);
			} else {
				$ALTsites = 0;
				$ALTpct = "0.00";
			}
			print "\t$ALTsites\t$ALTpct";
		}
	}
	print "\n";
}
