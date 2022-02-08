#!/usr/bin/perl

$die = "
sciMET_aggAlign.pl [alignment_stats_file] [file for separate run] ...

Each file must only include reports from iterations
of a single run, not multiple. (ie all iterations of
reads 1 and 2 from the same sequencing run)

";

if (!defined $ARGV[1]) {die $die};

$r1_total = 0; $r2_total = 0;
$r1_aligned = 0; $r2_aligned = 0;
$r1_CG = 0; $r1_CHG = 0; $r1_CHH = 0;
$r1_mCG = 0; $r1_mCHG = 0; $r1_mCHH = 0;
$r2_CG = 0; $r2_CHG = 0; $r2_CHH = 0;
$r2_mCG = 0; $r2_mCHG = 0; $r2_mCHH = 0;
$all_CG = 0; $all_CHG = 0; $all_CHH = 0;
$all_mCG = 0; $all_mCHG = 0; $all_mCHH = 0;

for ($i = 0; $i < @ARGV; $i++) {
	$r1_file_total = 0; $r2_file_total = 0;
	open IN, "$ARGV[0]";
	while ($l = <IN>) {
		chomp $l;
		if ($l =~ /^Bismark report for/) {
			if ($l =~ /R1/) {$read = 1} elsif ($l =~ /R2/) {$read = 2};
		} elsif ($l =~ /^Sequences analysed in total/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				if ($num > $r1_file_total) {$r1_file_total = $num};
			} else {
				if ($num > $r2_file_total) {$r2_file_total = $num};
			}
		} elsif ($l =~ /^Number of alignments with a unique/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_aligned += $num;
			} else {
				$r2_aligned += $num;
			}
		} elsif ($l =~ /^Total methylated C's in CpG context/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_mCG += $num; $all_mCG += $num;
				$r1_CG += $num; $all_CG += $num;
			} else {
				$r2_mCG += $num; $all_mCG += $num;
				$r2_CG += $num; $all_CG += $num;;
			}
		} elsif ($l =~ /^Total methylated C's in CHG context/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_mCHG += $num; $all_mCHG += $num;
				$r1_CHG += $num; $all_CHG += $num;
			} else {
				$r2_mCHG += $num; $all_mCHG += $num;
				$r2_CHG += $num; $all_CHG += $num;;
			}
		} elsif ($l =~ /^Total methylated C's in CHH context/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_mCHH += $num; $all_mCHH += $num;
				$r1_CHH += $num; $all_CHH += $num;
			} else {
				$r2_mCHH += $num; $all_mCHH += $num;
				$r2_CHH += $num; $all_CHH += $num;;
			}
		} elsif ($l =~ /^Total unmethylated C's in CpG context/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_CG += $num; $all_CG += $num;
			} else {
				$r2_CG += $num; $all_CG += $num;;
			}
		} elsif ($l =~ /^Total unmethylated C's in CHG context/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_CHG += $num; $all_CHG += $num;
			} else {
				$r2_CHG += $num; $all_CHG += $num;;
			}
		} elsif ($l =~ /^Total unmethylated C's in CHH context/) {
			@P = split(/\s/); $num = pop(@P);
			if ($read == 1) {
				$r1_CHH += $num; $all_CHH += $num;
			} else {
				$r2_CHH += $num; $all_CHH += $num;;
			}
		}
	} close IN;
	$r1_total += $r1_file_total;
	$r2_total += $r2_file_total;
}

$all_total = $r1_total+$r2_total;
$all_aligned = $r1_aligned+$r2_aligned;
$all_pct = sprintf("%.2f", ($all_aligned/$all_total)*100);
$r1_pct = sprintf("%.2f", ($r1_aligned/$r1_total)*100);
$r2_pct = sprintf("%.2f", ($r2_aligned/$r2_total)*100);
$r1_mCG_pct = sprintf("%.2f", ($r1_mCG/$r1_CG)*100);
$r1_mCHG_pct = sprintf("%.2f", ($r1_mCHG/$r1_CHG)*100);
$r1_mCHH_pct = sprintf("%.2f", ($r1_mCHH/$r1_CHH)*100);
$r2_mCG_pct = sprintf("%.2f", ($r2_mCG/$r2_CG)*100);
$r2_mCHG_pct = sprintf("%.2f", ($r2_mCHG/$r2_CHG)*100);
$r2_mCHH_pct = sprintf("%.2f", ($r2_mCHH/$r2_CHH)*100);
$all_mCG_pct = sprintf("%.2f", ($all_mCG/$r1_CG)*100);
$all_mCHG_pct = sprintf("%.2f", ($all_mCHG/$all_CHG)*100);
$all_mCHH_pct = sprintf("%.2f", ($all_mCHH/$all_CHH)*100);

print "#Compiled alignment report for $opt{'O'}

READ1 total: $r1_total\taligned: $r1_aligned ($r1_pct %)
   CG total: $r1_CG\tmCG: $r1_mCG ($r1_mCG_pct %)
  CHG total: $r1_CHG\tmCHG: $r1_mCHG ($r1_mCHG_pct %)
  CHH total: $r1_CHH\tmCHH: $r1_mCHH ($r1_mCHH_pct %)
  
READ2 total: $r2_total\taligned: $r2_aligned ($r2_pct %)
   CG total: $r2_CG\tmCG: $r2_mCG ($r2_mCG_pct %)
  CHG total: $r2_CHG\tmCHG: $r2_mCHG ($r2_mCHG_pct %)
  CHH total: $r2_CHH\tmCHH: $r2_mCHH ($r2_mCHH_pct %)
  
ALL   total: $all_total\taligned: $all_aligned ($all_pct %)
   CG total: $all_CG\tmCG: $all_mCG ($all_mCG_pct %)
  CHG total: $all_CHG\tmCHG: $all_mCHG ($all_mCHG_pct %)
  CHH total: $all_CHH\tmCHH: $all_mCHH ($all_mCHH_pct %)
";

exit;