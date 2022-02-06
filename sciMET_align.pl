#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:1:2:t:w:R:sXA:a:B:b:mk:l:", \%opt);

$hg38na = "/home/groups/oroaklab/refs/hg38/hs38d1_noalt_bismark/";
$mm10 = "/home/groups/oroaklab/refs/mm10/bismark/";
$hybrid = "/home/groups/oroaklab/refs/bismark_hybrid.hg19.mm10/";
$bowtie = "/home/groups/oroaklab/src/bowtie2/bowtie2-2.3.1/";
$threads = 1;
$r1_rounds = 4;
$r2_rounds = 2;
$r1_trim = 15;
$r2_trim = 15;

$die = "

sciMET_align.pl (options) -R [reference path] -O [output prefix] -1 [read1.trimmed.fq.gz] -2 [read2.trimmed.fq.gz]

Wrapper for bismark to run alignment of sciMETv2 reads.
reads 1 and 2 can be a list that is comma-separated.

Options:

-R   [STR]   Reference path (required)
             Shortcuts:
               hg38 = $hg38na
               mm10 = $mm10
               hyb = $hybrid
-O   [STR]   Output prefix (required)

-1   [STR]   Trimmed read 1*
-2   [STR]   Trimmed read 2*
             *Minimum of one read is required.
              Both required for sort & merge

-t   [INT]   Number of threads for bowtie2 to use
              Uses ~2x this number plus another ~2 threads total
              Will use same number for sort & merge if applicable

-A   [INT]   Number of alignment rounds for read 1 (def = $r1_rounds)
               Note: if 0, then no alignment will be performed.
-a   [INT]   Bases to trim for each subsequnt read 1 round (def = $r1_trim)
-B   [INT]   Number of alignment rounds for read 2 (def = $r2_rounds)
-b   [INT]   Bases to trim for each subsequent read 2 round (def = $r2_trim)

-w   [STR]   Path to bowtie2 (def = $bowtie)
-s           Do not proceed with sorting
-m           Do not proceed with merging (will merge if -1 and -2 are provided)
-X           Delete intermediate files (def = keep)

-k   [INT]   Pick up at specified read 1 alignment round
-l   [INT]   Pick up at specified read 2 alignment round
               Will pick up at trimming, specify the new
               unmapped fastq file to be trimmed. Must
               also specify a new -O.

";

if (!defined $opt{'R'}) {
	die "\nERROR: Provide a reference as -R\n$die";
} else {
	if ($opt{'R'} eq "hg38") {$ref = $hg38na}
	elsif ($opt{'R'} eq "mm10") {$ref = $mm10}
	elsif ($opt{'R'} eq "hyb") {$ref = $hybrid}
	else {$ref = $opt{'R'}};
}

if (!defined $opt{'1'} && !defined $opt{'2'}) {die "\nERROR: Specify at least one read as -1 or -2\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};
if (defined $opt{'w'}) {$bowtie = $opt{'w'}};
if (defined $opt{'e'}) {$trim_bases = $opt{'e'}};
if (defined $opt{'E'}) {$trim_rounds = $opt{'E'}};
if (defined $opt{'T'}) {$trim_reads = $opt{'T'}};

open LOG, ">$opt{'O'}.align.log";
$ts = localtime(time);
print LOG "$ts\tAlignment called.\n";


if (-d "$opt{'O'}.bams") {
	die "ERROR: $opt{'O'}.bams already exists! Delete or use a new output name!\n";
}

if (-d "$opt{'O'}.trim_reads") {
	die "ERROR: $opt{'O'}.trim_reads already exists! Delete or use a new output name!\n";
}

system("mkdir $opt{'O'}.bams");
system("mkdir $opt{'O'}.trim_reads");

if (defined $opt{'1'}) {
	$ts = localtime(time);
	print LOG "$ts Aligning read 1 for $r1_rounds rounds.\n";
	for ($round = 1; $round <= $r1_rounds; $round++) {
		if (defined $opt{'k'}) {$round = $opt{'k'}};
		if ($round > 1) {
			$ts = localtime(time);
			print LOG "$ts Trimming unaligned R1 reads by $r1_trim, round $round...\n";
			$prev_round = $round - 1;
			if (defined $opt{'k'}) {
				$trim = "seqtk trimfq -b 0 -e $r1_trim $opt{'1'} > $opt{'O'}.trim_reads/$opt{'O'}.R1.$prev_round.trimmed.fq";
			} else {
				$trim = "seqtk trimfq -b 0 -e $r1_trim $opt{'O'}.trim_reads/$opt{'O'}.R1.$prev_round.unmapped.fq.gz > $opt{'O'}.trim_reads/$opt{'O'}.R1.$prev_round.trimmed.fq";
			}
			print LOG "\t$trim\n";
			system($trim);
			$in_reads = "$opt{'O'}.trim_reads/$opt{'O'}.R1.$prev_round.trimmed.fq";
		} else {
			$in_reads = "$opt{'1'}";
		}
		
		$align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R1.$round --unmapped --pbat -p $threads $ref $in_reads >> $opt{'O'}.bismark.log 2>> $opt{'O'}.bismark.log";
		print LOG "\t$align\n";
		system($align);
		system("mv $opt{'O'}.R1.$round/*.bam $opt{'O'}.bams/$opt{'O'}.R1.$round.bam");
		system("mv $opt{'O'}.R1.$round/*.fq.gz $opt{'O'}.trim_reads/$opt{'O'}.R1.$round.unmapped.fq.gz");
		system("cat $opt{'O'}.R1.$round/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
		
		if (defined $opt{'X'}) {system("rm -f -R $opt{'O'}.R1.$round/")};
		
		if (!defined $opt{'s'}) {
			$sort = "samtools sort -@ $threads -n $opt{'O'}.bams/$opt{'O'}.R1.$round.bam > $opt{'O'}.bams/$opt{'O'}.R1.$round.nsrt.bam 2>> $opt{'O'}.align.log";
			print LOG "\t$sort\n";
			system($sort);
			if (defined $opt{'X'}) {system("rm -f $opt{'O'}.bams/$opt{'O'}.R1.$round.bam")};
		}
	}
}

if (defined $opt{'2'}) {
	$ts = localtime(time);
	print LOG "$ts Aligning read 2 for $r2_rounds rounds.\n";
	for ($round = 1; $round <= $r2_rounds; $round++) {
		if (defined $opt{'l'}) {$round = $opt{'l'}};
		if ($round > 1) {
			$ts = localtime(time);
			print LOG "$ts Trimming unaligned R2 reads by $r2_trim, round $round...\n";
			$prev_round = $round - 1;
			if (defined $opt{'l'}) {
				$trim = "seqtk trimfq -b 0 -e $r2_trim $opt{'2'} > $opt{'O'}.trim_reads/$opt{'O'}.R2.$prev_round.trimmed.fq";
			} else {
				$trim = "seqtk trimfq -b 0 -e $r2_trim $opt{'O'}.trim_reads/$opt{'O'}.R2.$prev_round.unmapped.fq.gz > $opt{'O'}.trim_reads/$opt{'O'}.R2.$prev_round.trimmed.fq";
			}
			print LOG "\t$trim\n";
			system($trim);
			$in_reads = "$opt{'O'}.trim_reads/$opt{'O'}.R2.$prev_round.trimmed.fq";
		} else {
			$in_reads = "$opt{'2'}";
		}
		
		$align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R2.$round --unmapped -p $threads $ref $in_reads >> $opt{'O'}.bismark.log 2>> $opt{'O'}.bismark.log";
		print LOG "\t$align\n";
		system($align);
		system("mv $opt{'O'}.R2.$round/*.bam $opt{'O'}.bams/$opt{'O'}.R2.$round.bam");
		system("mv $opt{'O'}.R2.$round/*.fq.gz $opt{'O'}.trim_reads/$opt{'O'}.R2.$round.unmapped.fq.gz");
		system("cat $opt{'O'}.R2.$round/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
		
		if (defined $opt{'X'}) {system("rm -f -R $opt{'O'}.R2.$round/")};
		
		if (!defined $opt{'s'}) {
			$sort = "samtools sort -@ $threads -n $opt{'O'}.bams/$opt{'O'}.R2.$round.bam > $opt{'O'}.bams/$opt{'O'}.R2.$round.nsrt.bam 2>> $opt{'O'}.align.log";
			print LOG "\t$sort\n";
			system($sort);
			if (defined $opt{'X'}) {system("rm -f $opt{'O'}.bams/$opt{'O'}.R2.$round.bam")};
		}
	}
}

open IN, "$opt{'O'}.alignment_stats.txt";
$r1_total = 0; $r2_total = 0;
$r1_aligned = 0; $r2_aligned = 0;
$r1_CG = 0; $r1_CHG = 0; $r1_CHH = 0;
$r1_mCG = 0; $r1_mCHG = 0; $r1_mCHH = 0;
$r2_CG = 0; $r2_CHG = 0; $r2_CHH = 0;
$r2_mCG = 0; $r2_mCHG = 0; $r2_mCHH = 0;
$all_CG = 0; $all_CHG = 0; $all_CHH = 0;
$all_mCG = 0; $all_mCHG = 0; $all_mCHH = 0;
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^Bismark report for/) {
		if ($l =~ /R1/) {$read = 1} elsif ($l =~ /R2/) {$read = 2};
	} elsif ($l =~ /^Sequences analysed in total/) {
		@P = split(/\s/); $num = pop(@P);
		if ($read == 1) {
			if ($num > $r1_total) {$r1_total = $num};
		} else {
			if ($num > $r2_total) {$r2_total = $num};
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
open OUT, ">$opt{'O'}.compiled_alignment.txt";
print OUT "#Compiled alignment report for $opt{'O'}

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
"; close OUT;

if (defined $opt{'1'} && defined $opt{'2'} && !defined $opt{'s'} && !defined $opt{'m'}) {
	$merge = "samtools merge -@ $threads -n $opt{'O'}.nsrt.bam $opt{'O'}.bams/*.nsrt.bam 2>> $opt{'O'}.align.log";
	print LOG "\t$merge\n";
	system($merge);
}

exit;