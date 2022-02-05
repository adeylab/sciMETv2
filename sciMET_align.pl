#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:1:2:t:w:R:sXA:a:B:b:m", \%opt);

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
print LOG "$ts\tAlignemnt called.\n";


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
		if ($round > 1) {
			$ts = localtime(time);
			print LOG "$ts Trimming unaligned R1 reads by $r1_trim, round $round...\n";
			$prev_round = $round - 1;
			$trim = "seqtk trimfq -b 0 -e $r1_trim $opt{'O'}.trim_reads/$opt{'O'}.R1.$prev_round.unmapped.fq.gz > $opt{'O'}.trim_reads/$opt{'O'}.R1.$prev_round.trimmed.fq";
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
		if ($round > 1) {
			$ts = localtime(time);
			print LOG "$ts Trimming unaligned R2 reads by $r2_trim, round $round...\n";
			$prev_round = $round - 1;
			$trim = "seqtk trimfq -b 0 -e $r2_trim $opt{'O'}.trim_reads/$opt{'O'}.R2.$prev_round.unmapped.fq.gz > $opt{'O'}.trim_reads/$opt{'O'}.R2.$prev_round.trimmed.fq";
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

if (defined $opt{'1'} && defined $opt{'2'} && !defined $opt{'s'} && !defined $opt{'m'}) {
	$merge = "samtools merge -@ $threads -n $opt{'O'}.nsrt.bam $opt{'O'}.bams/*.nsrt.bam 2>> $opt{'O'}.align.log";
	print LOG "\t$merge\n";
	system($merge);
}

exit;