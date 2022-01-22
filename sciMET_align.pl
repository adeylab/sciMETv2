#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:1:2:t:b:R:sXE:e:T:", \%opt);

$hg38na = "/home/groups/oroaklab/refs/hg38/hs38d1_noalt_bismark/";
$mm10 = "/home/groups/oroaklab/refs/mm10/bismark/";
$hybrid = "/home/groups/oroaklab/refs/bismark_hybrid.hg19.mm10/";
$bowtie = "/home/groups/oroaklab/src/bowtie2/bowtie2-2.3.1/";
$threads = 1;
$trim_rounds = 2;
$trim_bases = 15;
$trim_reads = 0;

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

-T   [0|1|2|12]  Trim and re-align unaliged reads in rounds.
              0 = none, 1 or 2 is only read 1 or 2, 12 = both (def = $trim_reads)
-E   [INT]   Will perform the specified number of rounds. (def = $trim_rounds)
-e   [INT]   Bases to end trim per round. (def = $trim_bases)

-b   [STR]   Path to bowtie2 (def = $bowtie)
-s           Do not proceed with sorting
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
if (defined $opt{'b'}) {$bowtie = $opt{'b'}};
if (defined $opt{'e'}) {$trim_bases = $opt{'e'}};
if (defined $opt{'E'}) {$trim_rounds = $opt{'E'}};
if (defined $opt{'T'}) {$trim_reads = $opt{'T'}};

open LOG, ">$opt{'O'}.align.log";
$ts = localtime(time);
print LOG "$ts\tAlignemnt called.\n";

if (defined $opt{'1'}) {
	$ts = localtime(time);
	print LOG "$ts Aligning read 1, Command:\n";
	$r1_align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R1 --unmapped --pbat -p $threads $ref $opt{'1'} >> $opt{'O'}.R1.align.log 2>> $opt{'O'}.R1.align.log";
	print LOG "\t$r1_align\n";
	system($r1_align);
	system("mv $opt{'O'}.R1/$opt{'1'}_unmapped_reads.fq.gz $opt{'O'}.R1.unmapped.fq.gz");
	system("cat $opt{'O'}.R1/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
	
	if (!defined $opt{'s'}) {
		$ts = localtime(time);
		print LOG "$ts Name sorting read 1...\n";
		$sort1 = "samtools sort -@ $threads -n $opt{'O'}.R1/*.bam > $opt{'O'}.R1.nsrt.bam 2>> $opt{'O'}.sort.log";
		print LOG "\t$sort1\n";
		system($sort1);
		if (defined $opt{'X'}) {system("rm -f -R $opt{'O'}.R1/")};
	}
	
	if ($trim_reads > 0) {
		for ($round = 1; $round <= $trim_rounds; $round++) {
			$ts = localtime(time);
			print LOG "$ts Trimming unaligned R1 reads by $trim_bases, round $round...\n";
			$trim = "seqtk trim -b 0 -e $trim_bases $opt{'O'}.R1.unmapped.fq.gz | gzip > $opt{'O'}.R1.unmapped.trim_$round.fq.gz";
			print LOG "\t$trim\n";
			system($trim);
			system("rm -f -R $opt{'O'}.R1.unmapped.fq.gz");
			
			$ts = localtime(time);
			print LOG "$ts Aligning after trim...\n";
			$r1_align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R1.trim_$round --unmapped --pbat -p $threads $ref $opt{'O'}.R1.unmapped.trim_$round.fq.gz >> $opt{'O'}.R1.align.log 2>> $opt{'O'}.R1.align.log";
			print LOG "\t$r1_align\n";
			system($r1_align);
			system("mv $opt{'O'}.R1.trim_$round/*_unmapped_reads.fq.gz $opt{'O'}.R1.unmapped.fq.gz");
			system("cat $opt{'O'}.R1.trim_$round/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
			
			if (!defined $opt{'s'}) {
				$ts = localtime(time);
				print LOG "$ts Name sorting read 1 trimmed, round $round...\n";
				$sort1 = "samtools sort -@ $threads -n $opt{'O'}.R1.trim_$round/*.bam > $opt{'O'}.R1.trim_$round.nsrt.bam 2>> $opt{'O'}.sort.log";
				print LOG "\t$sort1\n";
				system($sort1);
				if (defined $opt{'X'}) {system("rm -f -R $opt{'O'}.R1.trim_$round/")};
			}
			
		}
	}
	if (defined $opt{'X'}) {system("rm -f $opt{'O'}.R1.align.log $opt{'O'}.sort.log")};
}

if (defined $opt{'2'}) {
	$ts = localtime(time);
	print LOG "$ts Aligning read 2, Command:\n";
	$r2_align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R2 --unmapped -p $threads $ref $opt{'2'} >> $opt{'O'}.R2.align.log 2>> $opt{'O'}.R2.align.log";
	print LOG "\t$r2_align\n";
	system($r2_align);
	system("mv $opt{'O'}.R2/$opt{'2'}_unmapped_reads.fq.gz $opt{'O'}.R2.unmapped.fq.gz");
	system("cat $opt{'O'}.R2/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
	
	if (!defined $opt{'s'}) {
		$ts = localtime(time);
		print LOG "$ts Name sorting read 2...\n";
		$sort2 = "samtools sort -@ $threads -n $opt{'O'}.R2/*.bam > $opt{'O'}.R2.nsrt.bam 2>> $opt{'O'}.sort.log";
		print LOG "\t$sort2\n";
		system($sort2);
		if (defined $opt{'X'}) {
			system("rm -f -R $opt{'O'}.R2/");
		}
	}
	
	if ($trim_reads > 1) {
		for ($round = 1; $round <= $trim_rounds; $round++) {
			$ts = localtime(time);
			print LOG "$ts Trimming unaligned R1 reads by $trim_bases, round $round...\n";
			$trim = "seqtk trim -b 0 -e $trim_bases $opt{'O'}.R2.unmapped.fq.gz | gzip > $opt{'O'}.R2.unmapped.trim_$round.fq.gz";
			print LOG "\t$trim\n";
			system($trim);
			system("rm -f -R $opt{'O'}.R2.unmapped.fq.gz");
			
			$ts = localtime(time);
			print LOG "$ts Aligning after trim...\n";
			$r2_align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R2.trim_$round --unmapped --pbat -p $threads $ref $opt{'O'}.R2.unmapped.trim_$round.fq.gz >> $opt{'O'}.R2.align.log 2>> $opt{'O'}.R2.align.log";
			print LOG "\t$r2_align\n";
			system($r2_align);
			system("mv $opt{'O'}.R2.trim_$round/*_unmapped_reads.fq.gz $opt{'O'}.R2.unmapped.fq.gz");
			system("cat $opt{'O'}.R2.trim_$round/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
			
			if (!defined $opt{'s'}) {
				$ts = localtime(time);
				print LOG "$ts Name sorting read 2 trimmed, round $round...\n";
				$sort2 = "samtools sort -@ $threads -n $opt{'O'}.R2.trim_$round/*.bam > $opt{'O'}.R2.trim_$round.nsrt.bam 2>> $opt{'O'}.sort.log";
				print LOG "\t$sort2\n";
				system($sort2);
				if (defined $opt{'X'}) {system("rm -f -R $opt{'O'}.R2.trim_$round/")};
			}
			
		}
	}
	if (defined $opt{'X'}) {system("rm -f $opt{'O'}.R2.align.log $opt{'O'}.sort.log")};
}

exit;