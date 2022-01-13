#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:1:2:t:b:R:smX", \%opt);

$hg38na = "/home/groups/oroaklab/refs/hg38/hs38d1_noalt_bismark/";
$mm10 = "/home/groups/oroaklab/refs/mm10/bismark/";
$hybrid = "/home/groups/oroaklab/refs/bismark_hybrid.hg19.mm10/";
$bowtie = "/home/groups/oroaklab/src/bowtie2/bowtie2-2.3.1/";
$threads = 1;

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
              Will use 2x -t for sort & merge if applicable
-b   [STR]   Path to bowtie2 (def = $bowtie)
-s           Do not proceed with sorting
-m           Do not proceed with merging (if -1 and -2 and no -s)
              Recommended if there are sorted bams from other
              sequencing runs that need to be merged.
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

open LOG, ">$opt{'O'}.align.log";
$ts = localtime(time);
print LOG "$ts\tAlignemnt called.\n";

if (defined $opt{'1'}) {
	$ts = localtime(time);
	print LOG "$ts Aligning read 1, Command:\n";
	$r1_align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R1 --unmapped --pbat -p $threads $ref $opt{'1'} >> $opt{'O'}.R1.align.log 2>> $opt{'O'}.R1.align.log";
	print LOG "\t$r1_align\n";
	system($r1_align);
	if (defined $opt{'X'}) {system("rm -f $opt{'O'}.R1.align.log")};
	system("cat $opt{'O'}.R1/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
}

if (defined $opt{'2'}) {
	$ts = localtime(time);
	print LOG "$ts Aligning read 2, Command:\n";
	$r2_align = "bismark --path_to_bowtie $bowtie -o $opt{'O'}.R2 --unmapped -p $threads $ref $opt{'2'} >> $opt{'O'}.R2.align.log 2>> $opt{'O'}.R2.align.log";
	print LOG "\t$r2_align\n";
	system($r2_align);
	if (defined $opt{'X'}) {system("rm -f $opt{'O'}.R2.align.log")};
	system("cat $opt{'O'}.R2/*SE_report.txt >> $opt{'O'}.alignment_stats.txt");
}

if (defined $opt{'1'} && !defined $opt{'s'}) {
	$ts = localtime(time);
	print LOG "$ts Name sorting read 1...\n";
	$threads2 = $threads*2;
	$sort1 = "samtools sort -@ $threads2 -n $opt{'O'}.R1/*.bam > $opt{'O'}.R1.nsrt.bam 2>> $opt{'O'}.merge.log";
	print LOG "\t$sort1\n";
	system($sort1);
	if (defined $opt{'X'}) {
		system("rm -f -R $opt{'O'}.R1/");
	}
}

if (defined $opt{'2'} && !defined $opt{'s'}) {
	$ts = localtime(time);
	print LOG "$ts Name sorting read 2...\n";
	$threads2 = $threads*2;
	$sort2 = "samtools sort -@ $threads2 -n $opt{'O'}.R2/*.bam > $opt{'O'}.R2.nsrt.bam 2>> $opt{'O'}.merge.log";
	print LOG "\t$sort2\n";
	system($sort2);
	if (defined $opt{'X'}) {
		system("rm -f -R $opt{'O'}.R2/");
	}
}

if (defined $opt{'1'} && defined $opt{'2'} && !defined $opt{'s'} && !defined $opt{'m'}) {
	$ts = localtime(time);
	print LOG "$ts Merging...\n";
	$merge = "samtools merge -@ $threads2 -n $opt{'O'}.nsrt.bam $opt{'O'}.R1.nsrt.bam $opt{'O'}.R2.nsrt.bam 2>> $opt{'O'}.merge.log";
	print LOG "\t$merge\n";
	if (defined $opt{'X'}) {
		system("rm -f -R $opt{'O'}.merge.log");
	}
}

exit;