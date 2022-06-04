#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("O:1:2:t:R:sr:u", \%opt);

$hg38na = "/home/groups/oroaklab/refs/hg38/hs38d1_noalt_bismark/";
$mm10 = "/home/groups/oroaklab/refs/mm10/bismark/";
$hybrid = "/home/groups/oroaklab/refs/bismark_hybrid.hg19.mm10/";
$threads = 1;


$die = "

sciMET_align_pe.pl (options) -R [reference path] -O [output prefix] -1 [read1.trimmed.fq.gz] -2 [read2.trimmed.fq.gz]

Wrapper for bismark to run alignment of sciMETv2 reads.
reads can be a list that is comma-separated.

Options:

-R   [STR]   Reference path (required)
             Shortcuts:
               hg38 = $hg38na
               mm10 = $mm10
               hyb = $hybrid
-O   [STR]   Output prefix (required)

-1   [STR]   Trimmed read 1 (paired, req)
-2   [STR]   Trimmed read 2 (paired, req)

-t   [INT]   Number of threads for bowtie2 to use
              Uses ~2x this number plus another ~2 threads total
              Will use same number for sort. (def = $threads)

-u           Retain unmapped reads (def = discard)
-s           Do not proceed with sorting
-r   [STR]   Report alignment stats to slack channel
              Requires 'slack' as cli callable

";

$start_time = localtime(time);

if (!defined $opt{'R'}) {
	die "\nERROR: Provide a reference as -R\n$die";
} else {
	if ($opt{'R'} eq "hg38") {$ref = $hg38na}
	elsif ($opt{'R'} eq "mm10") {$ref = $mm10}
	elsif ($opt{'R'} eq "hyb") {$ref = $hybrid}
	else {$ref = $opt{'R'}};
}

if (!defined $opt{'1'} || !defined $opt{'2'}) {die "\nERROR: Reads 1 and 2 MUST be specified!\n$die"};
if (!defined $opt{'O'}) {die "\nERROR: Specify output as -O\n$die"};
if (defined $opt{'t'}) {$threads = $opt{'t'}};

open LOG, ">$opt{'O'}.align.log";
$ts = localtime(time);
print LOG "$ts\tAlignment called.\n";

if (defined $opt{'u'}) { # align swapping reads 1 and 2 and run as standard directional (faster than pbat)
	$pe_align_call = "bismark -p $threads --local --unmapped -o $opt{'O'} $ref -1 $opt{'2'} -2 $opt{'1'} 2>> $opt{'O'}.align.log";
} else {
	$pe_align_call = "bismark -p $threads --local -o $opt{'O'} $ref -1 $opt{'2'} -2 $opt{'1'} 2>> $opt{'O'}.align.log";
}
print LOG "Command: $pe_align_call\n";
system("$pe_align_call");
system("mv $opt{'O'}/*.bam $opt{'O'}.unsorted.bam");
system("mv $opt{'O'}/*_report.txt $opt{'O'}.align_report.txt");


if (!defined $opt{'s'}) {
	$sort = "samtools sort -T $opt{'O'}.TMP -m 4G -@ $threads -n $opt{'O'}.unsorted.bam > $opt{'O'}.nsrt.bam 2>> $opt{'O'}.align.log";
	$ts = localtime(time);
	print LOG "$ts\tSorting.\n";
	system("$sort");
}

$ts = localtime(time);
print LOG "$ts\tDone.\n";

if (defined $opt{'r'}) {
	$message = "Alignment complete for $opt{'O'}\nStart time: $start_time\nEnd time: $ts\nCall: $pe_align_call\n";
	system("slack -F $opt{'O'}.align_report.txt -c \"$message\" $opt{'r'} >/dev/null 2>/dev/null");
}

exit;