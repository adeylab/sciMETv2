#!/usr/bin/perl

BEGIN {
	use FindBin '$RealBin';
	push @INC, $RealBin;
}

# defaults
$adapters = "$RealBin/sciMETv2_adapters.fa";
$trimmomatic = "/home/users/oconneru/bin/Trimmomatic-0.38/trimmomatic-0.38.jar";
$min_RL = 30;
$threads = 1;

getopts("O:A:1:2:T:m:t:", \%opt);

$die = "

sciMETv2_trim.pl (options) -O output_prefix -1 read1.fq.gz -2 read2.fq.gz

Runs trimmomatic with sciMET parameters.

Options:

-A   [STR]   Adapter fastq (def = $adapters)
-O   [STR]   Output Prefix
-1   [STR]   Read1 fastq
-2   [STR]   Read2 fastq
-T   [STR]   Trimmomatic jar file path (def = $trimmomatic)
-m   [INT]   Min read length (def = $min_RL)
-t   [INT]   Threads to use (def = $threads)

";

if (!defined $opt{'O'}) {die "\nSpecify output prefix as -O\n$die"};
if (!defined $opt{'1'} && !defined $opt{'2'}) {die "\nSpecify at least one read as -1 or -2\n$die"};
if (defined $opt{'A'}) {$adapters = $opt{'A'}};
if (defined $opt{'T'}) {$trimmomatic = $opt{'T'}};
if (defined $opt{'m'}) {$min_RL = $opt{'m'}};
if (defined $opt{'t'}) {$threads = $opt{'t'}};

if (defined $opt{'1'}) {
	$trim_R1 = "java -Xmx8G -jar $trimmomatic SE -threads $threads $opt{'1'} $opt{'O'}.trimmed.R1.fq.gz ILLUMINACLIP:$adapters:2:30:10 MINLEN:$min_RL >> $opt{'O'}.trimmed.R1.log 2>> $opt{'O'}.trimmed.R1.log";
	print STDERR "Trimming read 1. Command:\n\t$trim_R1\n";
	system("$trim_R1");
}

if (defined $opt{'2'}) {
	$trim_R2 = "java -Xmx8G -jar $trimmomatic SE -threads $threads $opt{'2'} $opt{'O'}.trimmed.R2.fq.gz ILLUMINACLIP:$adapters:2:30:10 MINLEN:$min_RL >> $opt{'O'}.trimmed.R2.log 2>> $opt{'O'}.trimmed.R2.log";
	print STDERR "Trimming read 2. Command:\n\t$trim_R2\n";
	system("$trim_R2");
}

exit;