#!/usr/bin/perl

$die = "

sciMET_perCellAlignRate.pl \
	[read 1 AND 2 fastq file(s), post-trim - can be comma separated OR files from trim output on R1 and R2 stats - also comma sep; can be a mix] \
	[aligned bam file(s), pre-rmdup OR complexity.txt file(s), can be comma separated an a mix of types] \
	[list fo cellIDs to include] \
	[output file]

Calculate per-cell alignment rates.

";

if (!defined $ARGV[3]) {die $die};

open IN, "$ARGV[2]";
while ($l = <IN>) {
	chomp $l;
	$l =~ s/\s.+$//;
	$CELLS_ct{$l} = 0;
} close IN;

@FQS = split(/,/, $ARGV[0]);
foreach $fq (@FQS) {
	if ($fq =~ /gz$/) { # gzipped fastq
		open IN, "zcat $fq |";
		while ($tag = <IN>) {
			chomp $tag; $tag =~ s/^@//; $tag =~ s/:.+$//;
			if (defined $CELLS_ct{$tag}) {
				$CELLS_ct{$tag}++;
			}
			$null = <IN>; $null = <IN>; $null = <IN>;
		} close IN;
	} elsif ($fq =~ /txt$/) {
		open IN, "$fq";
		while ($l = <IN>) {
			chomp $l;
			($barc,$raw,$trim,$pct) = split(/\t/, $l);
			if (defined $CELLS_ct{$barc}) {
				$CELLS_ct{$barc}+=$trim;
			}
		}
	} else {
		print STDERR "ERROR: Cannot determine file type of $fq - skipping\n";
	}
}

@BAMS = split(/,/, $ARGV[1]);
foreach $bam (@BAMS) {
	if ($bam =~ /bam$/) {
		open IN, "samtools view $bam |";
		while ($l = <IN>) {
			chomp $l;
			@P = split(/\t/, $l);
			$barc = $P[0]; $barc =~ s/:.+//;
			if (defined $CELLS_ct{$barc}) {
				$CELLS_aln{$barc}++;
			}
		} close IN;
	} elsif ($bam =~ /complexity/) {
		open IN, "$bam";
		while ($l = <IN>) {
			chomp $l;
			($num,$barc,$aln,$uniq,$pct) = split(/\t/, $l);
			if (defined $CELLS_ct{$barc}) {
				$CELLS_aln{$barc}+=$aln;
			}
		}
	} else {
		print STDERR "ERROR: Cannot determine bam vs complexity for $bam! skipping!\n";
	}
}

open OUT, ">$ARGV[3]";
foreach $cellID (keys %CELLS_aln) {
	if ($CELLS_ct{$cellID}>0) {
		$pct = sprintf("%.2f", ($CELLS_aln{$cellID}/$CELLS_ct{$cellID})*100);
		print OUT "$cellID\t$CELLS_ct{$cellID}\t$CELLS_aln{$cellID}\t$pct\n";
	}
} close OUT;