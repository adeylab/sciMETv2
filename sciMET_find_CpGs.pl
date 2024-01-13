#!/usr/bin/perl

$die = "

sciMET_find_GpCs.pl [genome fasta].gz > [GC sites bed file]

";

if (!defined $ARGV[0]) {die $die};

if ($ARGV[0] =~ /gz$/) {
	open IN, "zcat $ARGV[0] |";
} else {
	open IN, "$ARGV[0]";
}
while ($l = <IN>) {
	chomp $l;
	if ($l =~ /^\>/) {
		$pos = 0;
		$chr = $l;
		$chr =~ s/^\>//;
		$chr =~ s/\s.+$//;
		$prevLast = "N";
	} else {
		@S = split(//, $l);
		$pos++;
		if ($S[0] =~ /g/i) {
			if ($prevLast =~ /c/i) {
				print "$chr\t".($pos-1)."\t".($pos-1)."\t-\n";
				print "$chr\t$pos\t$pos\t\+\n";
			}
		}
		for ($i = 1; $i < @S; $i++) {
			$pos++;
			if ($S[$i] =~ /g/i && $S[$i-1] =~ /c/i) {
				print "$chr\t".($pos-1)."\t".($pos-1)."\t-\n";
				print "$chr\t$pos\t$pos\t\+\n";
			}
		}
		$prevLast = $S[@S-1];
	}
} close IN;

exit;
