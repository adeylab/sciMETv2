#!/usr/bin/perl

$die = "

sciMET_find_DCM.pl [genome fasta].gz > [DCM sites bed file]

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
		$fivemer = "NNNNN";
	} else {
		@S = split(//, $l);	
		for ($i = 0; $i < @S; $i++) {
			$pos++;
			$fivemer =~ s/^.//;
			$fivemer .= uc($S[$i]);
			if ($fivemer eq "CCAGG" || $fivemer eq "CCTGG") {
				$mBase = $pos-3;
				print "$chr\t$mBase\t$mBase\t\+\n";
				$mBase = $pos-1;
				print "$chr\t$mBase\t$mBase\t\-\n";
			}
		}
	}
} close IN;

exit;
