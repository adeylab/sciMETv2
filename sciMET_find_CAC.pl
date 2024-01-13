#!/usr/bin/perl

$die = "

sciMET_find_CAC.pl [genome fasta].gz > [CAC sites bed file]

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
		$Nmer = "NNN";
	} else {
		@S = split(//, $l);	
		for ($i = 0; $i < @S; $i++) {
			$pos++;
			$Nmer =~ s/^.//;
			$Nmer .= uc($S[$i]);
			if ($Nmer eq "CAC") {
				$mBase = $pos-2;
				print "$chr\t$mBase\t$mBase\t\+\n";
			} elsif ($Nmer eq "GTG") {
				$mBase = $pos;
				print "$chr\t$mBase\t$mBase\t\-\n";
			}
		}
	}
} close IN;

exit;
