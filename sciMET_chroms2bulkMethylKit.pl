#!/usr/bin/perl

$die = "

sciMET_chroms2bulkMethylKit.pl [input chroms file] > [output methylKit file]

Chrom files must be gzipped and sorted!

";

if (!defined $ARGV[0]) {
	die $die;
}

opendir(CHR, "$ARGV[0]");
@CHRS = readdir(CHR);
closedir(CHR);

print "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n";

foreach $chrFile (sort @CHRS) {
	if ($chrFile =~ /gz$/) {
		open IN, "zcat $ARGV[0]/$chrFile |";
		$prevPos = "null";
		while ($l = <IN>) {
			chomp $l;
			($chr,$pos,$pos2,$cellID,$call) = split(/\t/, $l);
			if ($pos ne $prevPos) {
				if ($prevPos ne "null") {
					$cov = ($C+$mC);
					if ($cov >= 10) {
						$freqC = sprintf("%.2f", ($mC/$cov)*100);
						$freqT = sprintf("%.2f", ($C/$cov)*100);
						print "$chr.$pos\t$chr\t$pos\tF\t$cov\t$freqC\t$freqT\n";
					}
				}
				$prevPos = $pos;
				$C = 0; $mC = 0;
			}
			if ($call eq "Z") {$mC++} elsif ($call eq "z") {$C++};
		}
		close IN;
	}
}
