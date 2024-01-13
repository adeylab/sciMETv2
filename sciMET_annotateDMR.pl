#!/usr/bin/perl

#chr1    870000  872000  noCap.1 0.00736569998486068     -37.5
#chr1    13201   13800   Enhancer

$die = "

sciMET_annotateDMR.pl [bed with classes] [DMRs] > [annotation added to DMR file]

DMR file:

chr  start  end  sample  cluster   qval   methDiff

Out adds: class from bed file or 'None'

";

if (!defined $ARGV[1]) {die $die};

open IN, "intersectBed -a $ARGV[1] -b $ARGV[0] -wa -wb |";

while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	
	$site = "$P[0].$P[1]";
	
	$SITE_class{$site} = $P[10];
	
} close IN;

open IN, "$ARGV[1]";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	
	$site = "$P[0].$P[1]";
	
	if (defined $SITE_class{$site}) {
		$class = $SITE_class{$site};
	} else {$class = "None"};
	
	print "$P[0]\t$P[1]\t$P[2]\t$P[3]\t$P[4]\t$P[5]\t$P[6]\t$class\n";
	
} close IN;
