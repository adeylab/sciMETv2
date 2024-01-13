#!/usr/bin/perl

$die = "

sciMET_compareDMRs.pl [truth bed of DMRs] [bed from collated cluster DMRs]

";

if (!defined $ARGV[1]) {die $die};

open IN, "cat $ARGV[0] $ARGV[1] | awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}' | sortBed | mergeBed -o distinct -c 4 |";
while ($l = <IN>) {
	chomp $l;
	@P = split(/\t/, $l);
	@SET = split(/,/, $P[3]);
	$set = "";
	foreach $item (sort @SET) {
		$set .= "$item&";
	}
	$set =~ s/&$//;
	$SET_ct{$set}++;
} close IN;



$euler_list = "set<-euler(c(";

foreach $set (keys %SET_ct) {
	$euler_list .= "\"$set\" = $SET_ct{$set},";
}

$euler_list =~ s/,$//;

print "

eulerr list:

$euler_list))

";
