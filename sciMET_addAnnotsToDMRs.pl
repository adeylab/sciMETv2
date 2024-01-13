#!/usr/bin/perl

$die = "

sciMET_addAnnotsToDMRs.pl [DMR bed file with qval & cluster etc] [output file]

Adds columns for Gencode category, Genehancer class

";

if (!defined $ARGV[1]) {die $die};

system("perl /home/users/adey/adey-src/add_annot_to_bed.pl $ARGV[0] /home/groups/ravnica/refs/hg38/hg38.gencode_regions.bed 4 $ARGV[1].tmp top");
system("perl /home/users/adey/adey-src/add_annot_to_bed.pl $ARGV[1].tmp /home/groups/ravnica/refs/hg38/hg38_genehancer.class.bed 4 $ARGV[1] top");
system("rm -f $ARGV[1].tmp");
