#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("w:h:O:", \%opt);

$width = 5;
$height = 4;

$die = "

plotGeneByAnnot.pl geneByAnnot.meth.txt

Will make a dotplot of targets x annotations from the file
produced by getGeneMeth.pl

Options:
   -w   [INT]   Width of plot (def = $width)
   -h   [INT]   Height of plot (def = $height)
   -O   [STR]   Output prefix (def = input file prefix)
   
";

if (!defined $ARGV[0]) {die $die};

if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.txt$//;
}

open R, ">$opt{'O'}.plot.r";
print R "library(ggplot2)
IN<-read.table(\"\$ARGV[0]\")
PLT<-ggplot(data=IN) + theme_bw() +
	 geom_point(aes(V2,V1,size=V3,color=V4)) +
	 scale_color_distiller(palette = \"Spectral\") +
	 xlab(\"Cluster\") + ylab(\"Gene\") + labs(color=\"z-score\",size=\"mCG\")
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=$width,height=$height)
"; close R;

system("Rscript $opt{'O'}.plot.r");