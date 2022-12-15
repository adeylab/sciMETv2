#!/usr/bin/perl

use Getopt::Std; %opt = ();

getopts("w:h:O:A:G:HS:", \%opt);

$width = 5;
$height = 4;

$die = "

sciMET_plotGeneByAnnot.pl (options) geneByAnnot.meth.txt

Will make a dotplot of targets x annotations from the file
produced by sciMET_getGeneMeth.pl

If 1 or more gene/target names are provided it will plot only those.

Options:
   -w   [INT]   Width of plot (def = $width)
   -h   [INT]   Height of plot (def = $height)
   -O   [STR]   Output prefix (def = input file prefix)
   -H           CH methylation mode (colors / scaling; def = CG)
   -S   [MIN,MAX] Min and max bounds on scale (def = all)
   
   -G   [LST]   Comma-separated list of genes
   -A   [LST]   Comma-separated list of annotations
                 For -G and -A, will plot in specified order
   
";

if (!defined $ARGV[0]) {die $die};

if (defined $opt{'w'}) {$width = $opt{'w'}};
if (defined $opt{'h'}) {$height = $opt{'h'}};
if (defined $opt{'S'}) {($scaleMin,$scaleMax) = split(/,/, $opt{'S'})};

if (!defined $opt{'O'}) {
	$opt{'O'} = $ARGV[0];
	$opt{'O'} =~ s/\.txt$//;
}

if (defined $opt{'G'} || defined $opt{'A'} || defined $opt{'S'}) {
	$included = 0;
	
	if (defined $opt{'G'}) {
		@GENE_LIST = split(/,/, $opt{'G'});
		foreach $gene (@GENE_LIST) {$GENES{$gene} = 1};
	}
	
	if (defined $opt{'A'}) {
		@ANNOT_LIST = split(/,/, $opt{'A'});
		foreach $annot (@ANNOT_LIST) {$ANNOTS{$annot} = 1};
	}
	
	open OUT, ">$opt{'O'}.plot.txt";
	open IN, "$ARGV[0]";
	while ($l = <IN>) {
		chomp $l;
		($target,$annot,$meth,$zscore) = split(/\t/, $l);
		
		if (defined $GENES{$target} || !defined $opt{'G'}) {
			if (defined $ANNOTS{$annot} || !defined $opt{'A'}) {
				if (defined $opt{'S'}) {
					if ($meth < $scaleMin) {$meth = $scaleMin};
					if ($meth > $scaleMax) {$meth = $scaleMax};
					$l = "$target\t$annot\t$meth\t$zscore";
				}
				$included++;
				print OUT "$l\n";
			}
		}
	} close IN; close OUT;
	$infile = "$opt{'O'}.plot.txt";
	
	if ($included<1) {
		system("rm -f $opt{'O'}.plot.txt");
		die "ERROR: Genes listed not detected in gene x annot methylation file!";
	}
} else {
	$infile = $ARGV[0];
}

if (defined $opt{'G'}) {
	$order_genes = "IN\$V1 <- factor(IN\$V1, levels = c(";
	for ($i = (@GENE_LIST-1); $i >= 0; $i--) {
		$order_genes .= "\"$GENE_LIST[$i]\",";
	} $order_genes =~ s/,$//; $order_genes .= "))";
} else {$order_genes = ""};

if (defined $opt{'A'}) {
	$order_clust = "IN\$V2 <- factor(IN\$V2, levels = c(";
	for ($i = 0; $i < @ANNOT_LIST; $i++) {
		$order_clust .= "\"$ANNOT_LIST[$i]\",";
	} $order_clust =~ s/,$//; $order_clust .= "))";
} else {$order_clust = ""};

if (!defined $opt{'H'}) { #CG colors / scaling
open R, ">$opt{'O'}.plot.r";
print R "library(ggplot2)
IN<-read.table(\"$infile\")
$order_clust
$order_genes
PLT<-ggplot(data=IN) + theme_bw() +
	 geom_point(aes(V2,V1,size=V3,color=V4)) +
	 scale_color_gradient2(low=\"#ff0082\", mid=\"#bdbdbd\", high=\"#08519c\") +
	 xlab(\"Cluster\") + ylab(\"Gene\") + labs(color=\"z-score\",size=\"mCG\") +
	 theme(axis.text.x = element_text(angle=90,hjust=0))
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=$width,height=$height)
"; close R;
} else { #CH colors / scaling
open R, ">$opt{'O'}.plot.r";
print R "library(ggplot2)
IN<-read.table(\"$infile\")
$order_clust
$order_genes
PLT<-ggplot(data=IN) + theme_bw() +
	 geom_point(aes(V2,V1,size=V3,color=V4)) +
	 scale_color_gradient2(low=\"#49e679\", mid=\"#bdbdbd\", high=\"#08519c\") +
	 xlab(\"Cluster\") + ylab(\"Gene\") + labs(color=\"z-score\",size=\"mCH\") +
	 theme(axis.text.x = element_text(angle=90,hjust=0))
ggsave(plot=PLT,filename=\"$opt{'O'}.png\",width=$width,height=$height)
ggsave(plot=PLT,filename=\"$opt{'O'}.pdf\",width=$width,height=$height)
"; close R;
}

system("Rscript $opt{'O'}.plot.r");