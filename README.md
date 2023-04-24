Processing sciMET datasets
1) unidex to generate demultiplexed reads as fastq files and filter to valid barcodes (scMET mode)
3) sciMET_trim.pl - runs trim_galore wrapper
4) sciMET_align_BSBolt.pl - BSBolt paired-end alignment wrapper, outputs a name (cellID) sorted bam file
5) sciMET_rmdup_pe.pl - de-duplicates by barcode and outputs a name sorted deduplicated bam file
6) sciMET_BSBolt_extract.pl - extracts mC status in CG or CH context, outputs in bismark-style calls for each chromosome
7) sciMET_sortChroms.pl - wrapper to sort the bed files in the output and gzip them (then delete the non-gzipped bed files)
8) sciMET_meth2mtx.pl - generate matrixes for mC fraction, ratio, coverage, and score (-1 to 1 centered on 0 as global mC of the cell)
9) dimensionality reduction using irlba (R) (on ratio or score matrix), then UMAP and clustering on the irlba output