Processing sciMET datasets (including barnyard / species spike-in)

Barnyard portion
1) unidex to generate demultiplexed reads as fastq files and filter to valid barcodes (scMET mode); this trims the first 10bp of
     read 1 which includes the randomer ligation region.
3) sciMET_trim.pl - runs trimmomatic using sciMET adapters and in single-end mode.
4) sciMET_align.pl to align to hybrid reference which will sort & merge as well
     if files from other runs should be merged, run with -m to skip the merge and then run samtools merge with all name sorted bam files
5) Remove duplicates with sciMET_rmdup.pl - then scitools plot-complexity to assess complexity and cell read depth
6) scitools barnyard-compare to get humand and mouse called cells, use read cutoff based on complexity plot

Species-specific

6) filter trimmed fastq files to be only human or mouse cell reads using sciMET_speciesSplit.pl
7) repeat alignment, sorting, merging, and rmdup for species alignments
8) run sciMET_extract.pl on the rmdup & filtered bam file to create context 'chrom' folders then sort it with sciMET_sortChroms.pl
9) run sciMET_meth2mtx.pl using windows the chroms folder for CG and then CH (separate runs, will auto-detect)
10) filter the CH matrix based on coverage etc... using sciMET_filtMtx.pl
11) scitools irlba on the CH methylaiton matrix (not the ratio.mtx)
12) run umap / matrix-pg / etc... & plot using scitools

13) run getGeneMeth.pl using CH over gene bodies or CG over promoter regions uses the chroms folder from step 8 use the annot
    file for clusters to make aggregate methylaiton for clusters will also produce single-cell level methylation
    this can be run with a long list of genes (use -L option) to have a file that can be plotted from
	plot with plotGeneByAnnot.pl and it will plot all (if a smalls et was used for making the file) or a subset that can be provided

Merging datasets:
1) CH methylaiton matrixes can be merged using sciMET_mergeMtx.pl. Merge the cov and the methylation level matrixes so that
   you can filter the merged matrix after.
2) Tools that leverage the sorted chroms folder containg methylaitonc calls for cells can use a list of folders, including
   sciMET_meth2mtx.pl and getGeneMeth.pl. As long as cell names are different, it will read them all in. This can be done
   instead of merging matrixes etc.