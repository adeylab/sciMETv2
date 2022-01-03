Processing sciMET datasets (including barnyard / species spike-in)

Barnyard portion
1) unidex to generate demultiplexed reads as fastq files and filter to valid barcodes (scMET mode)
2) sciMET_trim.pl - runs trimmomatic using sciMET adapters and in single-end mode
3) sciMET_align.pl to align to hybrid reference which will sort & merge as well
4) scitools rmdup - plot complexity
5) scitools barnyard-compare to get humand and mouse called cells, use read cutoff based on complexity plot

Species-specific

6) filter trimmed fastq files to be only human or mouse cell reads using scitools split-fastq
7) repeat alignment, sorting, merging, and rmdup for species alignments
8) run sciMET_extract.pl on the rmdup & filtered bam file to create context 'chrom' folders then sort it with sciMET_sortChroms.pl
9) run CHmeth2mtx.pl using windows and CGmeth2mtx.pl using annotated regions with the chroms folder
10) scitools irlba on the CH methylaiton matrix (not the ratio.mtx)
11) remove the first dimension (head -n1 then tail -n49 to same file) and see if it improves (otherwise keep all 50)
12) run umap / matrix-pg / etc... & plot using scitools

13) run getGeneMeth.pl using CH over gene bodies or CG over promoter regions
    uses the chroms folder from step 8
    use the annot file for clusters to make aggregate methylaiton for clusters
    will also produce single-cell level methylation
