# R package containing processed data and analysis vignettes for sortChIC

## Installation instructions

Install R package using devtools.

`devtools::install_github("jakeyeung/sortchicAnalysis", ref = "main")`

## Example pipeline for processing fastq files to generating count tables and running dimensionality reduction


In `example_processing_pipeline/` there is a full example that takes fastq files, demultiplexes, trims, maps, tags bam files, generates count tables, and runs count tables through dimensionality reduction. 

Example scripts are ordered from 1 to 6, which are run in ascending order. SingleCellMultiOmics scripts used v0.1.25 from https://github.com/BuysDB/SingleCellMultiOmics/releases/tag/v0.1.25 (SCMO)


### (1) demultiplex usingdemux.py (SCMO)

1-run_demux.sh: runs demux.py from SCMO. Adapters, molecule barcodes, and cell barcodes are removed from reads and encoded into the fastq headers.

### (2) remove adapters using cutadapt

2-trim_fastq.sh: runs cutadapt to remove Illumina adapters

### (3) map trimmed fastqs using bwa

3-map-fastq.sh: uses bwa to map trimmed fastqs

### (4) tag mapped bam files using bamtagmultiome.py (SCMO)

4-sort_index_tag_bam.sh: uses SCMO to read fastq headers and mapping information to record PCR duplicates, cell barcodes, cut locations, and other meta information. See https://github.com/BuysDB/SingleCellMultiOmics/blob/master/TAGS.MD for a full description of different tags and their descriptions.

### (5) generate count tables using bamToCountTable.py

5-make_count_tables*.sh: reads tagged bams and outputs counts falling into genomic regions. Uses --filterXA to ignore alternative hits, --minMQ for read quality, --dedup to remove PCR duplicates, --r1only to count fragments rather than reads (i.e. there are often two reads per fragment), -blacklist to ignore reads falling in bad regions.

### (6) run dimensionality reduction using run_LDA_model.R

6-run_LDA*.sh: runs latent Dirichlet allocation on cleaned count matrix (i.e. bad cells removed, bad regions removed).


## Vignettes for exploring data 

- [Explore HSPC dataset and get quantifications for log2 fold change from HSPCs to mature cell fates](https://github.com/jakeyeung/sortchicAnalysis/blob/main/vignettes/explore_HSPC_dataset_quantify_log2fc_across_genesets.md)

- [Explore dual labeling outputs and plot the relationships between heterochromatin and cell types](https://github.com/jakeyeung/sortchicAnalysis/blob/main/vignettes/show_dual_labeling_outputs.md)

- [Visualize TF activities in single cells](https://github.com/jakeyeung/sortchicAnalysis/blob/main/vignettes/plot_TF_activities.md)

- [Visualize Lin,Sca1,Kit levels on each epigenomic landscape](https://github.com/jakeyeung/sortchicAnalysis/blob/main/vignettes/plot_PCA_allmerged_LSK.md)


## FAQ (Frequenty Asked Questions)

### Where do I find the sortChIC bam files and how do I look at them at single-cell resolution?

For genome-wide analysis, it is better to look at bam files rather than count matrices (which are better for dimensionality reduction and fitting statistical models).

These bams have been already split and merged by cell types.

H3K4me3: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5018603

H3K27me3: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5018601

H3K4me1: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5018604

H3K9me3: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5018607

To get cell barcode from bam, look at the bam tag "CB", details for more tags in the github:
https://github.com/BuysDB/SingleCellMultiOmics/blob/master/TAGS.MD

To extract the cell barcode from a read from bam, use `pysam` e.g. `read.get_tag("CB")`

Avoid the qChIC bam files, they are old.
