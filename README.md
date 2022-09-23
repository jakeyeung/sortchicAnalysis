# Repository of code used to do sortChIC processing and analysis

### processing_scripts: 

  - These are scripts used for processing data, mainly through through the high performance computing cluster

### rstudio_scripts:

  - These are interactive scripts that loads data onto memory and performs downstream analyses on an Rstudio server

### Installing dependencies 

- Some packages are install.packages() some are through bioconductor, except JFuncs is installed through github devtools::install_github("jakeyeung/JFuncs")

# Example of running pipeline:

In `example_processing_pipeline/` there is a full example that takes fastq files, demultiplexes, trims, maps, tags bam files, generates count tables, and runs count tables through dimensionality reduction. 

Example scripts are ordered from 1 to 6, which are run in ascending order. SingleCellMultiOmics scripts used v0.1.25 from https://github.com/BuysDB/SingleCellMultiOmics/releases/tag/v0.1.25 (SCMO)


## (1) demultiplex usingdemux.py (SCMO)

1-run_demux.sh: runs demux.py from SCMO. Adapters, molecule barcodes, and cell barcodes are removed from reads and encoded into the fastq headers.

## (2) remove adapters using cutadapt

2-trim_fastq.sh: runs cutadapt to remove Illumina adapters

## (3) map trimmed fastqs using bwa

3-map-fastq.sh: uses bwa to map trimmed fastqs

## (4) tag mapped bam files using bamtagmultiome.py (SCMO)

4-sort_index_tag_bam.sh: uses SCMO to read fastq headers and mapping information to record PCR duplicates, cell barcodes, cut locations, and other meta information. See https://github.com/BuysDB/SingleCellMultiOmics/blob/master/TAGS.MD for a full description of different tags and their descriptions.

## (5) generate count tables using bamToCountTable.py

5-make_count_tables*.sh: reads tagged bams and outputs counts falling into genomic regions. Uses --filterXA to ignore alternative hits, --minMQ for read quality, --dedup to remove PCR duplicates, --r1only to count fragments rather than reads (i.e. there are often two reads per fragment), -blacklist to ignore reads falling in bad regions.

## (6) run dimensionality reduction using run_LDA_model.R

6-run_LDA*.sh: runs latent Dirichlet allocation on cleaned count matrix (i.e. bad cells removed, bad regions removed).


# Vignettes for exploring data 

Check out `vignettes/explore_HSPC_dataset_quantify_log2fc_across_genesets.md` for an example of how to explore the HSPC data and access the log2 fold change estimates between cell types.


