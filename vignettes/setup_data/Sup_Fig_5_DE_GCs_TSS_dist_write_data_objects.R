# Jake Yeung
# Date of Creation: 2022-09-16
# File: ~/projects/sortchicAnalysis/vignettes/setup_data/Sup_Fig_5_DE_GCs_TSS_dist_write_data_objects.R
# /home/hub_oudenaarden/jyeung/projects/scChiC/scripts/rstudioserver_analysis/spikeins/dendrograms/4-annotate_high_bins_gc.R
# recreates
# /hpc/archive/hub_oudenaarden/userBackups/jyeung/data/scChiC/from_rstudioserver/pdfs_all/dendrograms_heatmaps/High_bins_gc_and_distance_to_gene.pdf



library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(usethis)


# Load high bins ----------------------------------------------------------

# indir.bins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/DE_downstream_analysis_BM_allmerged_H3K27me3_cleaned"
indir.bins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures/SupFig5"
dat.high.bins <- lapply(jmarks, function(jmark){
  fname <- paste0("High_bins_all_marks_padjcutoff_dists_to_TSS.annot_table.jmidbug_fixed.", jmark, ".2021-02-19.txt")
  inf <- file.path(indir.bins, fname)
  fread(inf)
})

# bins.lst <- lapply(dat.high.bins, function(jdat){
#   jdat$CoordOriginal
# })
dat_high_bins <- dat.high.bins
usethis::use_data(dat_high_bins)

load("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures/SupFig5/gcs_genomewide.RData", v=T)
gr_gc_dat_dedup <- gr.gc.dat[!duplicated(gr.gc.dat), ]

usethis::use_data(gr_gc_dat_dedup)
