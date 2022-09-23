# Jake Yeung
# Date of Creation: 2022-09-12
# File: ~/projects/sortchicAnalysis/vignettes/Supplemental_Figures/Sup_Fig_2_DE_analysis_relative_to_HSPCs_write_data_objects.R
# read RData, rename, then write to output


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(usethis)

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load data  --------------------------------------------------------------

# poisson fits
fits_all <- lapply(jmarks, function(jmark){
  inf <- file.path(indir, paste0("poisson_fit_TSS.", jmark, ".2020-12-12.newannot2.witherrors.TSS.RData"))
  load(inf, v=T)
  # return only fits
  return(jfits.lst = jfits.lst)
  # return all obs
  # return(jfits.lst = jfits.lst, dat.annots.filt.mark = dat.annots.filt.mark, ncuts.for.fit.mark = ncuts.for.fit.mark, jmat.mark = jmat.mark)
})

# gene sets
dat_genes <- GetGeneSets(inf.genes.k4me1 = "inst/extdata/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt",
                         inf.genes.k4me3 = "inst/extdata/celltype_specific_genes_defined_by_K4me3_TSS.txt")

# metadata

# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells"
# indir <- "/hpc/archive/hub_oudenaarden/userBackups/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells"
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures"
dat_metas <- lapply(jmarks, function(jmarktmp){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.shuffled.", jmarktmp, ".2021-02-19.txt")
  fread(file.path(indir, fname))
})
# cluster2col <- hash::hash(dat.metas[[1]]$cluster, dat.metas[[1]]$clustercol)


# Save output -------------------------------------------------------------

usethis::use_data(dat_metas, overwrite = FALSE)
usethis::use_data(dat_genes, overwrite = FALSE)
usethis::use_data(fits_all, overwrite = FALSE)



