# Jake Yeung
# Date of Creation: 2022-09-16
# File: ~/projects/sortchicAnalysis/vignettes/Supplemental_Figures/Sup_Fig_5_bins_TSS_and_GC.R
#

rm(list=ls())

library(sortchicAnalysis)


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

library(hash)
library(topicmodels)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

data(dat_high_bins)
data(gr_gc_dat_dedup)
data(metas_pretty_lst)


make.plots <- TRUE

binsize <- 50000

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures/SupFig5_plots"
outpdf <- file.path(outdir, paste0("High_bins_gc_and_distance_to_gene.", Sys.Date(), ".pdf"))

if (make.plots){
  pdf(file = outpdf, useDingbats = FALSE)
}

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# hubprefix <- "/home/jyeung/hub_oudenaarden"

# merge some celltypes
ctypes <- list("Eryths" = "Erythroid",
               "Bcells" = "Lymphoid",
               "NKs" = "Lymphoid",
               "Granulocytes" = "Myeloid",
               "Basophils" = "Myeloid",
               "pDCs" = "Lymphoid",
               "DCs" = "Myeloid",
               "HSPCs" = "HSPCs",
               "Erythroid" = "Erythroid",
               "Lymphoid" = "Lymphoid",
               "Myeloid" = "Myeloid")


merge.ctypes.by.lineage <- FALSE

dat.metas <- lapply(metas_pretty_lst, function(x){
  x %>%
    rowwise() %>%
    mutate(lineage = ctypes[[cluster]])
})


# cluster to col
cluster2col <- hash(dat.metas$H3K4me1$cluster, dat.metas$H3K4me1$clustercol)
cluster2col[["Erythroid"]] <- "#0072B2"
cluster2col[["Lymphoid"]] <- "#56B4E9"
cluster2col[["Myeloid"]] <- "#D55E00"



# Load bins ---------------------------------------------------------------


bins.lst <- lapply(dat_high_bins, function(jdat){
  jdat$CoordOriginal
})

gr.gc.dat <- gr_gc_dat_dedup
gr.gc.dat.filt.lst <- lapply(jmarks, function(jmark){
  jbins <- bins.lst[[jmark]]
  jout <- subset(gr.gc.dat, bname %in% jbins) %>%
    mutate(mark = jmark)
  jout <- jout[!duplicated(jout$bname), ]
})
gr.gc.dat.filt.long <- gr.gc.dat.filt.lst %>%
  bind_rows()
gr.gc.dat.filt.long$mark <- factor(gr.gc.dat.filt.long$mark, levels = jmarks)



nbins <- length(unique(gr.gc.dat.filt.long$bname))

ggplot(gr.gc.dat.filt.long, aes(x = mark, y = gc)) +
  geom_boxplot(outlier.alpha = 0.1) +
  ggtitle(paste("nbins:", nbins)) +
  theme_bw()  +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# distance to nearest gene
dat.high.bins.long <- dat_high_bins %>%
  bind_rows()
dat.high.bins.long$mark <- factor(dat.high.bins.long$mark, levels = jmarks)

nbins2 <- length(unique(dat.high.bins.long$CoordOriginal))
ggplot(dat.high.bins.long, aes(x = mark, y = abs(distanceToTSS) + 1)) +
  geom_boxplot(outlier.alpha = 0.1) +
  scale_y_log10() +
  ggtitle(paste("nbins:", nbins2)) +
  geom_hline(yintercept = 25000, linetype = "dotted") +
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dev.off()



