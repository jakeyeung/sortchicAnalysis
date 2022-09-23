# Jake Yeung
# Date of Creation: 2022-09-12
# File: ~/projects/sortchicAnalysis/vignettes/Supplemental_Figures/Sup_Fig_4_DE_for_H3K9me3.R
# recreates /hpc/archive/hub_oudenaarden/userBackups/jyeung/data/scChiC/from_rstudioserver/pdfs_all/H3K4me1_H3K9me3_differential_expression_outputs/other_marks/H3K9me3_bins_by_pval.gene_gene_correlations.H3K9me3_vs_other_marks_bins_labeled.LowInK9.TRUE.2020-12-20.bugfix.boxplots.pdf
# from script:
# /home/hub_oudenaarden/jyeung/projects/scChiC/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/10-DE_downstream_all_marks.K9me3_bins.R



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)
library(DescTools)



keeptop <- 150
low.in.k9 <- TRUE

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures/SupFig4"

# Get DE outputs ----------------------------------------------------------

# fits_lst_lst <- lapply(jmarks, function(jmark){
fits_lst_lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jinf <- file.path(hubprefix, paste0("poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.RData"))
  load(jinf, v=T)
  return(jfits.lst)
})



ctypes <- c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs")
ctypes.k9me3 <- c("Eryths", "Bcells", "Granulocytes", "HSPCs")
dat.metas <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("umaps_final_get_marks.", jmark, ".metadata.2020-12-17.txt"))
  fread(inf)
})

dat.metas <- lapply(jmarks, function(jmark){
  dat.metas.tmp <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x == "rep1old", "zold", "anew"))
  } else {
    dat.metas.tmp$jrep2 <- sapply(dat.metas.tmp$jrep, function(x) ifelse(x != "rep1old", "zold", "anew"))
  }
  return(dat.metas.tmp)
})


# jmetas.pretty.lst <- lapply(jmarks, function(jmark){
metas_pretty_lst <- lapply(jmarks, function(jmark){
  jmeta <- dat.metas[[jmark]]
  if (jmark != "H3K9me3"){
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes)
  } else {
    jmeta$cluster <- factor(jmeta$cluster, levels = ctypes.k9me3)
  }
  jmeta <- jmeta %>% arrange(cluster, jrep)
})


cells.keep.lst <- lapply(metas_pretty_lst, function(jdat){
  jdat$cell
})

# Get bins keep -----------------------------------------------------------




pvals.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- fits_lst_lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
})

pvals.lst2 <- lapply(fits_lst_lst$H3K9me3, function(x) x$pval)
k9.bins <- which(pvals.lst2 < 1e-10)


# Get k9 bins and plot  ---------------------------------------------------


params.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfits.lst <- fits_lst_lst[[jmark]]
  jnames <- names(jfits.lst); names(jnames) <- jnames
  params.dat.all <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xkeep <- grepl("^Cluster.*.Estimate$", x = names(x))
    jparams <- x[xkeep]
    data.frame(bin = jname, param = names(jparams), estimate = unlist(jparams), mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  if (jmark == "H3K9me3"){
    params.dat.all <- params.dat.all %>%
      mutate(param = gsub("Eryth", "Eryths", param),
             param = gsub("Lymphoid", "Bcells", param))
  }
  return(params.dat.all)
})



pval.k9.sub <- subset(pvals.lst$H3K9me3, pval < 1e-10) %>%
  arrange(desc(pval))

k9.bins.names <- unique(pval.k9.sub$bin)
ctypes.keep <- c("Eryths", "Bcells", "Granulocytes")
params.keep <- paste("Cluster", ctypes.keep, ".Estimate", sep = "")

params_dat_wide <- data.table::dcast(subset(params.lst$H3K9me3, bin %in% k9.bins.names), formula = bin ~ param, value.var = "estimate") %>%
  rowwise() %>%
  mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
         ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
         ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
         ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
         ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
         ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
         Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
         Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
         Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
         HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))



params.dat.wide.lst <- lapply(jmarks, function(jmark){
  jsub <- subset(params.lst[[jmark]], bin %in% k9.bins.names & param %in% params.keep) %>%
    group_by(bin) %>% filter(max(abs(estimate)) < 5)
  jdat <- data.table::dcast(jsub,
                            formula = bin ~ param, value.var = "estimate") %>%
    rowwise() %>%
    mutate(ClusterBcells.Estimate = ClusterBcells.Estimate / log(2),
           ClusterGranulocytes.Estimate = ClusterGranulocytes.Estimate / log(2),
           ClusterEryths.Estimate = ClusterEryths.Estimate / log(2),
           ClusterBcells.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterBcells.Estimate, ClusterGranulocytes.Estimate)),
           ClusterEryths.Estimate_ClusterGranulocytes.Estimate = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate)),
           ClusterBcells.Estimate_ClusterEryths.Estimate =  mean(c(ClusterBcells.Estimate, ClusterEryths.Estimate)),
           Bcells.effect = ClusterBcells.Estimate - ClusterEryths.Estimate_ClusterGranulocytes.Estimate,
           Eryths.effect = ClusterEryths.Estimate - ClusterBcells.Estimate_ClusterGranulocytes.Estimate,
           Granulocytes.effect = ClusterGranulocytes.Estimate - ClusterBcells.Estimate_ClusterEryths.Estimate,
           HSPCs.effect = mean(c(ClusterEryths.Estimate, ClusterGranulocytes.Estimate, ClusterBcells.Estimate)))
  # keep only effect cnames
  cnames.keep.i <- grep("effect$", colnames(jdat))
  cnames.new <- paste(colnames(jdat)[cnames.keep.i], jmark, sep = "_")
  colnames(jdat)[cnames.keep.i] <- cnames.new
  cnames.keep.bin.i <- grep("bin", colnames(jdat))
  cnames.keep.merged.i <- c(cnames.keep.bin.i, cnames.keep.i)
  jdat.filt <- jdat[, cnames.keep.merged.i]
  return(jdat.filt)
})



if (low.in.k9){
  jsort.hspcs <- params_dat_wide %>%
    group_by(bin) %>%
    # arrange(HSPCs.effect)
    arrange(desc(HSPCs.effect))
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]

  jsort.bcell <- params_dat_wide %>%
    group_by(bin) %>%
    # arrange(desc(Bcells.effect))
    arrange(Bcells.effect)
  jbins.bcell <- jsort.bcell$bin[1:keeptop]

  jsort.granu <- params_dat_wide %>%
    group_by(bin) %>%
    # arrange(desc(Granulocytes.effect))
    arrange(Granulocytes.effect)
  jbins.granu <- jsort.granu$bin[1:keeptop]

  jsort.eryth <- params_dat_wide %>%
    group_by(bin) %>%
    # arrange(descEryths.effect))
    arrange(Eryths.effect)
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
} else {
  jsort.hspcs <- params_dat_wide %>%
    group_by(bin) %>%
    arrange(HSPCs.effect)
  jbins.hspcs <- jsort.hspcs$bin[1:keeptop]

  jsort.bcell <- params_dat_wide %>%
    group_by(bin) %>%
    arrange(desc(Bcells.effect))
  jbins.bcell <- jsort.bcell$bin[1:keeptop]

  jsort.granu <- params_dat_wide %>%
    group_by(bin) %>%
    arrange(desc(Granulocytes.effect))
  jbins.granu <- jsort.granu$bin[1:keeptop]

  jsort.eryth <- params_dat_wide %>%
    group_by(bin) %>%
    arrange(desc(Eryths.effect))
  jbins.eryth <- jsort.eryth$bin[1:keeptop]
}


bins.keep <- c(jbins.eryth, jbins.bcell, jbins.granu, jbins.hspcs)
bins.keep.lst <- list("Eryths" = jbins.eryth,
                      "Bcells" = jbins.bcell,
                      "Granulocytes" = jbins.granu,
                      "HSPCs" = jbins.hspcs)
bnames <- names(bins.keep.lst); names(bnames) <- bnames

bins_keep <- bins.keep
bins_keep_lst <- bins.keep.lst


# Load RData  -------------------------------------------------------------

inf.mat.adj <- file.path(hubprefix, "batch_corrected_imputed_values.bins.all_marks.mat.namesfix.2020-12-20.H3K27me3rep2rep3reseq.RData")
load(inf.mat.adj, v=T)

mat.adj.lst <- lapply(mat.adj.lst, function(jmat){
  rownames(jmat) <- jmat$rname
  jmat$rname <- NULL
  return(jmat)
})

bins.common <- Reduce(f = intersect, x = lapply(mat.adj.lst, rownames))
bins.keep.common <- bins.keep[bins.keep %in% bins.common]

mat_lst <- lapply(jmarks, function(jmark){
  jmat <- mat.adj.lst[[jmark]][bins.keep.common, cells.keep.lst[[jmark]]]
  print(dim(jmat))
  jmat <- apply(jmat, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.99)))
  jmat <- t(apply(jmat, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.99))))
  return(jmat)
})



# Write outputs -----------------------------------------------------------

params_dat_wide_lst <- params.dat.wide.lst
params_lst <- params.lst
pvals_lst <- pvals.lst
# pvals_lst2 <- pvals.lst2

# usethis::use_data(pvals_lst2, overwrite=FALSE)
usethis::use_data(pvals_lst, overwrite=FALSE)
usethis::use_data(params_lst, overwrite=FALSE)
usethis::use_data(params_dat_wide_lst, overwrite=FALSE)
usethis::use_data(params_dat_wide, overwrite = FALSE)
usethis::use_data(bins_keep, overwrite=FALSE)
usethis::use_data(bins_keep_lst, overwrite=FALSE)
usethis::use_data(mat_lst, overwrite = FALSE)
usethis::use_data(metas_pretty_lst, overwrite = FALSE)



