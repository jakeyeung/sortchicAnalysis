# Jake Yeung
# Date of Creation: 2022-10-22
# File: ~/projects/sortchicAnalysis/vignettes/setup_data/Main_Fig_6_dual_incubation.R


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tricolore)

# Load clean metadata -----------------------------------------------------

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/processed_2022-09-12"

dat.meta.allmerged.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(indir, paste0("metadata_bonemarrow_allmerged_", jmark, ".txt"))
  fread(inf)
})

dat.meta.lsk.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(indir, paste0("metadata_bonemarrow_LSKstained_", jmark, ".txt"))
  fread(inf)
})


# LSK ---------------------------------------------------------------------

# plot the ternary plot

dat.meta.lsk.long <- lapply(dat.meta.lsk.lst, function(jdat){
  subset(jdat, select = c(cell, sca1_f_impute, ckit_f_impute, lin_f_impute, ctype.from.LL))
}) %>%
  bind_rows()

colors_and_legend <- tricolore::Tricolore(dat.meta.lsk.long,
                                          p1 = 'sca1_f_impute', p2 = 'ckit_f_impute', p3 = 'lin_f_impute',
                                          hue = 0.1, chroma = 1, lightness = 1, spread = 2)
print(colors_and_legend)


# All cell types ----------------------------------------------------------

dat.meta.colors <- subset(dat.meta.allmerged.lst$k4me1, select = c(ctype.from.LL, colcodenew))

m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.allmerged.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = colcodenew)) +
    geom_point() +
    theme_bw() +
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), legend.position = "bottom")
})

print(m.lst)

# Show PCAs ---------------------------------------------------------------

pca.out.lst <- lapply(jmarks, function(jmark){
  inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
  pca.out <- readRDS(inf.pca.out)
})

dat.pca.merge.wide.lst <- lapply(jmarks, function(jmark){
  inf.objs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_merged_arrows.", jmark, ".2022-05-05.RData")
  load(inf.objs, v=T) # dat.pca.merge.wide
  return(dat.pca.merge.wide)
})


jgrid.factor <- 2
jgrid.n <- 16
jalphadots <- 1
jalphaarrows <- 1
jsizedots <- 1
jsizearrows <- 0.9
janglearrows <- 45
jlengtharrows <- 2
# get PC factors
pc1.factors.lst <- list("k27me3" = -1,
                        "k9me3" = -1,
                        "k4me3" = 1,
                        "k4me1" = -1)
pc2.factors.lst <- list("k27me3" = 1,
                        "k9me3" = 1,
                        "k4me3" = -1,
                        "k4me1" = 1)

m.pca.arrows <- lapply(jmarks, function(jmark){

  pca.out <- pca.out.lst[[jmark]]
  dat.pca.merge.wide <- dat.pca.merge.wide.lst[[jmark]]

  varexpl1 <- round(summary(pca.out)$importance[2, 1], digits = 2)
  varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)

  pc1.factor <- pc1.factors.lst[[jmark]]
  pc2.factor <- pc2.factors.lst[[jmark]]

  dat.arrows.pca <- scchicFuncs::GetArrows(dat.pca.merge.wide, grid.n = jgrid.n, jfactor = jgrid.factor)

  dat.pca.merge.wide <- dat.pca.merge.wide %>%
    left_join(., subset(dat.meta.allmerged.lst[[jmark]], select = c(cell, colcodenew, ctype.from.LL)))

  m.pca <- ggplot(mapping = aes(x = pc1.factor * pc1.x, y = pc2.factor * pc2.x, xend = pc1.factor * pc1.y, yend = pc2.factor * pc2.y)) +
    geom_point(data = dat.pca.merge.wide, mapping = aes(color = colcodenew), alpha = jalphadots, size = jsizedots) +
    geom_segment(arrow = arrow(length=unit(jlengtharrows, "mm"), angle = janglearrows),
                 data = dat.arrows.pca, alpha = jalphaarrows, size = jsizearrows) +
    scale_color_identity() +
    theme_bw() +
    ggtitle(jmark) +
    xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
    ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.pca)

})



# Plot specific examples  -------------------------------------------------

dat.impute.lst <- lapply(jmarks[c("k4me1", "k4me3", "k27me3")], function(jmark){
  inf.impute <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects/dat_impute_tss_", jmark, ".2022-04-24.rds")
  readRDS(inf.impute)
})

indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
dat.trajs <- lapply(jmarks, function(jmark){
  inf.trajs <- file.path(indir.traj, paste0("trajs_outputs.", jmark, ".rds"))
  if (jmark == "k27me3"){
    inf.trajs <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix/trajs_outputs_batch_corrected.k27me3.2022-04-21.rds")
  }
  print(inf.trajs)
  dat.traj <- readRDS(inf.trajs)
  # remove Basophil from Granulocyte trajectory
  dat.traj$Granulocytes <- dat.traj$Granulocytes %>%
    rowwise() %>%
    mutate(is.ctype = ifelse(ctype.from.LL == "Basophils", FALSE, is.ctype))
  return(dat.traj)
})


dat.markers <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/dat_genes_sub_join_from_orig.rds")
print(unique(dat.markers$jset))
jsets <- unique(dat.markers$jset); names(jsets) <- jsets


# Filter genes  -----------------------------------------------------------

genes.keep <- dat.markers$gene

dat.impute.filt.lst <- lapply(dat.impute.lst, function(jdat){
  rows.keep.i <- which(rownames(jdat) %in% genes.keep)
  jdat[rows.keep.i, ]
})




