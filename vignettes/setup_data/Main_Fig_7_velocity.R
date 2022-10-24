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

# all marks together
colors_and_legend <- tricolore::Tricolore(dat.meta.lsk.long,
                                          p1 = 'sca1_f_impute', p2 = 'ckit_f_impute', p3 = 'lin_f_impute',
                                          hue = 0.1, chroma = 1, lightness = 1, spread = 2)
print(colors_and_legend$key)

# split by marks
m.lsk.lst <- lapply(jmarks, function(jmark){
  jsub <- dat.meta.lsk.lst[[jmark]]
  colors_and_legend <- tricolore::Tricolore(jsub, 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
  jsub$rgb <- colors_and_legend$rgb
  print(colors_and_legend$key)

  # add other cells as grey underneath
  m <- ggplot() +
    geom_point(mapping = aes(x = umap1, umap2), data = dat.meta.lst[[jmark]], color = "grey85") +
    geom_point(mapping = aes(x = umap1, umap2, color = rgb), data = jsub) +
    scale_color_identity() +
    ggtitle(jmark) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})


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
  # inf.pca.out <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/downstream_pca_out.", jmark, ".2022-05-05.rds")
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
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") +
    theme_bw() +
    ggtitle(jmark) +
    xlab(paste0("PC1 (", varexpl1 * 100, "%)")) +
    ylab(paste0("PC2 (", varexpl2 * 100, "%)")) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.pca)

})


# Get PCA explained -------------------------------------------------------

dat_var_expl_lst <- lapply(jmarks, function(jmark){
  pca.out <- pca.out.lst[[jmark]]
  varexplvec <- round(summary(pca.out)$importance[2, 1:10], digits = 2)
  # varexpl2 <- round(summary(pca.out)$importance[2, 2], digits = 2)
  dat.var.expl <- data.frame(varexplvec, stringsAsFactors = FALSE) %>%
    mutate(mark = jmark)
  dat.var.expl$PC <- rownames(dat.var.expl)
  return(dat.var.expl)
})


# Save outputs ------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/primetime_data_to_remake_figures/MainFig7_outputs"

pca_out_lst <- pca.out.lst
dat_pca_merge_wide_lst <- dat.pca.merge.wide.lst
dat_meta_colors <- dat.meta.colors

dat_meta_allmerged_lst <- dat.meta.allmerged.lst
dat_meta_lsk_lst <- dat.meta.lsk.lst

# outpca <- file.path(outdir, paste0("pca_out_lst.rda"))
outpcamerged <- file.path(outdir, paste0("dat_pca_merge_wide_lst.rda"))
outallmerged <- file.path(outdir, paste0("dat_meta_allmerged_lst.rda"))
outallmergedlsk <- file.path(outdir, paste0("dat_meta_lsk_lst.rda"))

# save(pca_out_lst, file = outpca)
save(dat_var_expl_lst, dat_pca_merge_wide_lst, dat_meta_colors, file = outpcamerged)
save(dat_meta_allmerged_lst, file = outallmerged)
save(dat_meta_lsk_lst, file = outallmergedlsk)





