# Jake Yeung
# Date of Creation: 2022-10-24
# File: ~/projects/sortchicAnalysis/vignettes/plot_PCA_allmerged_LSK.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(sortchicAnalysis)
library(tricolore)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

data("dat_meta_lsk_lst")
data("dat_meta_allmerged_lst")
data("dat_pca_merge_wide_lst") # dat_pca_merge_wide_lst and dat_var_expl_lst

# Show UMAPs with LSK -----------------------------------------------------


dat.meta.lsk.long <- lapply(dat_meta_lsk_lst, function(jdat){
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
  jsub <- dat_meta_lsk_lst[[jmark]]
  colors_and_legend <- tricolore::Tricolore(jsub, 'sca1_f_impute', 'ckit_f_impute', 'lin_f_impute')
  jsub$rgb <- colors_and_legend$rgb
  # print(colors_and_legend$key)

  # add other cells as grey underneath
  m <- ggplot() +
    geom_point(mapping = aes(x = umap1, umap2), data = dat_meta_allmerged_lst[[jmark]], color = "grey85") +
    geom_point(mapping = aes(x = umap1, umap2, color = rgb), data = jsub) +
    scale_color_identity() +
    ggtitle(jmark) +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  return(m)
})
print(m.lsk.lst)


# Show UMAPs with cell type -----------------------------------------------

dat.meta.colors <- subset(dat_meta_allmerged_lst$k4me1, select = c(ctype.from.LL, colcodenew))
m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat_meta_allmerged_lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = colcodenew)) +
    geom_point() +
    theme_bw() +
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") +
    ggtitle(jmark) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank(), legend.position = "bottom")
})
print(m.lst)



# Show PCA with velocities ------------------------------------------------



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

  dat.pca.merge.wide <- dat_pca_merge_wide_lst[[jmark]]

  varexpl1 <- subset(dat_var_expl_lst[[jmark]], PC == "PC1")$varexplvec[[1]]
  varexpl2 <- subset(dat_var_expl_lst[[jmark]], PC == "PC2")$varexplvec[[1]]

  pc1.factor <- pc1.factors.lst[[jmark]]
  pc2.factor <- pc2.factors.lst[[jmark]]

  dat.arrows.pca <- scchicFuncs::GetArrows(dat.pca.merge.wide, grid.n = jgrid.n, jfactor = jgrid.factor)

  dat.pca.merge.wide <- dat.pca.merge.wide %>%
    left_join(., subset(dat_meta_allmerged_lst[[jmark]], select = c(cell, colcodenew)))

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

})

print(m.pca.arrows)





