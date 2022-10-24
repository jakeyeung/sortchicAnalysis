# Jake Yeung
# Date of Creation: 2022-10-24
# File: ~/projects/sortchicAnalysis/vignettes/plot_TF_activities.R
#

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(heatmap3)

# Load motif activities ---------------------------------------------------

data("mat_motif_activity_k4me1")
data("dat_metas")

dat.merge <- as.data.frame(dat_metas$H3K4me1)

rownames(dat.merge) <- dat.merge$cell

cells.ordered <- rownames(mat_motif_activity_k4me1)
dat.merge <- dat.merge[cells.ordered, ]

dat.col <- dat.merge %>%
  dplyr::select(c(cell, clustercol))
colvec <- dat.col$clustercol

jmeth <- "ward.D2"
hm.out.transpose <- heatmap3(t(mat_motif_activity_k4me1), margins = c(5, 8), cexCol = 0.35, Colv = NA, Rowv = TRUE,
                             ColSideColors = colvec,
                             ColSideLabs = "celltype",
                             labCol = FALSE, scale = "row", revC = FALSE,
                             distfun = dist, hclustfun = hclust, method = jmeth)



# Show individual motifs on UMAP  -----------------------------------------


# examples of motifs
jmotif <- "Cebpb"
jmotif <- "Ebf1"
jmotif <- "Tal1"
jmotif <- "Gata1"
jmotif <- "Irf1"
jmotif <- "Hoxb5"
jmotif <- "Runx1"

jmotif <- "Erg"
jtitle <- paste("Motif activity", jmotif)
dat.merge.motifs <- data.frame(cell = rownames(mat_motif_activity_k4me1), motif_activity = mat_motif_activity_k4me1[, jmotif], stringsAsFactors = FALSE) %>%
  left_join(dat.merge)

m <- ggplot(dat.merge.motifs, aes_string(x = "umap1", y = "umap2", color = "motif_activity")) +
  geom_point(size = 0.75) +
  ggtitle(jtitle) +
  theme_bw() +
  scale_color_viridis_c() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(m)

