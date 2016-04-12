# Wolbachia infection rates - phylogenetic correction

set.seed(98523782)

setwd("~/Desktop/Projects/Wolbachia_rates")


library("phytools")

sessID <- sessionInfo()


Lep.tree <- read.nexus("Regier_data/journal.pone.0058568.s009.TRE")
Lep.tree
head(Lep.tree$tip.label)

plot.phylo(Lep.tree, cex = 0.4, no.margin = TRUE)


Lep.vcv <- vcvPhylo(tree = Lep.tree, anc.nodes = FALSE, model = "BM")
head(Lep.vcv)