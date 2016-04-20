# Wolbachia infection rates - phylogenetic correction

set.seed(98523782)

setwd("~/Desktop/Projects/Wolbachia_rates")


library("phytools")

sessID <- sessionInfo()


Lep.tree <- read.nexus("Regier_data/Fam_codes.TRE")
Lep.tree
head(Lep.tree$tip.label)


plot.phylo(Lep.tree, cex = 0.4, no.margin = TRUE)


Lep.vcv <- vcvPhylo(tree = Lep.tree, anc.nodes = FALSE, model = "BM")

length(Lep.tree$tip.label)
length(unique(Lep.tree$tip.label))

duplicated(Lep.tree$tip.label)
Lep.tree$tip.label[duplicated(Lep.tree$tip.label)]
get.off <- which(duplicated(Lep.tree$tip.label))

Lep.nodups <- drop.tip(Lep.tree, tip = get.off)
length(Lep.nodups$tip.label)
plot.phylo(Lep.nodups, cex = 0.4, no.margin = TRUE)
duplicated(Lep.nodups$tip.label)

Lep.vcv <- vcvPhylo(tree = Lep.nodups, anc.nodes = FALSE, model = "BM")

head(Lep.vcv)