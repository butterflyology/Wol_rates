# Wolbachia infection rates - phylogenetic correction

set.seed(98523782)

setwd("~/Desktop/Projects/Wolbachia_rates")


library("phytools")

sessID <- sessionInfo()


Lep.tree <- read.nexus("Regier_data/Fam_codes.TRE")
Lep.tree
head(Lep.tree$tip.label)

length(Lep.tree$tip.label)
length(unique(Lep.tree$tip.label))

duplicated(Lep.tree$tip.label)
Lep.tree$tip.label[duplicated(Lep.tree$tip.label)]
get.off <- which(duplicated(Lep.tree$tip.label))

Lep.nodups <- drop.tip(Lep.tree, tip = get.off)
length(Lep.nodups$tip.label)

# pdf(file = "Regier_fams.pdf", bg = "white")
plot.phylo(Lep.nodups, cex = 0.5, no.margin = TRUE)
# dev.off()
duplicated(Lep.nodups$tip.label)

Lep.vcv <- vcvPhylo(tree = Lep.nodups, anc.nodes = FALSE, model = "BM")
head(Lep.vcv)