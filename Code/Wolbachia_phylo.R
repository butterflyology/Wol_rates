# Wolbachia infection rates - phylogenetic correction

set.seed(98523782)

### setwd("~/Desktop/Projects/Wolbachia_rates")
setwd("~/Dropbox/Wol_rates")

library("phytools")

sessID <- sessionInfo()


Lep.tree <- read.nexus("Regier_data/Fam_codes.TRE")
Lep.tree
head(Lep.tree$tip.label)


#plot.phylo(Lep.tree, cex = 0.4, no.margin = TRUE)


#Lep.vcv <- vcvPhylo(tree = Lep.tree, anc.nodes = FALSE, model = "BM")

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

###### Get Weinert data ready
weinDat <- read.csv("Data/Weinert_data_cleaned.csv", stringsAsFactors = FALSE)

wol <- weinDat[weinDat$Bacterium == "Wolbachia" & weinDat$Order == "lepidoptera" & weinDat$Infected <= weinDat$Total & weinDat$Family != "Unid." & weinDat$Family != "Unid" & weinDat$Family != "Riodinidae",]

wol$spp <- paste(wol$Family,wol$species, sep = "_")
wol$fam <- wol$Family
wol <- wol[order(wol$spp),]
famCode <- casefold(substring(unique(wol$fam), 1,4), upper=TRUE)

unShared <- which(Lep.nodups$tip.label %in% famCode == FALSE)

finalTree <- drop.tip(Lep.nodups, tip = unShared)
save(finalTree, file="finalTree.R")

Lep.vcv <- vcv.phylo(phy = finalTree, corr = TRUE, model = "BM")
oPhy <- order(colnames(Lep.vcv))
Lep.vcv <- Lep.vcv[oPhy, oPhy]

save(Lep.vcv, file="Lep.vcv.R")