# Wolbachia infection rates - phylogenetic correction

set.seed(98523782)

### setwd("~/Desktop/Projects/Wolbachia_rates")
### setwd("~/Dropbox/Wol_rates")

library("phytools")
library("geiger")

sessID <- sessionInfo()


Lep.tree <- read.nexus("Regier_data/Fam_codes.TRE")
Lep.tree
head(Lep.tree$tip.label)


#plot.phylo(Lep.tree, cex = 0.4, no.margin = TRUE)



# Lep.vcv <- vcvPhylo(tree = Lep.tree, anc.nodes = FALSE, model = "BM")

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


lam <- 10^(-1:4)
cv <- sapply(lam, function(x) sum(attr(chronopl(Lep.nodups, lambda = x, CV = TRUE), "D2")))
plot(lam, cv, pch = 19, ylab = "cross-validation score", xlab = expression(paste(lambda)), las = 1, cex = 1.5) # lowest CV. Small lambda means that every branch gets its own rate, larger and you are more clock like. Best lambda is 1e3.

Lep.nodups.ultra <- chronopl(phy = Lep.nodups, lambda = 1e3, CV = TRUE, eval.max = 1e3, iter.max = 1e4)
is.ultrametric(Lep.nodups.ultra)
plot(Lep.nodups.ultra)

# write.nexus(Lep.nodups.ultra, file = "Regier_data/Lep.ultra.nex")
 Lep.ultra <- read.nexus("Regier_data/Lep.ultra.nex")


Lep.vcv <- vcv.phylo(phy = Lep.ultra, corr = TRUE, model = "BM")
oPhy <- order(colnames(Lep.vcv))
Lep.vcv <- Lep.vcv[oPhy, oPhy]

x <- famCode %in% colnames(Lep.vcv) != TRUE
famCode[x]


weinDat <- read.csv("Data/Weinert_data_cleaned.csv", stringsAsFactors = FALSE)

wol <- weinDat[weinDat$Bacterium == "Wolbachia" & weinDat$Order == "lepidoptera" & weinDat$Infected <= weinDat$Total & weinDat$Family != "Unid." & weinDat$Family != "Unid" & weinDat$Family != "Riodinidae",]

wol$spp <- paste(wol$Family,wol$species, sep = "_")
wol$fam <- wol$Family
wol <- wol[order(wol$spp),]
famCode <- casefold(substring(unique(wol$fam), 1,4), upper=TRUE)


unShared <- which(Lep.ultra$tip.label %in% famCode == FALSE)

finalTree <- drop.tip(Lep.ultra, tip = unShared)
 write.nexus(finalTree, file="finalTree.tre")

read.nexus("Data/finalTree.tre")


vcv.BM <- vcv.phylo(phy = finalTree, corr = TRUE, model = "BM")
oPhy <- order(colnames(vcv.BM))
vcv.BM <- vcv.BM[oPhy, oPhy]



names(wol)
head(Lep.vcv)


## rescale via OU model with different alphas
ou.1 <- rescale(x = finalTree, model = "OU", 0.1)
vcv.ou1 <- vcv.phylo(phy = ou.1, corr = TRUE, model = "BM")[oPhy,oPhy]

ou.5 <- rescale(x = finalTree, model = "OU", 0.5)
vcv.ou5 <- vcv.phylo(phy = ou.5, corr = TRUE, model = "BM")[oPhy,oPhy]

ou.9 <- rescale(x = finalTree, model = "OU", 0.9)
vcv.ou9 <- vcv.phylo(phy = ou.9, corr = TRUE, model = "BM")[oPhy,oPhy]

Lep.vcv <- list(vcv.BM=vcv.BM, vcv.ou1=vcv.ou1, vcv.ou5=vcv.ou5, vcv.ou9=vcv.ou9)
save(Lep.vcv, file="Data/Lep.vcv.ultra.R")


##### Need to simulate data


