# Data and code for infection prevalence and frequency paper by Marion and Hamm
library("rstan")
library("shinystan")
library("loo")
library("phytools")
setwd("~/Dropbox/Wol_rates")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# save(file = "Data/Wolbachia_output.RData", list = ls())
#load("Data/Lep.vcv.R")
load("Data/finalTree.R")

weinDat <- read.csv("Data/Weinert_data_cleaned.csv",stringsAsFactors=FALSE)


sessID <- sessionInfo()

wol <- weinDat[weinDat$Bacterium == "Wolbachia" 
               & weinDat$Order == "lepidoptera" 
               & weinDat$Infected <= weinDat$Total
               & weinDat$Family != "Unid." 
               & weinDat$Family != "Unid" 
               & weinDat$Family != "Riodinidae",]

wol$spp <- paste(wol$Family,wol$species, sep = "_")
wol$fam <- wol$Family
wol <- wol[order(wol$spp),]

140

temp <- unique.matrix(cbind(wol$fam,wol$spp)) 

wolDat <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[,1])))


wolDatNPhy <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = diag(1, nrow=length(unique(wol$fam))),
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[,1])))


fit <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = wolDat, iter = 500, chains = 4, seed = 4, control = list(adapt_delta = 0.96))


fitNP <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = wolDatNPhy, iter = 2500, chains = 4, seed = 4, control = list(adapt_delta = 0.96))

launch_shinystan(fit)

#fitSimp <- stan(file = "MultiLevModel/hammModMultLev.stan", data = wolDat, iter = 2500, chains = 4, seed = 4, control = list(adapt_delta = 0.96)) 

print(fitSimp, "thetaG")
print(fit, "infectG")


logPhy <- extract_log_lik(fit)
logSimp <- extract_log_lik(fit2)

aicPhy <- waic(logPhy)
aicSimp <- waic(logSimp)

compare(aicPhy, aicSimp)

#----------------------
## Phylogeny Plot
#----------------------

probS <- as.matrix(fit, "infectS")
medS <- apply(probS, 2, median)
probS <- probS[oS,]
oS <- order(medS)
plot(1:ncol(probS), medS[oS], ylim = c(0,1),type="n")
for(s in 1:ncol(probS)) {
  lines(rep(s,nrow(probS)), probS[,s], pch=16, col=adjustcolor("grey20", alpha.f=0.1))
}
points(1:length(medS), as.vector(medS[oS]), pch=16, col="red",cex=0.4)






