# Data and code for infection prevalence and frequency paper by Marion and Hamm
library("rstan")
library("shinystan")
library("loo")
library("yarrr")
setwd("~/Dropbox/Wol_rates")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# save(file = "Data/Wolbachia_output.RData", list = ls())
load("Data/Lep.vcv.R")


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



temp <- unique.matrix(cbind(wol$fam,wol$spp)) 

datNPhy <- list(nObs = dim(wol)[1],
                nSpp = length(unique(wol$spp)),
                nFam = length(unique(wol$fam)),
                phyCor = diag(1, nrow=length(unique(wol$fam))),
                Pos = wol$Infected, N = wol$Total, 
                Species = as.numeric(as.factor(wol$spp)),
                Family = as.numeric(as.factor(temp[,1])))


datBM <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv$vcv.BM,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[,1])))

datOU1 <- list(nObs = dim(wol)[1],
              nSpp = length(unique(wol$spp)),
              nFam = length(unique(wol$fam)),
              phyCor = Lep.vcv$vcv.ou1,
              Pos = wol$Infected, N = wol$Total, 
              Species = as.numeric(as.factor(wol$spp)),
              Family = as.numeric(as.factor(temp[,1])))

datOU5 <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv$vcv.ou5,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[,1])))

datOU9 <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv$vcv.ou9,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[,1])))

fitNPhy <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datNPhy, iter = 1000, chains = 4, seed = 4, control = list(adapt_delta = 0.96))

fitBM <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datBM, iter = 1000, chains = 4, seed = 4, control = list(adapt_delta = 0.96), pars=c("zFam", "zSpp"), include=FALSE)

fitOU1 <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datOU1, iter = 1000, chains = 4, seed = 4, control = list(adapt_delta = 0.96))

fitOU5 <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datOU5, iter = 1000, chains = 4, seed = 4, control = list(adapt_delta = 0.96))

fitOU9 <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datOU9, iter = 1000, chains = 4, seed = 4, control = list(adapt_delta = 0.96))


logNP <- extract_log_lik(fitNPhy)
logBM <- extract_log_lik(fitBM)
logOU1 <- extract_log_lik(fitOU1)
logOU5 <- extract_log_lik(fitOU5)
logOU9 <- extract_log_lik(fitOU9)

aicNP <- waic(logNP)
aicBM <- waic(logBM)
aicOU1 <- waic(logOU1)
aicOU5 <- waic(logOU5)
aicOU9 <- waic(logOU9)

x <- compare(aicNP, aicBM, aicOU1, aicOU5, aicOU9)


nPhy <- as.data.frame(fitNPhy, "infectG")
BM <- as.data.frame(fitBM, "infectG")
OU1 <- as.data.frame(fitOU1, "infectG")
OU5 <- as.data.frame(fitOU5, "infectG")
OU9 <- as.data.frame(fitOU9, "infectG")

prob <- rbind(nPhy, BM, OU1,OU5,OU9)$infectG
mod <- rep(c("nP", "BM", "OU1", "OU5", "OU9"), each=nrow(BM))
out <- data.frame(prob=prob,mod=mod)

upper <- function(x) {
  return(quantile(x,prob=0.975))
}

lower <- function(x) {
  return(quantile(x,prob=0.025))
}


par(mar = c(3,4,1.5,0.8))
pirateplot(prob~mod, data=out, ylim=c(0,1), line.fun=upper, pal="bugs", bar.o=0, line.o=0.7, bean.o=0.8, point.o=0.05, yaxt="l",bty="l", las=1)

pirateplot(prob~mod, data=out, ylim=c(0,1), line.fun=lower, pal="bugs", bar.o=0, line.o=0.7, bean.o=0, point.o=0, yaxt="l",bty="l", las=1, add=TRUE)

pirateplot(prob~mod, data=out, ylim=c(0,1), line.fun=median, pal="bugs", bar.o=0, line.o=0.7, bean.o=0, point.o=0, yaxt="l",bty="l", las=1, add=TRUE)



OU5 <- as.data.frame(fitOU5, "infectF")

meds <- apply(OU5,2,upper)
medOrd <- order(meds)
ouPost <- unlist(OU5[,medOrd])
families <- rep(colnames(Lep.vcv$vcv.ou9)[medOrd], each=nrow(OU5))
out <- data.frame(prob=ouPost, families=families)

pirateplot(prob~families, data=out, ylim=c(0,1), line.fun=upper, pal="ghostbusters", bar.o=0, line.o=0.7, bean.o=0.8, point.o=0.05, yaxt="l",bty="l", las=1)
           
pirateplot(prob~families, data=out, ylim=c(0,1), line.fun=lower, pal="ghostbusters", bar.o=0, line.o=0.7, bean.o=0, point.o=0, yaxt="l",bty="l", las=1, add=TRUE)
           
pirateplot(prob~families, data=out, ylim=c(0,1), line.fun=median, pal="ghostbusters", bar.o=0, line.o=0.7, bean.o=0, point.o=0, yaxt="l",bty="l", las=1, add=TRUE)







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





