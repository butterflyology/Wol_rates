library("rstan")
library("shinystan")
library(loo)
library(vioplot)
setwd("~/Dropbox/car_rates/Code")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# save(file = "Data/carbachia_output.RData", list = ls())
weinDat <- read.csv("Weinert_data_cleaned.csv")
car <- weinDat[weinDat$Bacterium == "Cardinium" & weinDat$Class == "insecta" & weinDat$Infected <= weinDat$Total,]
car$spp <- paste(car$Order, car$Family,car$species, sep = "_")
car$fam <- paste(car$Order, car$Family, sep="_")
car <- car[order(car$spp),]

temp <- unique.matrix(cbind(car$fam,car$spp)) 

carDat <- list(nObs = dim(car)[1],
               nSpp = length(unique(car$spp)),
               nFam = length(unique(car$fam)), 
               Pos = car$Infected, N = car$Total, 
               Species = as.numeric(as.factor(car$spp)),
               Family = as.numeric(as.factor(temp[,1])))

Pos <- car$Infected

fit <- stan(file = "hammModMultLev.stan", data = carDat, iter = 1000, chains = 4, seed = 4, pars=c('thetaFRaw','thetaSRaw','sigmaRaw'), include=FALSE) 
launch_shinystan(fit)

yNew <- t(extract(fit, pars="yNew")$yNew)
low <- apply(yNew, 1, quantile, prob=0.025)
hi <- apply(yNew, 1, quantile, prob=0.975)
med <- apply(yNew, 1, quantile, prob=0.5)

par(mfrow=c(1,1))
plot(med~Pos, type="n", ylim=c(0,max(hi)))
for(i in 1:length(med)) {
  lines(rep(Pos[i],2), c(hi[i],low[i]))
}
points(Pos, med, col="red", cex=0.8)
abline(a=0, b=1)


thetaS <- t(extract(fit, "thetaS")$thetaS)
medS <- apply(thetaS, 1, quantile, prob=0.5)
thetaS <- thetaS[order(medS),]

medS <- apply(thetaS, 1, quantile, prob=0.5)
upperS <- apply(thetaS, 1, quantile, prob=0.975)
lowerS <- apply(thetaS, 1, quantile, prob=0.025)

par(mfrow=c(1,3))
par(mar = c(3, 4, 3, 0.3))
plot(1:nrow(thetaS), medS, ylim=c(0,1), type="n", main="Species Level",ylab="theta",las=1)
for(i in 1:nrow(thetaS)) {
  lines(rep(i,2), c(upperS[i],lowerS[i]))
}
points(1:nrow(thetaS),medS, col="red", pch=".")



thetaF <- t(extract(fit, "thetaF")$thetaF)

medF <- apply(thetaF, 1, quantile, prob=0.5)
thetaF <- thetaF[order(medF),]

medF <- apply(thetaF, 1, quantile, prob=0.5)
upperF <- apply(thetaF, 1, quantile, prob=0.975)
lowerF <- apply(thetaF, 1, quantile, prob=0.025)


plot(1:length(medF), medF, ylim=c(0,1), type="n", main="Family Level",las=1, ylab="")
for(i in 1:nrow(thetaF)) {
  lines(rep(i,2), c(upperF[i],lowerF[i]))
}
points(1:nrow(thetaF),medF, col="red", pch=".")


thetaG <- extract(fit, "thetaG")$thetaG
boxplot(thetaG, main="Insecta Level", las=1, ylab="")















