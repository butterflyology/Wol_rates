<<<<<<< HEAD
set.seed(8247811)

library("rstan")
library("shinystan")
library("loo")
library("vioplot")
setwd("~/Desktop/Projects/Wolbachia_rates")

# save(list = ls(), file = "MultiLevModel/MultLevMod.RData")
# load("MultiLevModel/MultLevMod.RData")

=======
library("rstan")
library("shinystan")
library(loo)
library(vioplot)
setwd("~/Dropbox/Wol_rates/Code")
>>>>>>> e5bb42f49ed2a81e60bb53e8311d7d8bfcde9c0b
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# save(file = "Data/Wolbachia_output.RData", list = ls())
<<<<<<< HEAD
weinDat <- read.csv("Data/Weinert_data_cleaned.csv")
wol <- weinDat[weinDat$Bacterium == "Wolbachia" & weinDat$Class == "insecta" & weinDat$Infected <= weinDat$Total, ]
wol$spp <- paste(wol$Order, wol$Family, wol$species, sep = "_")
wol$fam <- paste(wol$Order, wol$Family, sep = "_")
wol <- wol[order(wol$spp), ]

temp <- unique.matrix(cbind(wol$fam, wol$spp)) 
=======
weinDat <- read.csv("Weinert_data_cleaned.csv")
wol <- weinDat[weinDat$Bacterium == "Wolbachia" & weinDat$Class == "insecta" & weinDat$Infected <= weinDat$Total,]
wol$spp <- paste(wol$Order, wol$Family,wol$species, sep = "_")
wol$fam <- paste(wol$Order, wol$Family, sep="_")
wol <- wol[order(wol$spp),]

temp <- unique.matrix(cbind(wol$fam,wol$spp)) 
>>>>>>> e5bb42f49ed2a81e60bb53e8311d7d8bfcde9c0b

wolDat <- list(nObs = dim(wol)[1],
                nSpp = length(unique(wol$spp)),
                nFam = length(unique(wol$fam)), 
                Pos = wol$Infected, N = wol$Total, 
                Species = as.numeric(as.factor(wol$spp)),
                Family = as.numeric(as.factor(temp[,1])))

Pos <- wol$Infected
<<<<<<< HEAD
fit <- stan(file = "MultiLevModel/hammModMultLev.stan", data = wolDat, iter = 1000, chains = 4, seed = 4, pars = c('thetaFRaw','thetaSRaw','sigmaRaw'), include = FALSE) 
launch_shinystan(fit)

yNew <- t(extract(fit, pars = "yNew")$yNew)
low <- apply(yNew, 1, quantile, prob = 0.025)
hi <- apply(yNew, 1, quantile, prob = 0.975)
med <- apply(yNew, 1, quantile, prob = 0.5)

# pdf(file = "MultiLevModel/Hammertime_Wol.pdf", bg = "white")
par(mfrow = c(1, 1))
plot(med ~ Pos, type = "n", ylim = c(0, max(hi)), las = 1, ylab = "Median", xlab = "Positive")
abline(a=0, b=1, lty = 2, lwd = 3, col = rgb(0, 0, 0, 0.5))
for(i in 1:length(med)) {
  lines(rep(Pos[i], 2), c(hi[i], low[i]), lwd = 2)
}
points(Pos, med, col = rgb(1, 0, 0, 0.5), cex = 0.8)


thetaS <- t(extract(fit, "thetaS")$thetaS)
medS <- apply(thetaS, 1, quantile, prob = 0.5)
thetaS <- thetaS[order(medS), ]

medS <- apply(thetaS, 1, quantile, prob = 0.5)
upperS <- apply(thetaS, 1, quantile, prob = 0.975)
lowerS <- apply(thetaS, 1, quantile, prob = 0.025)

par(mfrow=c(1, 3))
par(mar = c(3, 4, 3, 0.3))
plot(1:nrow(thetaS), medS, ylim=c(0,1), type="n", main = "Species Level", ylab = expression(paste(theta)), las = 1)
for(i in 1:nrow(thetaS)) {
  lines(rep(i, 2), c(upperS[i], lowerS[i]), col = rgb(0, 0, 0, 0.1), lwd = 0.8)
}
points(1:nrow(thetaS), medS, col = "red", pch = 23, cex = 0.3)
=======

fit <- stan(file = "hammModMultLev.stan", data = wolDat, iter = 1000, chains = 4, seed = 4, pars=c('thetaFRaw','thetaSRaw','sigmaRaw'), include=FALSE) 
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
>>>>>>> e5bb42f49ed2a81e60bb53e8311d7d8bfcde9c0b



thetaF <- t(extract(fit, "thetaF")$thetaF)

<<<<<<< HEAD
medF <- apply(thetaF, 1, quantile, prob = 0.5)
thetaF <- thetaF[order(medF), ]

medF <- apply(thetaF, 1, quantile, prob = 0.5)
upperF <- apply(thetaF, 1, quantile, prob = 0.975)
lowerF <- apply(thetaF, 1, quantile, prob = 0.025)


plot(1:length(medF), medF, ylim = c(0, 1), type = "n", main = "Family Level", las = 1, ylab = "")
for(i in 1:nrow(thetaF)) {
  lines(rep(i,2), c(upperF[i],lowerF[i]), col = rgb(0, 0, 0, 0.6), lwd = 1)
}
points(1:nrow(thetaF),medF, col = "red", pch = 23, cex = 0.3)


thetaG <- extract(fit, "thetaG")$thetaG
boxplot(thetaG, main = "Insecta Level", las = 1, ylab = "", pch = 19, col = "grey", cex = 1.5) # what does Insect level mean? The Class or Orders?
# dev.off()
# Order level? 
=======
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















>>>>>>> e5bb42f49ed2a81e60bb53e8311d7d8bfcde9c0b
