# Data and code for infection prevalence and frequency paper by Marion and Hamm
library("rstan")
library("shinystan")
library("loo")
library("yarrr")
library("rethinking")

sessID <- sessionInfo()
 setwd("~/Dropbox/Wol_rates")
# setwd("~/Desktop/Projects/Wolbachia_rates")
source("Code/pplots.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# save(file = "Data/Wolbachia_output.RData", list = ls())
# load("Data/Lep.vcv.R")
load("Data/Lep.vcv.ultra.R")
# source("~/Dropbox/Wol-Leps/MultiLevModel/stanMods.R")

weinDat <- read.csv("Data/Weinert_data_cleaned.csv", stringsAsFactors = FALSE)



wol <- weinDat[weinDat$Bacterium == "Wolbachia" 
               & weinDat$Order == "lepidoptera" 
               & weinDat$Infected <= weinDat$Total
               & weinDat$Family != "Unid." 
               & weinDat$Family != "Unid" 
               & weinDat$Family != "Riodinidae", ]

wol$spp <- paste(wol$Family, wol$species, sep = "_")
wol$fam <- wol$Family
wol <- wol[order(wol$spp), ]
str(wol)
sum(wol$Total)
sum(wol$Infected)

temp <- unique.matrix(cbind(wol$fam, wol$spp)) 

datNPhy <- list(nObs = dim(wol)[1],
                nSpp = length(unique(wol$spp)),
                nFam = length(unique(wol$fam)),
                phyCor = diag(1, nrow=length(unique(wol$fam))),
                Pos = wol$Infected, N = wol$Total, 
                Species = as.numeric(as.factor(wol$spp)),
                Family = as.numeric(as.factor(temp[, 1])))


datBM <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv$vcv.BM,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[, 1])))

datOU1 <- list(nObs = dim(wol)[1],
              nSpp = length(unique(wol$spp)),
              nFam = length(unique(wol$fam)),
              phyCor = Lep.vcv$vcv.ou1,
              Pos = wol$Infected, N = wol$Total, 
              Species = as.numeric(as.factor(wol$spp)),
              Family = as.numeric(as.factor(temp[, 1])))

datOU5 <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv$vcv.ou5,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[, 1])))

datOU9 <- list(nObs = dim(wol)[1],
               nSpp = length(unique(wol$spp)),
               nFam = length(unique(wol$fam)),
               phyCor = Lep.vcv$vcv.ou9,
               Pos = wol$Infected, N = wol$Total, 
               Species = as.numeric(as.factor(wol$spp)),
               Family = as.numeric(as.factor(temp[, 1])))

fitNPhy <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datNPhy, iter = 5000, chains = 4, seed = 4, control = list(adapt_delta = 0.96), pars=c("zFam", "zSpp", "sppLogit", "yNew", "infectS"), include = FALSE, thin=5)

fitBM <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datBM, iter = 5000, chains = 4, seed = 4, control = list(adapt_delta = 0.96), pars=c("zFam", "zSpp", "sppLogit", "yNew", "infectS"), include = FALSE, thin=5)

fitOU1 <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datOU1, iter = 5000, chains = 4, seed = 4, control = list(adapt_delta = 0.96), pars=c("zFam", "zSpp", "sppLogit", "yNew", "infectS"), include = FALSE, thin=5)

fitOU5 <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datOU5, iter = 5000, chains = 4, seed = 4, control = list(adapt_delta = 0.96), pars=c("zFam", "zSpp", "sppLogit", "yNew", "infectS"), include = FALSE, thin=5)

fitOU9 <- stan(file = "Code/MultilevelLeps/famLevPhylo.stan", data = datOU9, iter = 5000, chains = 4, seed = 4, control = list(adapt_delta = 0.96), pars=c("zFam", "zSpp", "sppLogit", "yNew", "infectS"), include = FALSE, thin=5)

# mods <- list(fitNPhy=fitNPhy, fitBM=fitBM, fitOU1=fitOU1, fitOU5=fitOU5, fitOU9 = fitOU9)
# save(mods, file="Code/MultilevelLeps/stanMods.R")
load("~/Dropbox/Wol-Leps/MultiLevModel/stanMods.R")
 fitNPhy <- mods$fitNPhy
 fitBM <- mods$fitBM
 fitOU1 <- mods$fitOU1
 fitOU5 <- mods$fitOU5
 fitOU9 <- mods$fitOU9
# 
# logNP <- extract_log_lik(fitNPhy)
# logBM <- extract_log_lik(fitBM)
# logOU1 <- extract_log_lik(fitOU1)
# logOU5 <- extract_log_lik(fitOU5)
# logOU9 <- extract_log_lik(fitOU9)
# 
# aicNP <- waic(logNP)
# aicBM <- waic(logBM)
# aicOU1 <- waic(logOU1)
# aicOU5 <- waic(logOU5)
# aicOU9 <- waic(logOU9)
# compare(aicNP, aicBM, aicOU1, aicOU5, aicOU9)

x <- compare(fitNPhy, fitBM, fitOU1, fitOU5, fitOU9)
wts <- data.frame(mod = dimnames(x@output)[[1]], wts = round(x@output$weight,3))

# unaveraged
nPhy <- as.data.frame(fitNPhy, "infectG")$infectG
BM   <- as.data.frame(fitBM, "infectG")$infectG
OU1  <- as.data.frame(fitOU1, "infectG")$infectG
OU5  <- as.data.frame(fitOU5, "infectG")$infectG
OU9  <- as.data.frame(fitOU9, "infectG")$infectG

avg <- rowSums(cbind(
  as.matrix(fitNPhy, "infectG") * wts[wts$mod == "fitNPhy",2], 
  as.matrix(fitBM, "infectG") * wts[wts$mod == "fitBM",2],  
  as.matrix(fitOU1, "infectG") * wts[wts$mod == "fitOU1",2],
  as.matrix(fitOU5, "infectG") * wts[wts$mod == "fitOU5",2],
  as.matrix(fitOU9, "infectG") * wts[wts$mod == "fitOU9",2]))

upper <- function(x) {
  return(quantile(x, prob = 0.975))
}

lower <- function(x) {
  return(quantile(x, prob = 0.025))
}


prob <- c(nPhy, BM, OU1, OU5, OU9, avg)
mod <- rep(c("NP","BM", "0.1", "0.5","0.9", " Model"), each = length(BM))
out <- data.frame(prob = prob, mod = mod)

#quartz(width=3.4, height=4, bg="white")

par(mar = c(2.9, 3.4, 0.16, 0))
par(xpd=TRUE)
pPlot1(prob ~ mod, data = out, ylim = c(0, 1), line.fun = median, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.8, point.o = 0.03, yaxt = "l", bty = "l", las = 1, ylab="", col.lab="white", point.pch=16)


pPlot3(prob ~ mod, data = out, ylim = c(0, 1), line.fun = lower, line.lwd=2.5, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, yaxt = "l", bty = "l", las = 1, add = TRUE, ylab="" )

pPlot3(prob ~ mod, data = out, ylim = c(0, 1), line.fun = upper, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, yaxt = "l", bty = "l", las = 1, add = TRUE, line.lwd=2.5)

mtext("Infection probability", 2, cex=1.11, line = 2.37, font=1)
segments(x0 = 2.7, y0 = -0.15, x1 = 5.3, y1=-0.15, lwd=3)
text(x=3.99, y=-0.189, expression(paste("OU ", bold(alpha))), cex=1.1, font=2)

text(x=6.1, y=-0.138, "\nAvg.", cex=1, font=1)

quartz.save(type = "pdf",file = "infection probability.pdf", bg="white")


nPhy <- as.data.frame(fitNPhy, "infectF")
BM   <- as.data.frame(fitBM, "infectF")
OU1  <- as.data.frame(fitOU1, "infectF")
OU5  <- as.data.frame(fitOU5, "infectF")
OU9  <- as.data.frame(fitOU9, "infectF")

famList <- list(
  as.matrix(fitNPhy, "infectF") * wts[wts$mod == "fitNPhy",2], 
  as.matrix(fitBM, "infectF") * wts[wts$mod == "fitBM",2],  
  as.matrix(fitOU1, "infectF") * wts[wts$mod == "fitOU1",2],
  as.matrix(fitOU5, "infectF") * wts[wts$mod == "fitOU5",2],
  as.matrix(fitOU9, "infectF") * wts[wts$mod == "fitOU9",2])

avgF <- as.data.frame(Reduce('+', famList))


fams <- unique(wol$Family)
meds <- apply(avgF,2,upper)

meds <- apply(avgF, 2, upper)

medOrd <- order(meds)
ouPost <- unlist(avgF[, medOrd])
families <- rep(fams[medOrd], each = nrow(avgF))
out <- data.frame(prob = ouPost, families = families)

xlabel <- rep(" ", ncol(OU9))


quartz(width=7.09, height=6, bg="white")
par(mar = c(6.6, 3.4, 0.17, 0))
par(xpd=FALSE)

# pdf(file = "Images/Fam_freqs.pdf", bg = "white")
par(las = 2)
pirateplot(prob ~ families, data = out, ylim = c(0, 1), line.fun = upper, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.8, point.o = 0.05, yaxt = "l", bty = "n", las = 2, xaxt = "n", xlab = "Lepidoptera family", ylab = "Infection frequency")
           
pirateplot(prob ~ families, data = out, ylim = c(0, 1), line.fun = lower, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, add = TRUE)
           
pirateplot(prob ~ families, data = out, ylim = c(0, 1), line.fun = median, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, add = TRUE)
# dev.off()


pPlot2(prob ~ families, data = out, ylim = c(0, 1), line.fun = median, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.8, point.o = 0.05, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", point.pch=16)
    
pPlot3(prob ~ families, data = out, ylim = c(0, 1), line.fun = lower, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.0, point.o = 0, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", line.lwd=2.5, add=TRUE)

pPlot3(prob ~ families, data = out, ylim = c(0, 1), line.fun = upper, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", line.lwd=2.5, add=TRUE)
       

mtext("Infection probability", 2, cex=1.11, line = 2.37, font=1)
mtext("Family", 1, cex=1.11, line = 5.5, font=1)

quartz.save(type = "pdf", file="familyProb.pdf")


#----------------------
## Phylogeny Plot
#----------------------

probS <- as.matrix(fitOU1, "infectS")
medS <- apply(probS, 2, median)
probS <- probS[oS, ]
oS <- order(medS)
plot(1:ncol(probS), medS[oS], ylim = c(0, 1),type = "n")
for(s in 1:ncol(probS)) {
  lines(rep(s, nrow(probS)), probS[, s], pch = 16, col = adjustcolor("grey20", alpha.f = 0.1))
}
points(1:length(medS), as.vector(medS[oS]), pch = 16, col = "red", cex = 0.4)


# barplot of sample size by family
unique(wol$spp)
length(unique(wol$spp))

g <- table(wol$Family)
# g <- as.list(g)
length(g)

pdf("Fams_bar.pdf", bg = "white", width=3.39, height=4)

quartz(width=3.4, height=4, bg="white")
#par(pin = c(3.39, 3.5))
par(ps=10)
par(mar = c(6.1, 3, 0.4, 0))

barplot(sort(g), las = 2, ylim = c(0, 350), cex.names = 0.9, ylab = "Samples", axes=FALSE, space=0.3)
 #dev.off()



# pdf("Images/Fams_bar.pdf", bg = "white")
par(mar = c(6.1, 4, 1, 1))
barplot(sort(g), las = 2, ylim = c(0, 350), cex.names = 0.9, ylab = "Samples")
# dev.off()


# correcting for species diversity by dividing by number of species based on van Nieukerken et al. 2001
gprime <- as.table(c(g[1] / 185, g[2] / 49, g[3] / 113, g[4] / 9655, g[6] / 660, g[7] / 20, g[8] / 24569, g[9] / 339, g[10] / 4700, g[11] / 23002, g[12] / 1866, g[13] / 36, g[14] / 606, g[15] / 1952, g[16] / 5201, g[17] / 11772, g[18] / 3800, g[19] / 6152, g[20] / 570, g[21] / 1164, g[22] / 150, g[23] / 1318, g[24] / 5921, g[25] / 2349, g[26] / 1463, g[27] / 10387, g[28] / 686))

# pdf("Images/Prop_fams_bar.pdf", bg = "white")
par(mar = c(6.1, 4.5, 1, 1))

barplot(sort(gprime), las = 2, cex.names = 0.9, ylab = "Proportion of family represented", ylim = c(0, 0.7))
# dev.off()


barplot(sort(gprime), las = 2, cex.names = 0.9, ylab = "Proportion of family represented", ylim = c(0, 1.0))
# dev.off()

