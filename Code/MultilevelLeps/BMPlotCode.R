library("yarrr")
library("rethinking")
library("stargazer")

sessID <- sessionInfo()
setwd("~/Dropbox/Wol_rates")
# setwd("~/Desktop/Projects/Wolbachia_rates")
source("Code/pplots.R")
load("~/Dropbox/Wol-Leps/MultiLevModel/stanMods.R")
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
unique(wol$spp)
length(unique(wol$spp))
length(unique(wol$Family))

dim(wol[wol$Infected >= 0, ])
dim(wol[wol$Infected >= 1, ])
dim(wol[wol$Infected == 0, ])

pos <- wol[wol$Infected >= 1, ]
length(unique(pos$spp))
neg <- wol[wol$Infected == 0, ]
length(unique(neg$spp))
tot <-wol[wol$Infected >= 0, ]
length(unique(tot$Family))

temp <- unique.matrix(cbind(wol$fam, wol$spp)) 

       Pos = wol$Infected, N = wol$Total, 
              Species = datNPhy <- list(nObs = dim(wol)[1],
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
       as.numeric(as.factor(wol$spp)),
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



fitBM <- mods$fitBM

# nPhy <- as.data.frame(fitNPhy, "infectF")
BM   <- as.data.frame(fitBM, "infectF")
# OU1  <- as.data.frame(fitOU1, "infectF")
# OU5  <- as.data.frame(fitOU5, "infectF")
# OU9  <- as.data.frame(fitOU9, "infectF")

fams <- unique(wol$Family)
meds <- apply(BM,2,upper)

meds <- apply(BM, 2, upper)

medOrd <- order(meds)
ouPost <- unlist(BM[, medOrd])
families <- rep(fams[medOrd], each = nrow(BM))
out <- data.frame(prob = ouPost, families = families)

xlabel <- rep(" ", ncol(BM))


# quartz(width=7.09, height=6, bg="white")

# pdf(file = "Images/Fam_freqs.pdf", bg = "white")



par(mar = c(7, 3.4, 0.17, 0))
pPlot2(prob ~ families, data = out, ylim = c(0, 1), line.fun = median, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.8, point.o = 0.05, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", point.pch=19)

pPlot3(prob ~ families, data = out, ylim = c(0, 1), line.fun = lower, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.0, point.o = 0, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", line.lwd=2.5, add=TRUE, point.pch=19)

pPlot3(prob ~ families, data = out, ylim = c(0, 1), line.fun = upper, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", line.lwd=2.5, add=TRUE, point.pch=19)


mtext("Infection probability", 2, cex=1.11, line = 2.37, font=1)
mtext("Family", 1, cex=1.11, line = 6, font=1)
