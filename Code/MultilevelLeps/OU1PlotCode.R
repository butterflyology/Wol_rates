library("yarrr")
library("rethinking")
library("stargazer")

sessID <- sessionInfo()
setwd("~/Dropbox/Wol_rates")
# setwd("~/Desktop/Projects/Wolbachia_rates")
source("Code/pplots.R")
load("~/Dropbox/Wol-Leps/MultiLevModel/stanMods.R")

fitOU1 <- mods$fitOU1


# nPhy <- as.data.frame(fitNPhy, "infectF")
# BM   <- as.data.frame(fitBM, "infectF")
OU1  <- as.data.frame(fitOU1, "infectF")
# OU5  <- as.data.frame(fitOU5, "infectF")
# OU9  <- as.data.frame(fitOU9, "infectF")

fams <- unique(wol$Family)
meds <- apply(OU1,2,upper)

meds <- apply(OU1, 2, upper)

medOrd <- order(meds)
ouPost <- unlist(OU1[, medOrd])
families <- rep(fams[medOrd], each = nrow(OU1))
out <- data.frame(prob = ouPost, families = families)

xlabel <- rep(" ", ncol(OU1))


# quartz(width=7.09, height=6, bg="white")
par(mar = c(7, 3.4, 0.17, 0))
pPlot2(prob ~ families, data = out, ylim = c(0, 1), line.fun = median, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.8, point.o = 0.05, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", point.pch=16)

pPlot3(prob ~ families, data = out, ylim = c(0, 1), line.fun = lower, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0.0, point.o = 0, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", line.lwd=2.5, add=TRUE)

pPlot3(prob ~ families, data = out, ylim = c(0, 1), line.fun = upper, pal = "ghostbusters", bar.o = 0, line.o = 0.7, bean.o = 0, point.o = 0, yaxt = "l", bty = "l", las = 1, xaxt = "n", xlab = "", ylab="", line.lwd=2.5, add=TRUE)


mtext("Infection probability", 2, cex=1.11, line = 2.37, font=1)
mtext("Family", 1, cex=1.11, line = 6, font=1)