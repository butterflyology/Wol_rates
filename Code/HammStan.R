# Data and code for infection prevalence and frequency paper by Marion and Hamm


library("rstan")
library("shinystan")

setwd("~/Desktop/Projects/Wolbachia_rates")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load("Data/Wolbachia_rates.Rdata")

fit <- stan(file = "Wol_Stan/HammMod.stan", data = W.data, iter = 10000, chains = 4, seed = 4) 
launch_shinystan(fit)

Ntot <- aggregate(N ~ Species, FUN = sum)$N
Nstud <- as.vector(table(Species))

sppInfect <- monitor(extract(fit, pars = "infectNew", permuted = FALSE), print = FALSE)[, c("mean")]

hist(sppInfect,breaks=20)
plot(density(sppInfect))
boxplot(sppInfect)
mean(sppInfect);median(sppInfect)
sd(sppInfect)


#####
##### data from Weinert et al 2015
#####

