library(rstan)
library(shinystan)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
setwd("~/Dropbox/BayesClass/Lecture6")

load("Lecture6Data.R")

datE1 <- lecture6Data$datE1

dat1 <- list(nObs = nrow(datE1), 
            x = as.numeric(relevel(datE1$mates,"oneMate"))-1,
            obs=datE1$mass)

modE1 <- stan(file="mod1.stan", dat=dat1, iter=1000, chains=4, seed=3)

print(modE1)

post <- as.data.frame(modE1, pars=c("alpha", "beta"))

oneMating <- post$alpha
multMating <- post$alpha + post$beta
boxplot(oneMating, multMating)

dif <- oneMating-multMating
hdi95 <- quantile(dif, probs=c(0.025, 0.975))
hist(dif)
abline(v=hdi95, col="red", lwd=3)

datE2 <- lecture6Data$datE2

none <- ifelse(datE2$treats == "none", 1,0)
one <- ifelse(datE2$treats == "one", 1,0)
three <- ifelse(datE2$treats == "three", 1,0)
six <- ifelse(datE2$treats == "six", 1,0)
xMat <- as.matrix(cbind(none, one, three, six))

dat2 <- list(nObs = nrow(datE2), nTreats=length(unique(datE2$treats)), xMat=xMat, obs = datE2$algalMass)

modE2 <- stan(file="mod2.stan", data=dat2, iter=1000, chains=4, seed=3)

post2 <- as.matrix(modE2)[,1:5]
means <- sweep(post2[,-1], 1, post2[,1], "+")
vals <- cbind(post2[,1], means)
colnames(vals) <- c("control","none","one","three","six")

boxplot(vals)

diffCont <- vals[,1]-vals[,2]
quant <- quantile(diffCont, probs=c(0.025, 0.975))
hist(diffCont)
abline(v=quant,col="red",lwd=3)

contGreat <- ifelse(vals[,1]>vals[,2], 1, 0)

sum(contGreat)/nrow(vals)


dat2.1 <- list(nObs = nrow(datE2),
             nTreats=length(unique(datE2$treats)), 
             treats = rep(1:5, times=as.vector(table(datE2$treats))), 
             obs = datE2$algalMass)

modE2.1 <- stan(file="mod2.1.stan", data=dat2.1, iter=1000, chains=4, seed=3)

post2.1 <- as.matrix(modE2.1)[,1:5]
colnames(post2.1) ~ unique(dat$treats)
boxplot(post2.1)

modE2.2 <- stan(file="mod2.1.stan", data=dat2.1, iter=1000, chains=4, seed=3)

post2.2 <- as.matrix(modE2.2)[,1:5]
colnames(post2.2) ~ unique(dat$treats)
boxplot(post2.2)

