# Wolbachia junk folder


betabinexch0 <- function(theta, data){
	eta = theta[, 1]
	K = theta[, 2]
	y = data[, 1]
	n = data[, 2]
	N = length(y)
	val = O * K
	for(i in 1:N){
		val = val + lbeta(K * eta + y[i], K * (1 - eta) + n[i] -y[i])
		val = val - N * lbeta(K * eta, K * (1 - eta))
		val = val - 2 * log(1 + K) - log(eta) - log(1 - eta)
		retrun(val)
		}
}

data(cancermortality)

mycontour(betabinexch0, c(0.0001, 0.003, 1, 20000), cancermortality)


betabinom <- function(theta, data){
    eta <-  theta[1]
    K <-  theta[2]
    y <-  data[, 1]
    n <-  data[, 2]
    N <-  length(y)
    logf <-  function(y, n, K, eta){
    	lbeta(K * eta + y, K * (1 - eta) + n - y) - lbeta(K * eta, K * (1 - eta))
    	}
    val <-  sum(logf(y, n, K, eta))
    val <-  val - 2 * log(1 + K) - log(eta) - log(1 - eta)
    return(val)
}

mycontour(betabinom, c(0.0001, 0.003, 1, 20000), cancermortality)


mycontour(betabinexch0, c(0.0001, 0.003, 1, 20000), cancermortality, xlab = expression(paste(eta)), ylab = expression(paste(kappa))) #same if using #betabinom function I modified

mycontour(betabinexch, c(-8, -4.5, 3, 16.5), cancermortality, xlab = expression(paste("logit ", eta)), ylab = expression(paste("log ", kappa)))


# does this assume a common probability of infection?
mycontour(betabinexch0, c(0.2, 0.35, 0.4, 1.4), Wol.rates, xlab = expression(paste(eta)), ylab = expression(paste(kappa)), las = 1, main = "Beta-binomial model")

mycontour(betabinexch, c(-1.5, -0.5, -0.7, 0.5), Wol.rates, xlab = expression(paste("logit ", eta)), ylab = expression(paste("log ", kappa)), las = 1, main = "Transformed parameters")

#####
##### Bayes influence 
#####
# observation sensitivity in beta-binomial model
data(cancermortality)

start <- array(c(-7, 6), c(1, 2))
fit <- laplace(betabinexch, start, cancermortality)
tpar <- list(m = fit$mode, var = 2 * fit$var, df = 4)
theta <- sir(betabinexch, tpar, 1000, cancermortality)
intervals <- bayes.influence(theta, cancermortality)


#####
##### BayesPref
#####
library("bayespref")

head(WD)
summary(WD)
str(WD)

# Don't run this. It is very similar to the JAGS model but generates a giant data file (100 Mb)
# wb1 <- bayesPref(WD, mcmcL = 1e4, pops = FALSE)
str(wb1)
p1 <- prefPlot(wb1[[1]])
cw <- credibleIntervals(wb1[[1]], burn = 1e3)


burn.low <- floor(0.1 * (length(wb1[[1]]$PopPref[1, ])))
post.burn <- wb1[[1]]$PopPref[1, ][burn.low:length(wb1[[1]]$PopPref[1, ])]
str(post.burn)


effectiveSize(post.burn) # 237
autocorr.plot(post.burn, lag.max = 200)



Pos <- W.stuff$y
N <- W.stuff$n
Species <- W.stuff$Species
 
 
 
W.data <- list(Nsp = length(levels(WD.sp$Species)),  #Give your the number of species
 Nobs = length(WD.sp$Positive),  #Give you the numbers of observations
 Species = as.numeric(WD.sp$Species), #Give you the species ID in numeric form. 
 N = WD.sp$n, 
 Pos = WD.sp$Positive)
 
Mod.2 <- function(){
#prior
   sigma ~ dunif(0, 10)
   tau <- pow(sigma, -2)
   for (j in 1:Nsp){
    alpha[j] ~ dnorm(0, tau)
   }
 #likelihood
   for(i in 1:Nobs){
     Pos[i] ~ dbinom(theta, N[i])
     logit(theta[]) <- alpha[Species[i]]  
   }
} 
 
my.params <- c("alpha", "sigma")
my.inits <- function(){list(alpha = runif(W.data$Nsp, -2, 2), sigma = runif(1, 0.01, 5))}
 
 
 
 ##and you can run the new model
fit.1 <- jags(data = W.data, inits = my.inits , my.params, n.chains = 4,  n.iter = 1e4, n.burnin = 1e3, model.file = Mod.2)

##### Model 1, this mode is not correct


Pos <- W.data$Pos
N <- W.data$N
Species <- length(W.data$Pos)

Mod.1 <- function(){
	#prior
	theta ~ dbeta(1, 1)
	#likelihood
	for(i in 1:Species){
		Pos[i] ~ dbinom(theta, N[i])
	}
}

params <- c("theta")
W.data <- list(Pos = Pos, N = N, Species = Species)

fit.1 <- jags(data = W.data, inits = NULL, params, n.chains = 4, n.iter = 1e4, n.burnin = 1e3, model.file = Mod.1)
fit.1
str(fit.1)

fit.1.mcmc <- as.mcmc(fit.1)
summary(fit.1.mcmc)
xyplot(fit.1.mcmc)
autocorr.plot(fit.1.mcmc) # Nice!
effectiveSize(fit.1.mcmc) # 4e3, very nice

hist(fit.1$BUGSoutput$sims.list$theta, col = "grey", las = 1, xlab = expression(paste(theta)), xlim = c(0.30, 0.35), breaks = 20, main = "")
median(fit.1$BUGSoutput$sims.list$theta)
quantile(fit.1$BUGSoutput$sims.list$theta, probs = c(0.025, 0.975)) # This is not right. Identical to bn.conf.exact()

Linden <- function(){
	y[i] ~ dbinom(theta[i], N[i])
	logit(theta[i]) <- alpha0 + alpha[i]
	alpha[i] ~ dnorm(0, tau)
}
Lin.params <- c("theta", "alpha")
Lind.data <- list(Pos = Pos, N = N)

Lin.fit <- jags(data = Lind.data, inits = NULL, Lin.params, n.chains = 4, n.iter = 1e4, n.burnin = 1e3, model.file = Linden)
fit.1
str(fit.1)



Mod.1 <- function(){
#prior
  for (j in 1:Nsp){
   alpha[j] ~ dbeta(1, 1)
  }
#likelihood
  for(i in 1:Nobs){
    Pos[i] ~ dbinom(theta[i], N[i])
    theta[i] <- alpha[Species[i]]
   }
} 

params <- c("theta", "alpha")

fit.1 <- jags(data = W.data, inits = NULL, params, n.chains = 4, n.iter = 1e4, n.burnin = 1e3, model.file = Mod.1)

str(fit.1)


hist(fit.1$BUGSoutput$mean$theta, xlim = c(0, 1), breaks = 20)
hist(fit.1$BUGSoutput$mean$alpha, xlim = c(0, 1), breaks = 20)


# Beta-binomial model for overdispersion applied to Wolbachia infection rates

library("rjags")
library("R2jags")
library("coda")
library("lattice")
library("rstan")

set.seed(87635378)
setwd("~/Desktop/Projects/Wolbachia_rates/")
source("DBDA2E-utilities.R")
# load("Wolbachia_rates.RData")
# save.image("Wolbachia_rates.RData")
(SID <- sessionInfo())


#####
##### Importing and setting up data
#####

##### Import data set by species
WD.sp <- read.delim("Rates3.txt", sep = "\t", header = TRUE)
WD.sp <- WD.sp[complete.cases(WD.sp), ]
WD.sp$index <- as.factor(as.numeric(WD.sp$Species))
str(WD.sp)
summary(WD.sp)
head(WD.sp)

WD.sp[WD.sp$Species == "Acraea_encedon", ]

##### Import data set by Genus
WD.gn <- read.delim("Rates2.txt", sep = "\t", header = TRUE)
WD.gn <- WD.gn[complete.cases(WD.gn), ]
WD.gn$species <- NULL
WD.gn$index <- as.factor(as.integer(WD.gn$Genus))
str(WD.gn)
head(WD.gn)

#####
##### Descriptive statistics
#####

##### making a proportion parameter: logit(p) = log(p / (1 - p))
prop.inf <- WD.sp[, 3] / WD.sp[, 2]
hist(prop.inf, col = "grey", breaks = 20, ylim = c(0, 500), las = 1, xlab = "Proportion infected / sample")

logit.prop <- log(prop.inf / (1 - prop.inf))
hist(logit.prop, col = "grey", xlab = "logit proportion infected", xlim = c(-4.1, 4.1), las = 1)

##### Exact confidence interval for binomial data

bn.conf.exact <- function(x, n, conf.level){
  p <- x / n
	alpha <- 1 - conf.level
	alpha <- rep(alpha, length = length(p))
	alpha2 <- 0.5 * alpha
	x1 <- x == 0
	x2 <- x == n
	lb <- ub <- x
	lb[x1] <- 1
	ub[x2] <- n[x2] - 1
	lcl <- 1 - qbeta(1 - alpha2, n + 1 - x, lb)
	ucl <- 1 - qbeta(alpha2, n - ub, x + 1)
	if(any(x1))
	lcl[x1] <- rep(0, sum[x1])
	if(any(x2))
	ucl[x2] <- rep(1, sum[x1])
	out <- matrix(c(p, lcl, ucl), nrow = 1, ncol = 3, dimnames=list("data", c("mean", "lower", "upper")))
	print(out)
}

bn.conf.exact(1, 2, conf.level = 0.95) #show the addition
bn.conf.exact(100, 200, conf.level = 0.95)


colSums(WD.sp[, 2:4])
bn.conf.exact(sum(WD.sp[, 3]), sum(WD.sp[, 2]), conf.level = 0.95)

bn.conf.exact(sum(WD.sp[WD.sp$Species == "Acraea_encedon", 3]), sum(WD.sp[WD.sp$Species == "Acraea_encedon", 2]), conf.level = 0.95)

#####
##### hierarchical Bayesian fun time
#####

Mod.1 <- function(){
#prior
  for (j in 1:Nsp){
   alpha[j] ~ dbeta(1, 1)
  }
#likelihood
  for(i in 1:Nobs){
    Pos[i] ~ dbinom(theta[i], N[i])
    theta[i] <- alpha[Species[i]]
   }
} 

Pos <- WD.sp$Positive
N <- WD.sp$n
Nobs <- length(WD.sp$Positive)
Nsp <- length(unique(WD.sp$index))
Species <- as.numeric(WD.sp$index) # Need an index to associate Nobs with Nsp, this does not work

params <- c("theta", "alpha")
W.data <- list(Pos = Pos, Nobs = Nobs, Nsp = Nsp, Species = Species, N = N)

fit.1 <- jags(data = W.data, inits = NULL, params, n.chains = 4, n.iter = 1e4, n.burnin = 1e3, model.file = Mod.1)

str(fit.1)


hist(fit.1$BUGSoutput$mean$theta, xlim = c(0, 1), breaks = 20)
hist(fit.1$BUGSoutput$mean$alpha, xlim = c(0, 1), breaks = 20)


library("shinystan")

modelString = "
	data {
	int<lower=0> Nobs;
	int<lower=0> Nsp;
	int<lower=0> Pos[Nobs];
	int<lower=0> N[Nobs];
	int<lower=0> Species[Nobs];
}

parameters {
	vector[Nsp] piSRaw;
	real<lower=0> sigmaG;
	real<lower=0, upper=1> piG;
}

transformed parameters {
	vector<lower=0, upper=1>[Nsp] piS;
	
  for(s in 1:Nsp){
  	piS[s] <- inv_logit(logit(piG) + sigmaG*piSRaw[s]);
  }
  
  }

model {
	piSRaw ~ normal(0,1);
	sigmaG ~ cauchy(0,2);
	
	for(i in 1:Nobs){
		Pos[i] ~ binomial(N[i], piS[Species[i]]);
	}
}

generated quantities {
	vector<lower=0>[Nobs] yNew;
	
	for(n in 1:Nobs){
		yNew[n] <- binomial_rng(N[n], piS[Species[n]]);
	}
}
"
sm <- stan_model(model_code = modelString)
shinystan(sm)
W.data

fit <- sampling(sm, data = W.data, chains = 4)
head(fit)

pis <- extract(fit, "piS")
str(pis)
plot(density(pis$piS))
hist(pis$piS)

fm <- as.mcmc(fit)
str(fit)

my_shiny <- as.shinystan(fit)

launch_shinystan(my_shiny)

