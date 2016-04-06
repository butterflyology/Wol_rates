# Beta-binomial model for overdispersion applied to Wolbachia infection rates

library("coda")
library("lattice")
library("rstan")
library("shinystan")
library("dplyr")
library("geiger")

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

table(WD.sp$Species)

# break down the infection freqs by species

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


by.sp <- WD.sp %>% group_by(Species) %>% summarise(N = sum(n), sum(Positive))
tail(by.sp)
dim(by.sp)

pro.inf <- by.sp[, 3] / by.sp[, 2]
# pdf(file = "Prop-spec.pdf", bg = "white")
hist(pro.inf[, 1], col = "grey", main = "Proportion infected by species", xlab = "Infection frequency", las = 1, ylim = c(0, 200), xlim = c(0, 1), breaks = 20)
# dev.off()

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

bn.conf.exact(100, 200, conf.level = 0.95)
bn.conf.exact(1, 2, conf.level = 0.95)

colSums(WD.sp[, 2:4])
bn.conf.exact(sum(WD.sp[, 3]), sum(WD.sp[, 2]), conf.level = 0.95)

bn.conf.exact(sum(WD.sp[WD.sp$Species == "Acraea_encedon", 3]), sum(WD.sp[WD.sp$Species == "Acraea_encedon", 2]), conf.level = 0.95)


#####
##### hierarchical Bayesian fun time
#####

Pos <- WD.sp$Positive
N <- WD.sp$n
Nobs <- length(WD.sp$Positive)
Nsp <- length(unique(WD.sp$index))
Species <- as.numeric(WD.sp$index) # Need an index to associate Nobs with Nsp, this does not work


W.data <- list(Pos = Pos, Nobs = Nobs, Nsp = Nsp, Species = Species, N = N)
W.data

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
# shinystan(sm)


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


#####
##### By family
#####

Fams <- read.csv("LepFams.csv", header = TRUE)
head(Fams)
str(Fams)

Fams <- Fams[complete.cases(Fams), ]
Fams$index <- as.factor(as.numeric(Fams$Family))
summary(Fams)
head(Fams)

dim(Fams[Fams$index == 16, ])
summary(Fams)

dim(Fams[Fams$Family == "Papilionidae", ])

(Fams[length(Fams$Family == 1), ])
table(Fams$Family)


#####
##### 27 Jan 2016 code
#####

fams <- read.csv("~/Desktop/Projects/Wolbachia_rates/LepFams.csv")

LepTree <- read.nexus("./Regier_data/journal.pone.0058568.s009.TRE")
plot(LepTree, cex = 0.2)

lam <- 10^(-1:4)
cv <- sapply(lam, function(x) sum(attr(chronopl(LepTree, lambda = x, CV = TRUE), "D2")))

plot(lam, cv, pch = 19, ylab = "cross-validation score", xlab = expression(paste(lambda)), las = 1, cex = 1.5) # lowest CV. Small lambda means that every branch gets its own rate, larger and you are more clock like. 