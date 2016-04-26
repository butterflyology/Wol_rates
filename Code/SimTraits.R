### Simulating trait evolution on a tree. By J. M. Beulieu and C.A. Hamm
set.seed(87624781)


library("geiger")
library("phytools")
library("boot")

invLogit <- function(x){ 
	out <- exp(x) / (1 + exp(x)
	return(out)
	}


BrownianEvolveParameters <- function(phy, start.value, rate, inverselogit = TRUE){
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1
    evolved.parameter <- integer(ntips)
    evolved.parameter[ROOT] <- start.value
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    edge.length <- phy$edge.length
    for(i in N:1) {
        evolved.parameter[des[i]] <- rnorm(1, mean = evolved.parameter[anc[i]], sd = sqrt(edge.length[i] * rate))
    }
    if(inverselogit == TRUE){
        return(inv.logit(evolved.parameter[1:ntips]))
    }else{
        return(evolved.parameter[1:ntips])
    }
}


t1 <- pbtree(b = 0.5, d = 0.0001, n = 100, nsim = 1, type = "cont")
plot(t1)

bt1 <- BrownianEvolveParameters(phy = t1, start.value = 0, rate = 0.0001, inverselogit = TRUE)
bt1
hist(bt1)

OUEvolveParameters <- function(phy, alpha, sigma.sq, mean, inverselogit = TRUE){
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1
    evolved.parameter <- integer(ntips + phy$Nnode)
    evolved.parameter[ROOT] <- mean
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    edge.length <- phy$edge.length
    for(i in N:1) {
        evolved.parameter[des[i]] = evolved.parameter[anc[i]] * exp(-alpha * edge.length[i]) + (mean) * (1 - exp(-alpha * edge.length[i])) + sigma.sq * rnorm(1, 0, 1) * sqrt((1 -exp(-2 * alpha * edge.length[i])) / (2 * alpha))
    }
    if(inverselogit == TRUE){
        return(inv.logit(evolved.parameter[1:ntips]))
    }else{
        return(evolved.parameter[1:ntips])
    }
}


ou1 <- OUEvolveParameters(phy = t1, alpha = 0.9, sigma.sq = 0.005, mean = 0.9, inverselogit = TRUE)
ou1
hist(ou1)


