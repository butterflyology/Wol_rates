
######################################################################################################################################
######################################################################################################################################
### Simulating trait evolution on a tree.
######################################################################################################################################
######################################################################################################################################

BrownianEvolveParameters <- function(phy, start.value, rate, logspace = TRUE){
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1
    evolved.parameter <- integer(ntips + phy$Nnode)
    evolved.parameter[ROOT] <- start.value
    anc <- phy$edge[, 1]
    des <- phy$edge[, 2]
    edge.length <- phy$edge.length
    for(i in N:1) {
        evolved.parameter[des[i]] <- rnorm(1, mean = evolved.parameter[anc[i]], sd = sqrt(edge.length[i] * rate))
    }
    if(logspace == TRUE){
        return(exp(evolved.parameter))
    }else{
        return(evolved.parameter)
    }
}

t1 <- pbtree(b = 0.5, d = 0.001, n = 100, nsim = 1, type = "cont")
plot(t1)

bt1 <- BrownianEvolveParameters(phy = t1, start.value = 0, rate = 0.0001)
bt1


OUEvolveParameters <- function(phy, alpha, sigma.sq, mean, logspace=TRUE){
    phy <- reorder(phy, "postorder")
    ntips <- length(phy$tip.label)
    N <- dim(phy$edge)[1]
    ROOT <- ntips + 1
    evolved.parameter <- integer(ntips + phy$Nnode)
    evolved.parameter[ROOT] <- mean
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    edge.length <- phy$edge.length
    for(i in N:1) {
        evolved.parameter[des[i]] = evolved.parameter[anc[i]] * exp(-alpha*edge.length[i])+(mean)*(1-exp(-alpha*edge.length[i]))+sigma.sq*rnorm(1,0,1)*sqrt((1-exp(-2*alpha*edge.length[i]))/(2*alpha))
    }
    if(logspace==TRUE){
        return(exp(evolved.parameter))
    }else{
        return(evolved.parameter)
    }
}
