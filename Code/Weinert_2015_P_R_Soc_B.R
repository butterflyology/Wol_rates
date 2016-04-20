######################################################################################################################
######################################################################################################################
########## These R scripts were used to produce results in Weinert, Araujo-Jnr, Ahmed and Welch              #########
########## "The incidence of bacterial endosymbionts in terrestrial arthropods" (Proc. R. Soc. B.)           ######### 
##########  Eq. numbers refer to the Supplementary Methods  (c) 2014.                                        ######### 
###################################################################################################################### 
###################################################################################################################### 

###################################################################################################################### 
##### Calculate a beta binomial log likelihood term (see eq. 11)
lnBetabinomial <- function(mu, rho, n, k) {
  a <- get.a(mu, rho)
  b <- get.b(mu, rho)
  if( any(!is.finite(unique(c(a, b)))) || any(unique(c(a, b)) <= 0) ) {
    return(-1e10)    
  }
  x <- lchoose(n, k) + lbeta(a + k, b + n - k) - lbeta(a, b)
  return(x)
}
##### Beta binomial term with a pooled sample  (see eq. 3)
lnPooledBetabinomial <- function(mu, rho, n, k, m)
{
  a <- get.a(mu, rho)
  b <- get.b(mu, rho)
  if( any(!is.finite(unique(c(a, b)))) || any(unique(c(a, b)) <= 0) ) { 
    return(-1e10)    
  }
  ans <- numerical_integration(a, b, n, k, m)
  return(ans) 
}
###################################################################################################################### 
### Numerical integration for estimating incidence with pooled samples;  (see eqs. 3, 12 and 13)
numerical_integration <- function(a, b, n, k, m) 
{
	A <- cbind(n, k, m)
  	a2 <- apply(A, MARGIN = 1, FUN = function(d, a, b, mfunc, bounds) {
  		integrate(doubexpfunc, lower = bounds[1], upper = bounds[2], a = a, b = b, n = d[1], k = d[2], m = d[3], stop.on.error = FALSE) 
  		}, a = a, b = b, mfunc = mfunc, bounds = c(-Inf, Inf))
	ans <- unlist(lapply(a2, FUN = function(z) { 
		z$value
		}))
	return(log(ans))
}
### The log of the hyperbolic cosine
logcosh <- function(t){
  a <- vector('numeric', length(t))
  w <- which(t < 0)
  if(length(w) > 0)
    a[w] <- log1p(exp(2 * t[w])) - log(2) - t[w]
  w <- which(t >= 0)
  if(length(w) > 0)
    a[w] <- log1p(exp(-2 * t[w])) - log(2) + t[w]
  w <- which(abs(t) < 50)
  if(length(w) > 0)
    a[w] <- log(cosh(t[w])) 
  a
}
### Double expontential transformation of the integral
doubexpfunc <- function(x, a, b, n, k, m){
  # Useful
  p2 <- pi / 2
  p2sinhx <- p2 * sinh(x)
  l2 <- log(2)
  logcoshp2sinhx <- logcosh(p2sinhx)
  
  # precalculate log(x) and log(1-x)
  logx   <- plogis(2 * p2sinhx, log.p = TRUE) # sigmoidal function
  log1mx <- -p2sinhx - logcoshp2sinhx - l2
  if(any(na.omit(log1mx) > 0))
    log1mx[which(log1mx > 0)] <- 0
  # Calculate log(f(phi(x)))
  ans <- k*log1p(-exp(m * log1mx)) + log1mx * (m * (n - k) + b - 1) + logx * (a - 1) - lbeta(a, b)
  # Calculate f(phi(x)) * dphi(x)/dx
  ans <- exp(ans + log(p2) + logcosh(x) - 2 * logcoshp2sinhx - l2 + lchoose(n, k))
  
  # Check for non-finite values
  if(any(!is.finite(ans)))
    ans[which(!is.finite(ans))] <- 0
  
  if(any(x == 0) && k > 0)
    ans[which(x == 0)] <- 0
  ans
}
### The untransformed integral
intfunc <- function(x, a, b, n, k, m){
	logx   <- log(x)
	log1mx <- log(1 - x)
	ans <- exp(k * log1p(-exp(m * log1mx)) + log1mx * (m * (n - k) + b - 1) + logx * (a - 1) - lbeta(a, b) + lchoose(n, k))
	ans
}
###################################################################################################################### 
### Summation method for estimating incidence with pooled samples. (see eq. 4)
sfunc <- function(Z, a, b)
{
	n <- Z[1]; k <- Z[2]; m <- Z[3]
	i <- 0:k
	s <- log(sum((-1)^i * exp(lchoose(k, i) + lbeta(a, b + m * (n - k + i)) - lbeta(a, b))))
	s <- s + lchoose(n, k)
}
summation_method <- function(a, b, n, k, m){
	A <- cbind(n, k, m)
	ans <- apply(A, MARGIN = 1,FUN = sfunc, a = a, b = b)
	if(any(!is.finite(ans)))
    	ans[which(!is.finite(ans))] <- -1e5
  	return(ans)
}

###################################################################################################################### 
### Useful functions to transform between parameterisations (see, e.g., eqs. 8-10, 15-16) 
##### Calculate a (alpha) from mu and rho
get.a <- function(mu, rho) {
  a <- ((1 / rho) - 1) * mu
  return(a)
}
##### Calculate b (beta) from mu and rho
get.b <- function(mu, rho){
  b <- ((1 / rho) - 1) * (1 - mu)
  return(b)
}
##### Calculate integral over a Beta dist from c to 1
get.x <- function(mu, rho, c){
    return(1 - pbeta(c, get.a(mu, rho), get.b(mu, rho))) 
}
### Calculate F, mu and rho
get.F <- function(mu, rho, pp, ga) {
  h <- pp * ga + mu * (1 - pp)
  return(1 + ((rho - 1) * (1 - pp) * mu * (1 - mu)) / (h * (1 - h)))
}
get.mu <- function(a, b) {
  return(a / (a + b))
}
get.rho <- function(a, b) {
  return(1 / (1 + a + b))
}

###################################################################################################################### 
##### The likelihood function for the 2-parameter Beta distribution model (eqs. 2,5); without summation           ####
###################################################################################################################### 
lnL.Beta <- function(mu, rho, n, k, m){
  # If we have no infected pooled samples return the standard LnL
  inf.pools <- which(m > 1 & k > 0)
  if(length(inf.pools) == 0)
    return(lnBetabinomial(mu, rho, n * m, k))
  
  lnL.terms <- vector(mode = 'numeric', length = length(n))
  lnL.terms[-inf.pools] <- lnBetabinomial(mu, rho, n[-inf.pools] * m[-inf.pools], k[-inf.pools])
  lnL.terms[inf.pools]  <- lnPooledBetabinomial(mu, rho, n[inf.pools], k[inf.pools], m[inf.pools])
  return(lnL.terms)
}
###################################################################################################################### 
##### The likelihood function for the 4-parameter doubly-inflated beta distribution model (eqs. 5, 14)            ####
###################################################################################################################### 
lnL.Infl <- function(par,n,k,m,a,b,gamma,phi){
  
  mu  <- get.mu(a, b)
  rho <- get.rho(a, b)
 
  w1 <- which(k == n)
  w0 <- which(k == 0)
  
  excl <- c(w0, w1)
  le <- length(excl)
  if(le > 0) {	
    L <- 0
    if(le < length(n))
    {
      if(phi == 1)
        return(-1e5)
      
      L <- L + length(n[-excl]) * log(1 - phi) + sum(lnL.Beta(mu, rho, n[-excl], k[-excl], m[-excl]))
    }
    if(any(w0))
    {
      L <- L + sum(log(phi * (1 - gamma) + (1 - phi) * exp(lnL.Beta(mu, rho, n[w0], k[w0], m[w0]))))
    }	
    if(any(w1))
    {
      L <- L + sum(log(phi * gamma + (1 - phi) * exp(lnL.Beta(mu, rho, n[w1], k[w1], m[w1]))))
    }
    if(any(!is.finite(L))){ 
    	L[which(!is.finite(L))] <- -1e5 
    	}
    return(L)
  } else  {
  s <- sum(lnL.Beta(mu, rho, n, k, m))
  if(any(!is.finite(s))){ 
  	s[which(!is.finite(s))] <- -1e5 
  	}
  return(s)
  }
}