data {
  int<lower=0> Nobs; // number of observations
  int<lower=0> Nsp; // number of species
  int<lower=0> Pos[Nobs]; // data (number of positive infections)
  int<lower=0> N[Nobs]; // Total number of trials per sample
  int<lower=0> Species[Nobs]; // species as integers (1-301)
}

parameters {
// Vector of standard normal deviates for piS noncentering trick  
  vector[Nsp] piSRaw; 

// global std dev hyperprior	 for normal                    
  real<lower=0> sigmaG; 

// global prob of infection across spp; implies beta(1,1) == U(0,1) prior  
  real<lower=0, upper=1> piG;  
}

transformed parameters {
  vector<lower=0, upper=1>[Nsp] piS; 

/* spp-lev prob of infection. piS is dist. as the inverse logit of a normal(logit(piG),sigmaG). The math wankery is a noncentering trick (sometimes called the "Matt trick") that reduces sampling autocorrelation between mean and variance that often occurs. It's easier for the sampler to explore parameter space of a N(0,1). */ 
  
  for(s in 1:Nsp){
  	piS[s] <- inv_logit(logit(piG) + sigmaG * piSRaw[s]);
  }
  
}

model {
  piSRaw ~ normal(0,1); 
  sigmaG ~ cauchy(0,2); //Folded cauchy dist. for global sd.
	
// Likelihood; Positive infection is binomially distributed w/ spp-level probabilities
  
  for(i in 1:Nobs){
    Pos[i] ~ binomial(N[i], piS[Species[i]]);
  }
}

generated quantities {
/* New data for posterior predictive check simulated according to our model. If model does a good job then the new data should be very close to actual data. */ 
  
  vector<lower=0>[Nobs] yNew;
  vector<lower=0,upper=1>[Nsp] infectNew;	
  real<lower=0, upper=1> infectG;
  
  infectG <- bernoulli_rng(piG);
  
  for(s in 1:Nsp){
  	infectNew[s] <- bernoulli_rng(piS[s]);
  }
  
  for(n in 1:Nobs){
	yNew[n] <- binomial_rng(N[n], piS[Species[n]]);
  }
}