data {
  int<lower=0> nObs; // number of observations
  int<lower=0> nSpp; // number of species
  int<lower=0> nFam;
  int<lower=0> Pos[nObs]; // data (number of positive infections)
  int<lower=0> N[nObs]; // Total number of trials per sample
  int<lower=0> Species[nObs]; // species as integers (1-301)
  int<lower=0> Family[nSpp];
}

parameters {
  // Vector of standard normal deviates for piS noncentering trick  
  vector[nFam] thetaFRaw;
  vector[nSpp] thetaSRaw; 
  
  // global std dev hyperprior	 for normal                    
  real<lower=0> sigmaG; 
  #vector<lower=0>[nFam] sigmaF;
  real<lower=0> sigmaF;
  // global prob of infection across spp; implies beta(1,1) == U(0,1) prior  
  real<lower=0, upper=1> thetaG;  
}

transformed parameters {
  vector<lower=0, upper=1>[nFam] thetaF;
  vector<lower=0, upper=1>[nSpp] thetaS; 
  
 for(f in 1:nFam) 
   thetaF[f] <- inv_logit(logit(thetaG) + sigmaG * thetaFRaw[f]);
 
 
  for(s in 1:nSpp) {
    thetaS[s] <- inv_logit(logit(thetaF[Family[s]]) + 
      sigmaF * thetaSRaw[s]);
  }    
  
}

model {
  sigmaG ~ cauchy(0,5);
  thetaFRaw ~ normal(0,1);
  thetaSRaw ~ normal(0,1); 
  sigmaF ~ cauchy(0,5); //Folded cauchy dist. for global sd.

  
// Likelihood; Positive infection is binomially distributed w/ spp-level probabilities

  for(i in 1:nObs) {
    Pos[i] ~ binomial(N[i], thetaS[Species[i]]);
  }
}

generated quantities {
/* New data for posterior predictive check simulated according to our model. If model does a good job then the new data should be very close to actual data. */ 

  vector[nObs] log_lik;
  vector<lower=0>[nObs] yNew;
  vector<lower=0,upper=1>[nSpp] infectNewS;	
  vector<lower=0, upper=1>[nFam] infectNewF;
  real<lower=0, upper=1> infectG;
  
  infectG <- bernoulli_rng(thetaG);

  for(f in 1:nFam) 
    infectNewF[f] <- bernoulli_rng(thetaF[f]);
  

  for(s in 1:nSpp)
    infectNewS[s] <- bernoulli_rng(thetaS[s]);


  for(n in 1:nObs) {
    yNew[n] <- binomial_rng(N[n], thetaS[Species[n]]);
    log_lik[n] <- binomial_log(Pos[n], N[n], thetaS[Species[n]]);
  }
}