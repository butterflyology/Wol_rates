data {
  int<lower=0> nObs; // number of observations
  int<lower=0> nSpp; // number of species
  int<lower=0> nFam;
  int<lower=0> Pos[nObs]; // data (number of positive infections)
  int<lower=0> N[nObs]; // Total number of trials per sample
  int<lower=0> Species[nObs]; // species as integers (1-301)
  int<lower=0> Family[nSpp];
  corr_matrix[nFam] phyCor;
}

transformed data {
  matrix[nFam, nFam] d;
  
  d <- diag_matrix(rep_vector(1, nFam));

//   matrix[nFam, nFam] pCor;
//   
//   pCor <- cholesky_decompose(phyCor);
}

parameters {
  // Vector of standard normal deviates for piS noncentering trick  
  real<lower=0, upper=1> lambda;
  vector[nFam] zFam;
  vector[nSpp] zSpp; 
  vector<lower=0>[nFam] sigmaF;
  // global std dev hyperprior	 for normal                    
  real<lower=0> sigmaG; 
  #vector<lower=0>[nFam] sigmaF;
  #real<lower=0> sigmaF;
  // global prob of infection across spp; implies beta(1,1) == U(0,1) prior  
  real globalLogit;  
}

transformed parameters {
  matrix[nFam, nFam] pCor;
  vector[nFam] famLogit;
  vector[nSpp] sppLogit; 
  
  {
  matrix[nFam, nFam] lamCor;
  
  lamCor <- lambda * (phyCor - d) + d;
  pCor <- cholesky_decompose(lamCor);
  }
  
  famLogit <- globalLogit +(diag_pre_multiply(rep_vector(sigmaG, nFam), pCor) * zFam);
  
  
  for(s in 1:nSpp) 
    sppLogit[s] <- famLogit[Family[s]] + 
      sigmaF[Family[s]] * zSpp[s];
}

model {
  lambda ~ beta(2,2);
  globalLogit ~ normal(0,1);
  sigmaG ~ cauchy(0,5);
  zFam ~ normal(0,1);
  zSpp ~ normal(0,1); 
  sigmaF ~ cauchy(0,5); //Folded cauchy dist. for global sd.
  
  
  // Likelihood; Positive infection is binomially distributed w/ spp-level probabilities
  
  for(i in 1:nObs) 
    Pos[i] ~ binomial_logit(N[i], sppLogit[Species[i]]);
}

generated quantities {
  /* New data for posterior predictive check simulated according to our model. If model does a good job then the new data should be very close to actual data. */
    
    vector[nObs] log_lik;
  vector<lower=0>[nObs] yNew;
  vector[nFam] infectF;
  vector[nSpp] infectS;
  real infectG;
  
  infectG <- inv_logit(globalLogit);
  
  
  for(f in 1:nFam)
    infectF[f] <- inv_logit(famLogit[f]);
  
  
  for(s in 1:nSpp)
    infectS[s] <- inv_logit(sppLogit[s]);
  
  
  for(n in 1:nObs) {
    yNew[n] <- binomial_rng(N[n], infectS[Species[n]]);
    log_lik[n] <- binomial_logit_log(Pos[n], N[n], sppLogit[Species[n]]);
  }
}
