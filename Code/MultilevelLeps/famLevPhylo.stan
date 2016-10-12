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
  matrix[nFam, nFam] pCor;
  
  pCor = cholesky_decompose(phyCor);
}

parameters {
  // Vector of standard normal deviates for piS noncentering trick  
  vector[nFam] zFam;
  vector[nSpp] zSpp; 
  vector<lower=0>[nFam] sigmaF;
  real<lower=0> tau;

  real<lower=0> sigmaG; 

  real globalLogit;  
}

transformed parameters {
  vector[nFam] famLogit;
  vector[nSpp] sppLogit; 
 
  famLogit = globalLogit +(diag_pre_multiply(rep_vector(sigmaG, nFam), pCor) * zFam);

  
  for (s in 1:nSpp)
     sppLogit[s] = famLogit[Family[s]] + sigmaF[Family[s]] * zSpp[s];
    
}

model {
  globalLogit ~ normal(0,5);
  sigmaG ~ student_t(3,0,5);
  zFam ~ normal(0,1);
  zSpp ~ normal(0,1); 
  sigmaF ~ student_t(3,0,tau); //Folded cauchy dist. for global sd.
  tau ~ student_t(3,0,1);
  
  // Likelihood; Positive infection is binomially distributed w/ spp-level probabilities
  for (n in 1:nObs)
    Pos[n] ~ binomial_logit(N[n], sppLogit[Species[n]]);
}

generated quantities {
  /* New data for posterior predictive check simulated according to our model. If model does a good job then the new data should be very close to actual data. */

  vector[nObs] log_lik;
  int yNew[nObs];
  vector[nFam] infectF;
  vector[nSpp] infectS;
  real infectG;
  
     
  infectG = inv_logit(globalLogit);
 
  for(f in 1:nFam) 
    infectF[f] = inv_logit(famLogit[f]);
    
  
  for(s in 1:nSpp) 
    infectS[s] = inv_logit(sppLogit[s]);
    
  
    for(n in 1:nObs) {
      yNew[n] = binomial_rng(N[n], infectS[Species[n]]);
      log_lik[n] = binomial_logit_log(Pos[n], N[n], sppLogit[Species[n]]);
    }
}
