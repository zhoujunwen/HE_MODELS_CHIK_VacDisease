data {
     int<lower=0> Nobs;
     int Npos[Nobs];
     int Ntotal[Nobs];
     int Age[Nobs];
     int <lower=1>Ymax;
     matrix[Nobs, Ymax] AgeExpoMatrix;
}

parameters {
  row_vector<lower=0>[Ymax] foi;
  real<lower=0> sigma;
  real<lower=0, upper=1> p;   // Mixture weight for the bimodal prior
  real mu1;                   // Mean of the first normal distribution
  real<lower=0> sigma1;       // Standard deviation of the first normal distribution
  real mu2;                   // Mean of the second normal distribution
  real<lower=0> sigma2;       // Standard deviation of the second normal distribution
}
// original model 
// parameters {
//    row_vector<lower=0>[Ymax] foi;
//    real<lower=0> sigma;
// }


transformed parameters {
  real P[Nobs];
  real ScalerDotProduct[Nobs];
  
  for (i in 1:Nobs){
   ScalerDotProduct[i] = dot_product(AgeExpoMatrix[i,], foi);
   P[i] = 1 - exp(-ScalerDotProduct[i]);
 }
}

model {
  for (i in 1:Nobs)
    Npos[i] ~ binomial(Ntotal[i], P[i]);
  // Bimodal prior for foi using target +=
  for (ind in 1:Ymax)
    target += log_mix(p,
                      normal_lpdf(foi[ind] | mu1, sigma1),
                      normal_lpdf(foi[ind] | mu2, sigma2));
}
// original model 
// model {
//   for (i in 1:Nobs)
//     Npos[i] ~ binomial(Ntotal[i], P[i]) ;
//     sigma ~ cauchy(0, 1);
  
//   for(i in 2:Ymax)
//     foi[i] ~ normal(foi[i - 1], sigma);
//     foi[1] ~ normal(0, 1);
// }


generated quantities{
  vector[Nobs] Npos_sim;
  vector[Nobs] P_sim;
  vector[Nobs] logLikelihood;
  for(i in 1:Nobs){
    Npos_sim[i] = binomial_rng(Ntotal[i], P[i]);
    P_sim[i] = Npos_sim[i] / Ntotal[i];
    logLikelihood[i] = binomial_lpmf(Npos[i] | Ntotal[i], P[i]);
  }
}
