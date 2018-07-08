functions {
  
} 
data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix
  int<lower=1> P;
  // last time at-risk time
  vector<lower=0>[N] y;
  // indicator of event observation (0 = censored)
  vector<lower=0,upper=1>[N] dy;
  // covariate matrix (has no intercept)
  matrix[N, P] xmat;
  // prior means for log(kappa) and log(alpha) parameters
  real log_kappa_pmean;
  real log_alpha_pmean;
}

parameters {
  vector[P] beta_raw;
  // log_alpha_raw scaled weird to make U(-2, 2) init more reasonable
  real log_alpha_raw; 
  // aux = kappa^(-1/alpha) or exp(-log(kappa)/alpha)
  // so kappa = aux^(-alpha)
  real<lower=0> aux;
}

transformed parameters {
  vector[P] beta = beta_raw * 2.5; // for N(0, 2.5) prior
  real log_alpha = log_alpha_pmean + log_alpha_raw ; /// 5; // ncp
  real<lower=0> alpha = exp(log_alpha);
  real<lower=0> kappa = aux^(-alpha);
  real log_kappa = log(kappa); // N(log_kappa_pmean, 1) prior
}

model {
  vector[N] lp = xmat * beta + log_kappa;
  vector[N] scale = exp(-lp/alpha);
  
  // priors
  target += normal_lpdf(beta_raw | 0, 1);
  //target += normal_lpdf(log_alpha_raw | 0, 5); 
  //target += normal_lpdf(log_kappa | log_kappa_pmean, 1);
  
  // likelihood
  for (n in 1:N) {
    if (dy[n] == 1) {
      target += weibull_lpdf(y[n] | alpha, scale[n]);  
    } else {
      target += weibull_lccdf(y[n] | alpha, scale[n]);
    }
  }
}

