functions {
  
} 
data {
  // number of observations
  int<lower=0> N;
  // number of columns in design matrix
  int<lower=1> P;
  // last time at-risk
  vector<lower=0>[N] y;
  // indicator of event observation
  vector<lower=0,upper=1>[N] dy;
  // covariate matrix 
  matrix[N, P] xmat;
  // prior means for log(kappa) and log(alpha) parameters
  real log_kappa_pmean;
  real log_alpha_pmean;
}

parameters {
  vector[P] beta_raw;
  real log_alpha_raw;
  real log_kappa_raw;
}

transformed parameters {
  vector[P] beta = beta_raw * 2.5;
  real log_kappa = log_kappa_pmean + log_kappa_raw; //* (log(100) / 2);
  real log_alpha = log_alpha_pmean + log_alpha_raw / 5;
  real<lower=0> kappa = exp(log_kappa);
  real<lower=0> alpha = exp(log_alpha);
}

model {
  vector[N] lp = xmat * beta + log_kappa;
  vector[N] scale = exp(-lp/alpha);

  // priors
  target += normal_lpdf(beta_raw | 0, 1);
  target += normal_lpdf(log_alpha_raw | 0, 5);
  target += normal_lpdf(log_kappa_raw | 0, 1);
  
  // likelihood
  for (n in 1:N) {
    if (dy[n] == 1) {
      target += weibull_lpdf(y[n] | alpha, scale[n]);  
    } else {
      target += weibull_lccdf(y[n] | alpha, scale[n]);
    }
  }
}

