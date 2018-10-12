functions {
 
} 

data {
  int<lower=0> N;
  vector<lower=0, upper=1>[N] z;
  int<lower=1> P;
  matrix[N, P] x;
  vector<lower=0>[N] yr;
  vector<lower=0>[N] yt;
  int<lower=0,upper=1> dyr[N];
  int<lower=0,upper=1> dyt[N];
  int<lower=0,upper=1> use_priors;
  vector[6] log_alpha_pmean;
  vector[6] log_kappa_pmean;
  real<lower=0> sigma_pa;
  real<lower=0> sigma_pb;
}

transformed data {
  vector[N] type;
  int start_i[N];
  real log_kappa_psd = log(100) / 2;
  real log_alpha_psd = 2;
  
  for (n in 1:N) {
    start_i[n] = (z[n] == 1) ? 4 : 1;
    if (dyr[n] == 0 && dyt[n] == 0) {
      type[n] = 1; // type 1: observe neither event
    } else if (dyr[n] == 1 && dyt[n] == 0) {
      type[n] = 2; // type 2: observe non-terminal but terminal censored
    } else if (dyr[n] == 0 && dyt[n] == 1) {
      type[n] = 3; // type 3: observed terminal with no prior non-terminal
    } else if (dyr[n] == 1 && dyt[n] == 1) {
      type[n] = 4; // type 4: both non-terminal and terminal observed
    }
  } // end looping over n
}

parameters {
  // vectors of regression parameters
  matrix[P, 3] beta;
  
  // shape parameters (the one in exponent of time)
  // alpha > 1 -> hazard increases over time, more clumping
  vector[6] log_alpha_raw; 
  
  // scale parameters (without shift to make expected event time near 1)
  // bigger kappa -> faster event occurrence
  vector[6] log_kappa_raw;
  
  // variance of frailties
  real<lower=0> sigma;
}

transformed parameters {
  vector[6] log_kappa = log_kappa_pmean + log_kappa_raw * log_kappa_psd;
  vector[6] log_alpha = log_alpha_pmean + log_alpha_raw * log_alpha_psd;
  vector[6] kappa = exp(log_kappa);
  vector[6] alpha = exp(log_alpha);
  real inv_sigma = 1 / sigma;
}

model {
  matrix[N, 3] lp = x * beta;
  vector[N] ll;
  int i;
  real lp1;
  real lp2;
  real lp3;
  
  if (use_priors == 1) {
    to_vector(beta) ~ normal(0, 2.5);
    log_alpha_raw ~ normal(0, 1);
    log_kappa_raw ~ normal(0, 1);
    // TODO(LCOMM): convince myself once and for all that following line
    // is correct...
    // prior mean for sigma is sigma_pb / (sigma_pa - 1) for sigma_pa > 1
    // prior mode for sigma is sigma_pb / (sigma_pa + 1) for sigma_pa > 1
    // for large and equal (sigma_pa, sigma_pb), variance is roughly 1 / sigma_pa
    sigma ~ inv_gamma(sigma_pa, sigma_pb);
  }
    
  // likelihood
  for (n in 1:N) {
    i = start_i[n];
    lp1 = lp[n, 1] + log_kappa[i];
    lp2 = lp[n, 2] + log_kappa[i + 1];
    lp3 = lp[n, 3] + log_kappa[i + 2];
    
    if (type[n] == 1) {
      target += -inv_sigma * log1p(sigma * 
                (-weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
                 -weibull_lccdf(yt[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1]))));
    } else if (type[n] == 2) {
      target += (weibull_lpdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
                -weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i]))) +
                -(inv_sigma + 1) * log1p(sigma * 
                (-weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
                 -weibull_lccdf(yr[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1])) + 
                 -weibull_lccdf(yt[n] - yr[n] | alpha[i + 2], exp(-(lp3)/alpha[i + 2]))));
    } else if (type[n] == 3) {
      target += weibull_lpdf(yt[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1])) +
                -weibull_lccdf(yt[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1])) +
                -(inv_sigma + 1) * log1p(sigma * 
                (-weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
                 -weibull_lccdf(yt[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1]))));
    } else if (type[n] == 4) {
      target += log1p(sigma) + 
                (weibull_lpdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
                -weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i]))) +
                (weibull_lpdf(yt[n] - yr[n] | alpha[i + 2], exp(-(lp3)/alpha[i + 2])) +
                -weibull_lccdf(yt[n] - yr[n] | alpha[i + 2], exp(-(lp3)/alpha[i + 2]))) +
                -(inv_sigma + 2) * log1p(sigma * 
                (-weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
                 -weibull_lccdf(yr[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1])) + 
                 -weibull_lccdf(yt[n] - yr[n] | alpha[i + 2], exp(-(lp3)/alpha[i + 2]))));
    }
  }
}

