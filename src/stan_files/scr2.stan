functions {
  vector ll_scr2(vector yr, vector yt, int[] dyr, int[] dyt, vector tx, 
                  vector alpha, vector kappa, 
                  matrix x_1, matrix x_2, matrix x_3, 
                  vector omega1, vector omega2, vector omega3, 
                  vector beta, vector gamma_obs) {
  int N = num_elements(yr);
  vector[N] ll;
  vector[N] lp1_base;
  vector[N] lp2_base;
  vector[N] lp3_base;
  real lp1;
  real lp2;
  real lp3;
  int i;
  
  // linear predictors for gamma of 0
  lp1_base = x_1 * omega1;
  lp2_base = x_2 * omega2;
  lp3_base = x_3 * omega3;
  
  for (n in 1:N) {
    i = (tx[n] == 1) ? 4 : 1; // model index start
    lp1 = lp1_base[n] + beta[i] * gamma_obs[n] + kappa[i];
    lp2 = lp2_base[n] + beta[i + 1] * gamma_obs[n] + kappa[i + 1];
    lp3 = lp3_base[n] + beta[i + 2] * gamma_obs[n] + kappa[i + 2];
  
    if (dyr[n] == 0 && dyt[n] == 0) {
      // type 1: observe neither event
      ll[n] = weibull_lccdf(yr[n] | alpha[i], exp(-(lp1)/alpha[1])) + 
              weibull_lccdf(yt[n] | alpha[i + 1], exp(-(lp2)/alpha[2]));
    } else if (dyr[n] == 1 && dyt[n] == 0) {
      // type 2: observe non-terminal but terminal censored
      ll[n] = weibull_lpdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
              weibull_lccdf(yr[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1])) + 
              weibull_lccdf(yt[n] - yr[n] | alpha[i + 2], exp(-(lp3)/alpha[i + 2]));
    } else if (dyr[n] == 0 && dyt[n] == 1) {
      // type 3: observed terminal with no prior non-terminal
      ll[n] = weibull_lccdf(yr[n] | alpha[1], exp(-(lp1)/alpha[i])) +
              weibull_lpdf(yt[n] | alpha[2], exp(-(lp2)/alpha[i + 1]));
    } else if (dyr[n] == 1 && dyt[n] == 1) {
      // type 4: both non-terminal and terminal observed
      ll[n] = weibull_lpdf(yr[n] | alpha[i], exp(-(lp1)/alpha[i])) +
              weibull_lccdf(yr[n] | alpha[i + 1], exp(-(lp2)/alpha[i + 1])) + 
              weibull_lpdf(yt[n] - yr[n] | alpha[i + 2], exp(-(lp3)/alpha[i + 2]));
    }
  }
  
  return ll;
  }  
  
} 
data {
  // number of observations
  int<lower=0> N;
  // binary treatment indicator
  vector<lower=0, upper=1>[N] tx;
  // number of columns in 1st design matrix, including intercept
  int<lower=1> P_1;
  // design matrix for non-terminal model
  matrix[N, P_1] x_1;
  // number of columns in 2nd design matrix (no intercept)
  int<lower=1> P_2;
  // design matrix for terminal model w/o non-terminal history
  matrix[N, P_2] x_2;
  // number of columns in 3rd design matrix (no intercept)
  int<lower=1> P_3;
  // design matrix for terminal model with non-terminal history (no intercept)
  matrix[N, P_3] x_3;
  // observed non-terminal time
  vector<lower=0>[N] yr;
  // indicator of event observation for non-terminal event
  int<lower=0,upper=1> dyr[N];
  // observed terminal time
  vector<lower=0>[N] yt;
  // indicator of event observation for terminal event
  int<lower=0,upper=1> dyt[N];
}

transformed data {
  real rho = 0.5; // covariance for arm-specific frailties
}

parameters {
  // vectors of regression parameters
  vector[P_1] omega1;
  vector[P_2] omega2;
  vector[P_3] omega3;

  // shape parameters (the one in exponent of time)
  // alpha > 1 -> hazard increases over time, more clumping
  vector<lower=0>[6] alpha;
  
  // scale parameters
  // bigger kappa -> faster event occurrence
  vector[6] kappa;
  
  // frailty coefficients
  vector[6] beta;
  
  // arm-specific frailties (for the observed arm only)
  vector[N] gamma_obs;
}

model {
  vector[N] ll;
  
  // TODO(LCOMM): add other priors
  target += normal_lpdf(gamma_obs | 0, 1);
  
  // likelihood
  ll = ll_scr2(yr, yt, dyr, dyt, tx, alpha, kappa, x_1, x_2, x_3,
               omega1, omega2, omega3, beta, gamma_obs);
               
  target += sum(ll);
}

