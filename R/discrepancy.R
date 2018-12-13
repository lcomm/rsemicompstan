#' Extract single iteration implied DGP from a Stan fit object
#' 
#' Used for simulating replicate data for discrepancy measures
#' 
#' @param b (Post-warmup) MCMC iteration number
#' @param stan_fit Stan fit object
#' @return List of data generation parameters structured like result from 
#' \code{\link{return_dgp_parameters}}
#' @export
extract_params_b <- function(b, stan_fit) {
  betas  <- as.array(extract(stan_fit, pars = "beta"))[[1]][b, , ]
  beta1  <- betas[ , 1]
  beta2  <- betas[ , 2]
  beta3  <- betas[ , 3]
  alphas <- as.array(extract(stan_fit, pars = "alpha"))[[1]][b, ]
  kappas <- as.array(extract(stan_fit, pars = "kappa"))[[1]][b, ]
  sigma  <- as.array(extract(stan_fit, pars = "sigma"))[[1]][b]
  dgps   <- list(treated = list(beta1 = beta1,
                                beta2 = beta2,
                                beta3 = beta3,
                                alpha1 = alphas[4],
                                alpha2 = alphas[5],
                                alpha3 = alphas[6],
                                kappa1 = kappas[4],
                                kappa2 = kappas[5],
                                kappa3 = kappas[6]),
                 control = list(beta1 = beta1,
                                beta2 = beta2,
                                beta3 = beta3,
                                alpha1 = alphas[1],
                                alpha2 = alphas[2],
                                alpha3 = alphas[3],
                                kappa1 = kappas[1],
                                kappa2 = kappas[2],
                                kappa3 = kappas[3]),
                 sigma = sigma)
  return(dgps)
}



#' Extract DGP implied by summary measures of a Stan fit object
#' 
#' Not totally sure this is useful yet
#' 
#' @param stan_fit Stan model fit object
#' @param summary_measure Summary measure to extract. Defaults to "mean"
#' @return Parameter list like that from \code{\link{return_dgp_parameters}}
#' @export
extract_params_fitsummary <- function(stan_fit, summary_measure = "mean") {
  betas  <- matrix(summary(stan_fit, pars = "beta")$summary[ , summary_measure],
                   ncol = 3, byrow = TRUE)
  beta1 <- betas[ , 1]
  beta2 <- betas[ , 2]
  beta3 <- betas[ , 3]
  alphas <- summary(stan_fit, pars = "alpha")$summary[ , summary_measure]
  kappas <- summary(stan_fit, pars = "kappa")$summary[ , summary_measure]
  sigma  <- summary(stan_fit, pars = "sigma")$summary[ , summary_measure]
  dgps   <- list(treated = list(beta1 = beta1,
                                beta2 = beta2,
                                beta3 = beta3,
                                alpha1 = alphas[4],
                                alpha2 = alphas[5],
                                alpha3 = alphas[6],
                                kappa1 = kappas[4],
                                kappa2 = kappas[5],
                                kappa3 = kappas[6]),
                 control = list(beta1 = beta1,
                                beta2 = beta2,
                                beta3 = beta3,
                                alpha1 = alphas[1],
                                alpha2 = alphas[2],
                                alpha3 = alphas[3],
                                kappa1 = kappas[1],
                                kappa2 = kappas[2],
                                kappa3 = kappas[3]),
                 sigma = sigma)
  return(dgps)
}



#' Function to make a replicate data set based on the bth iteration of a stan 
#' fit object
#' 
#' Used for discrepancy measure calculations. Matches the Z, design matrix, and 
#' censoring times of the original data
#' 
#' @param b Iteration number to extract parameter values
#' @param seed Seed for simulating the new data values
#' @param stan_fit Stan model fit object
#' @param z Length-n treatment assignment vector
#' @param x n x P design matrix (no intercept)
#' @param cens_times Potential censoring times to match
#' @return Data set simulated under parameters from bth posterior draw
#' @export
make_repdata <- function(b = 1, seed = 123, stan_fit, z, x, cens_times, ...) {
  n <- length(z)
  params <- extract_params_b(b = b, stan_fit = stan_fit)
  rep_dat <- simulate_from_param(n = n, seed = seed + b * n, 
                                 params = params, 
                                 data_match = list(z = z, x = x, cens_times = cens_times),
                                 censor = TRUE, ...)
  return(rep_dat)
}



#' Calculate KM survivals in data separately by treatment arm
#' 
#' Useful for calculation of the discrepancy metrics
#' 
#' @param eval_t Length-K vector of times to evaluation K-M curve at
#' @param dat Data frame containing variables yt, dyt, and z
#' @return K x 3 data frame with eval_t, and K-M estimates at each eval_t in
#' the control (KM0) and treated (KM1) arms
#' @export
calculate_km_rep <- function(eval_t, dat) {
  km_fit0 <- survival::survfit(survival::Surv(yt, dyt) ~ 1, 
                               data = dat[dat$z == 0, ])
  survest0 <- stepfun(km_fit0$time, c(1, km_fit0$surv))
  rep_0 <- survest0(eval_t)
  
  km_fit1 <- survival::survfit(survival::Surv(yt, dyt) ~ 1, 
                               data = dat[dat$z == 1, ])
  survest1 <- stepfun(km_fit1$time, c(1, km_fit1$surv))
  rep_1 <- survest1(eval_t)
  res <- data.frame(eval_t = eval_t, KM0 = rep_0, KM1 = rep_1)
  return(res)
}



#' Calculate always-alive fraction in a 
#' 
#' @param eval_t Length-K vector of times for evaluating AA fraction
#' @param dat Data frame containing yt0_imp, yt1_imp, dyt0_imp, dyt1_imp
#' @return K-row data frame with time of evaluation and fraction always-alive
#' (according to KM to account for censoring at the final t)
#' @export
calculate_frac_aa_rep <- function(eval_t, dat) {
  dat$di <- pmin(dat$yt0_imp, dat$yt1_imp)
  dat$deltai <- pmax(dat$dyt0_imp, dat$dyt1_imp)
  km_fit_aa <- survival::survfit(survival::Surv(di, deltai) ~ 1, 
                                 data = dat)
  survest_aa <- stepfun(km_fit_aa$time, c(1, km_fit_aa$surv))
  rep_aa <- survest_aa(eval_t)
  res <- data.frame(eval_t = eval_t, AA = rep_aa)
  return(res)
}


#' Calculate the always-alive discrepancy indicator for a single MCMC
#' iteration (b)
#' 
#' @param eval_t Length-K vector of evaluation times
#' @param frac_aa_obs K x B Result of \code{\link{calculate_frac_aa_rep}} for 
#' the observed data set used to fit the model
#' @param b MCMC iteration number to calculate discrepancy for
#' @param seed Seed for replicate
#' @param stan_fit Stan fit object to extract bth posterior draw from
#' @param z Observed treatment assignment
#' @param x Observed design matrix
#' @param cens_times Censoring times for KM calculation
#' @return K x 2 matrix of booleans for whether the observed values exceed the
#' replicate data KM values for that arm
calculate_frac_aa_disc_b <- function(eval_t, frac_aa_obs, b, seed, stan_fit, z, 
                                     x, cens_times) {
  rep_dat <- make_repdata(b = b, seed = seed, stan_fit = stan_fit, 
                          z = z, x = x, cens_times = cens_times, 
                          observed = FALSE, add_imp = TRUE)
  
  frac_aa_rep <- calculate_frac_aa_rep(eval_t = eval_t, dat = rep_dat)
  # Unlike KM, which is constant across b, compare to bth posterior predicted in observed
  test_frac_aa_disc_b <- (frac_aa_obs[ , b] > frac_aa_rep)[ , -1] # omit eval_t column
  return(test_frac_aa_disc_b)
}


#' Calculate the Kaplan-Meier discrepancy indicator for a single MCMC
#' iteration (b)
#' 
#' @param eval_t Length-K vector of evaluation times
#' @param km_obs Result of \code{\link{calculate_km_rep}} for the
#' observed data set used to fit the model
#' @param b MCMC iteration number to calculate discrepancy for
#' @param seed Seed for replicate
#' @param stan_fit Stan fit object to extract bth posterior draw from
#' @param z Observed treatment assignment
#' @param x Observed design matrix
#' @param cens_times Censoring times for KM calculation
#' @return K x 2 matrix of booleans for whether the observed values exceed the
#' replicate data KM values for that arm
calculate_km_disc_b <- function(eval_t, km_obs, b, seed, stan_fit, z, x, 
                                cens_times) {
  rep_dat <- make_repdata(b = b, seed = seed, stan_fit = stan_fit, 
                          z = z, x = x, cens_times = cens_times)
  
  km_rep <- calculate_km_rep(eval_t = eval_t, dat = rep_dat)
  test_km_disc_b <- (km_obs > km_rep)[ , -1] # omit eval_t column
  return(test_km_disc_b)
}



#' Kaplan-Meier discrepancy test statistic function 
#' 
#' Vectorized along the MCMC iterations ("b")
v_calculate_km_disc_b <- Vectorize(FUN = calculate_km_disc_b, 
                                   vectorize.args = "b",
                                   SIMPLIFY = FALSE)


#' Always-alive discrepancy test statistic function 
#' 
#' Vectorized along the MCMC iterations ("b")
v_calculate_frac_aa_disc_b <- Vectorize(FUN = calculate_frac_aa_disc_b, 
                                        vectorize.args = "b",
                                        SIMPLIFY = TRUE)

#' Calculate matrix of Kaplan-Meier marginal survival discrepancy test
#' statistics based on a model fit
#' 
#' @param eval_t Length-K vector of evaluation times
#' @param stan_res Resulting fit object from \code{\link{run_scr_replicate}}
#' @param cens_times Vector of censoring times
#' @param seed Seed for simulating the posterior predictive data
#' @param subsamp Number of MCMC iterations to sample for the discrepancy 
#' calculation. Default is 500. If \code{subsamp = FALSE}, no subsampling will 
#' be done and the discrepancy measure will be calculated for each MCMC scan.
#' @return K x 2 matrix of test statistics for each arm at each eval_t
#' @export
calculate_km_disc <- function(eval_t, stan_res, cens_times, seed, subsamp = 500) {
  
  # Extract important components
  obs_dat <- stan_res$dat
  stan_fit <- stan_res$stan_fit
  B <- length(extract(stan_fit, par = "lp__")[[1]])
  if (subsamp == FALSE || (subsamp >= B)) {
    b_samp <- 1:B
  } else {
    b_samp <- round(seq(1, B, length.out = subsamp))
  }
  
  # Calculate KM for observed data
  km_obs <- calculate_km_rep(eval_t = eval_t, dat = obs_dat)
  
  # Calculate for the iterations of the discrepancy data
  z <- obs_dat$z
  x <- stan_res$xmat
  all_disc <- v_calculate_km_disc_b(eval_t = eval_t, 
                                    km_obs = km_obs, 
                                    b = b_samp, seed = seed, 
                                    stan_fit, z, x, cens_times)
  
  # Calculate mean of discrepancy test statistics across the iterations
  km_disc <- apply(simplify2array(all_disc), c(1, 2), FUN = mean)
  
  return(km_disc)
}



#' Calculate vector of always-alive discrepancy test statistics based on a 
#' model fit
#' 
#' @param eval_t Length-K vector of evaluation times
#' @param stan_res Resulting fit object from \code{\link{run_scr_replicate}}
#' @param cens_times Vector of censoring times
#' @param seed Seed for simulating the posterior predictive data
#' @param subsamp Number of MCMC iterations to sample for the discrepancy 
#' calculation. Default is 500. If \code{subsamp = FALSE}, no subsampling will 
#' be done and the discrepancy measure will be calculated for each MCMC scan.
#' @return Length-K vector of test statistics for always-alive at each eval_t
#' @export
calculate_frac_aa_disc <- function(eval_t, stan_res, cens_times, seed, subsamp = 500) {
  
  # Extract important components
  obs_dat <- stan_res$dat
  stan_fit <- stan_res$stan_fit
  B <- length(extract(stan_fit, par = "lp__")[[1]])
  if (subsamp == FALSE || (subsamp >= B)) {
    b_samp <- 1:B
  } else {
    b_samp <- round(seq(1, B, length.out = subsamp))
  }
  
  # Calculate fraction always-alive for observed data
  z <- obs_dat$z
  x <- stan_res$xmat
  pp_obs <- posterior_predict_sample(stan_fit = stan_fit, 
                                     yr = obs_dat$yr, 
                                     yt = obs_dat$yt, 
                                     dyr = obs_dat$dyr, 
                                     dyt = obs_dat$dyt, 
                                     z = z, 
                                     xmat = x)[ , , b_samp]
  
  # Stupid hacky thing to make b_samp slicing work later
  frac_aa_obs <- matrix(NA, nrow = length(eval_t), ncol = B)
  frac_aa_obs[ , b_samp] <- apply(pp_obs, MARGIN = 3, 
                                  FUN = v_calculate_frac_aa, eval_t = eval_t)
  
  # Calculate for the iterations of the discrepancy data
  all_disc <- v_calculate_frac_aa_disc_b(eval_t = eval_t, 
                                         frac_aa_obs = frac_aa_obs, 
                                         b = b_samp, seed = seed, 
                                         stan_fit, z, x, cens_times)
  
  # Calculate mean of discrepancy test statistics across the iterations
  frac_aa_disc <- rowMeans(all_disc)
  return(frac_aa_disc)
}
