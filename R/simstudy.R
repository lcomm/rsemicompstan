#'
#' Very lightly adapted from SemiCompRisks::simID function. Only difference
#' is that you can pass in censoring times 
#' 
#' 
modified_simID <- function(id = NULL, x1, x2, x3, 
                           beta1.true, beta2.true, beta3.true, 
                           alpha1.true, alpha2.true, alpha3.true, 
                           kappa1.true, kappa2.true, kappa3.true, 
                           theta.true, SigmaV.true = NULL, 
                           cens, cens_times = NULL) {
  if (!is.null(id) & is.null(SigmaV.true)) {
    stop("SigmaV.true must be given to simulate correlated data")
  }
  else {
    n <- dim(x1)[1]
    p1 <- dim(x1)[2]
    p2 <- dim(x2)[2]
    p3 <- dim(x3)[2]
    if (theta.true > 0) {
      gamma.true <- rgamma(n, 1/theta.true, 1/theta.true)
    }
    if (theta.true == 0) {
      gamma.true <- rep(1, n)
    }
    if (is.null(id)) {
      LP1 <- as.vector(beta1.true %*% t(x1))
      LP2 <- as.vector(beta2.true %*% t(x2))
      LP3 <- as.vector(beta3.true %*% t(x3))
    }
    if (!is.null(id)) {
      J <- length(unique(id))
      nj <- as.vector(table(id))
      Vmat <- mvrnorm(J, rep(0, 3), SigmaV.true)
      LP1 <- as.vector(beta1.true %*% t(x1) + rep(Vmat[, 
                                                       1], nj))
      LP2 <- as.vector(beta2.true %*% t(x2) + rep(Vmat[, 
                                                       2], nj))
      LP3 <- as.vector(beta3.true %*% t(x3) + rep(Vmat[, 
                                                       3], nj))
    }
    Rind <- NULL
    R <- rweibull(n, shape = alpha1.true, scale = exp(-(log(kappa1.true) + 
                                                          LP1 + log(gamma.true))/alpha1.true))
    D <- rweibull(n, shape = alpha2.true, scale = exp(-(log(kappa2.true) + 
                                                          LP2 + log(gamma.true))/alpha2.true))
    yesR <- R < D
    if (any(yesR)) {
      D[yesR] <- R[yesR] + rweibull(sum(yesR), shape = alpha3.true, 
                                    scale = exp(-(log(kappa3.true) + LP3[yesR] + log(gamma.true[yesR]))/alpha3.true))
    }
    delta1 <- rep(NA, n)
    delta2 <- rep(NA, n)
    y1 <- R
    y2 <- D
    if (is.null(cens_times)) {
      Cen <- runif(n, cens[1], cens[2])  
    } else {
      stopifnot(length(cens_times) == n)
      Cen <- cens_times
    }
    ind01 <- which(D < R & D < Cen)
    y1[ind01] <- D[ind01]
    delta1[ind01] <- 0
    delta2[ind01] <- 1
    ind10 <- which(R < D & R < Cen & D >= Cen)
    y2[ind10] <- Cen[ind10]
    delta1[ind10] <- 1
    delta2[ind10] <- 0
    ind00 <- which(R >= Cen & D >= Cen)
    y1[ind00] <- Cen[ind00]
    y2[ind00] <- Cen[ind00]
    delta1[ind00] <- 0
    delta2[ind00] <- 0
    ind11 <- which(R < Cen & D < Cen & R < D)
    delta1[ind11] <- 1
    delta2[ind11] <- 1
    ret <- data.frame(cbind(y1, delta1, y2, delta2))
    return(ret)
  }
}
  


#' Return data generation parameters for different scenarios
#' 
#' @param scenario
#' 1 = CORRECT specification scenario;
#' 2 = LOGNORMAL FRAILTY; frailties are lognormally distributed with same 
#' variance as in (1)
#' 3 = UNEQUAL FRAILTY COEFFICIENTS; frailties are gamma-distributed with same 
#' variance as in (1) but frailty has different exponent in every model
#' @return Named list of data generation parameters
#' @export
return_dgp_parameters <- function(scenario) {
  
  # Set reference rates
  strong_sigma <- 1.2
  control <- list(beta1 = c(0.34402, -0.14652, 0.35407, 0.33832, -0.11605),
                  beta2 = c(-0.22843, -0.10205, 0.43752, 0.58252, 0.33016),
                  beta3 = c(0.11858, -0.10549, 0.46426, 0.55319, 0.07844),
                  alpha1 = 0.96949,
                  alpha2 = 1.19239,
                  alpha3 = 0.91350,
                  kappa1 = 0.01591,
                  kappa2 = 0.00322,
                  kappa3 = 0.00857)
  treated <- list(beta1 = c(0.12413, -0.48149, 0.64025, 0.24664, -0.17149),
                  beta2 = c(-0.38046, -0.57772, 0.6031, 0.57706, 0.21149),
                  beta3 = c(0.16526, -0.1467, 0.46675, 0.35684, 0.15129),
                  alpha1 = 1.04497,
                  alpha2 = 1.30856,
                  alpha3 = 0.94227,
                  kappa1 = 0.0147,
                  kappa2 = 0.0041,
                  kappa3 = 0.01352)
  
  # Frailty variance
  sigma <- 1.44256
  
  # Frailty importances (not used in model fit but used to induce misspecification)
  control$omega1 <- control$omega2 <- control$omega3 <- 1
  treated$omega1 <- treated$omega2 <- treated$omega3 <- 1
  
  # Make basic lists
  p1 <- list(treated = treated, control = control, sigma = sigma, 
             distn = "gamma")
  p2 <- list(treated = treated, control = control, sigma = sigma, 
             distn = "lognormal")
  p3 <- list(treated = treated, control = control, sigma = sigma, 
             distn = "gamma")
  
  # Overwrite frailty coeffients for second misspecification type
  p3$control$omega1 <- exp(-0.1) 
  p3$control$omega2 <- exp(-0.2)
  p3$control$omega3 <- exp(0)
  p3$treated$omega1 <- exp(0.15)
  p3$treated$omega2 <- exp(0.2)
  p3$treated$omega3 <- exp(0)
  
  p <- switch(scenario,
              "1" = p1,
              "2" = p2,
              "3" = p3)
  return(p)
}



#' Simulate data for based on a parameter output list
#' 
#' @param n Sample size
#' @param seed Seed
#' @param params Parameter list like from \code{\link{return_dgp_parameters}}
#' @param data_match (Optional) list containing \code{z}, \code{x}, and 
#' \code{cens_times} to match; useful for simulating replicate data sets for 
#' discrepancy measures
#' @param censor Whether to impose censoring. Defaults to \code{TRUE}. (If 
#' \code{FALSE}, applies uniform censoring at starting at 2 * 99.999 percentile 
#' to ensure essentially no censoring takes place.)
#' @param cens_times Length-n vector of right censoring times. Ignored if 
#' \code{data_match} is not \code{NULL}.
#' @param observed Whether to only return data set of observed potential outcomes
#' or to include both sets. Defaults to \code{TRUE}.
#' @param add_imp Whether to add \code{_imp} suffix to the potential outcomes
#' @return Data frame
#' @export
simulate_from_param <- function(n = 5000, seed = 123, params, data_match = NULL,
                                censor = TRUE, cens_times = NULL,
                                observed = TRUE, add_imp = FALSE) {
  set.seed(seed)
  P <- length(params$treated$beta1)
  if (is.null(data_match)) {
    z <- rbinom(n, size = 1, prob = 0.5)
    if (P < 5) {
      x <- matrix(rnorm(n * P), ncol = P)  
    } else {
      x_bin <- matrix(rbinom(n = n * 4, size = 1, p = 0.5), nrow = n)
      x_con <- matrix(rnorm(n = n * (P - 4)), nrow = n)
      x <- scale(cbind(x_bin, x_con), center = TRUE, scale = FALSE)
    }
    colnames(x) <- paste0("X", 1:P)
  } else {
    z <- data_match$z
    x <- data_match$x
    cens_times <- data_match$cens_times
  }
  
  frailty <- simulate_frailty(n = n, sigma = params$sigma, 
                              distn = params$distn)
  x1 <- x2 <- x3 <- cbind(x, log(frailty))
  if (censor) {
    cens_lb <- get_eval_t() * 0.5
    cens_ub <- cens_lb * 6  
  } else {
    cens_lb <- max(qweibull(p = 0.99999, 
                            shape = params$treated$alpha1, 
                            scale = exp(-(log(params$treated$kappa1) +
                                            x1 %*% c(params$treated$beta1, 
                                                     params$treated$omega1)) / 
                                          params$treated$alpha1)),
                   qweibull(p = 0.99999, 
                            shape = params$treated$alpha2, 
                            scale = exp(-(log(params$treated$kappa2) +
                                            x1 %*% c(params$treated$beta2, 
                                                     params$treated$omega2)) / 
                                          params$treated$alpha2)),
                   qweibull(p = 0.99999, 
                            shape = params$control$alpha1, 
                            scale = exp(-(log(params$control$kappa1) +
                                            x1 %*% c(params$control$beta1,
                                                     params$control$omega1)) / 
                                          params$control$alpha1)),
                   qweibull(p = 0.99999, 
                            shape = params$control$alpha2, 
                            scale = exp(-(log(params$control$kappa2) +
                                            x1 %*% c(params$control$beta2,
                                                     params$control$omega2)) / 
                                          params$control$alpha2))) * 2
    cens_ub <- 2 * cens_lb
  }
  treated <- modified_simID(id = NULL, 
                            x1 = x1, x2 = x2, x3 = x3,
                            beta1.true = c(params$treated$beta1, 
                                           params$treated$omega1),
                            beta2.true = c(params$treated$beta2,
                                           params$treated$omega2),
                            beta3.true = c(params$treated$beta3,
                                           params$treated$omega3),
                            alpha1.true = params$treated$alpha1, 
                            alpha2.true = params$treated$alpha2,
                            alpha3.true = params$treated$alpha3,
                            kappa1.true = params$treated$kappa1, 
                            kappa2.true = params$treated$kappa2, 
                            kappa3.true = params$treated$kappa3,
                            theta.true = 0,
                            SigmaV.true = NULL,
                            cens = c(cens_lb, cens_ub),
                            cens_times = cens_times)
  control <- modified_simID(id = NULL, 
                            x1 = x1, x2 = x2, x3 = x3,
                            beta1.true = c(params$control$beta1,
                                           params$control$omega1),
                            beta2.true = c(params$control$beta2,
                                           params$control$omega2),
                            beta3.true = c(params$control$beta3,
                                           params$control$omega3),
                            alpha1.true = params$control$alpha1, 
                            alpha2.true = params$control$alpha2,
                            alpha3.true = params$control$alpha3,
                            kappa1.true = params$control$kappa1, 
                            kappa2.true = params$control$kappa2, 
                            kappa3.true = params$control$kappa3,
                            theta.true = 0,
                            SigmaV.true = NULL,
                            cens = c(cens_lb, cens_ub),
                            cens_times = cens_times)
  dat <- control * NA
  dat[z == 0, ] <- control[z == 0, ]
  dat[z == 1, ] <- treated[z == 1, ]
  po_names <- c("yr", "dyr", "yt", "dyt")
  if (add_imp) {
    names(dat) <- paste0(po_names, "_imp")
  } else {
    names(dat) <- po_names
  }
  if (observed) {
    return(cbind(dat, x, z = z))  
  } else {
    if (add_imp) {
      colnames(control) <- paste0(po_names, "0_imp")
      colnames(treated) <- paste0(po_names, "1_imp")
    } else {
      colnames(control) <- paste0(po_names, "0")
      colnames(treated) <- paste0(po_names, "1")
    }
    return(cbind(dat, control, treated, x, z = z, frailty = frailty))
  }
}



#' Simulate data for different scenarios
#' 
#' @param n Sample size
#' @param seed Seed
#' @param scenario Scenario; see \code{\link{return_dgp_parameters}}
#' @param censor Whether to impose censoring. Defaults to TRUE. (If FALSE, 
#' applies uniform censoring at starting at 2 * 99.999 percentile to ensure
#' essentially no censoring takes place.)
#' @param cens_times Vector of censoring times to impose. Ignored if 
#' \code{CENSOR = FALSE} and overridden if \code{data_match = TRUE}.
#' @param observed Whether to only return data set of observed potential outcomes
#' or to include both sets. Defaults to TRUE.
#' @param add_imp Whether to add _imp suffix to the potential outcomes
#' @return Data frame
#' @export
simulate_scenario <- function(n = 5000, seed = 123, scenario, censor = TRUE, 
                              cens_times = NULL,
                              observed = TRUE, add_imp = FALSE) {
  params <- return_dgp_parameters(scenario = as.character(scenario))
  res <- simulate_from_param(n = n, seed = seed, params = params, censor = censor, 
                             cens_times = cens_times,
                             observed = observed, add_imp = add_imp)
  return(res)
}



#' Simulate "true" data set for operating characteristics check
#' 
#' @param scenario Scenario; 
#' @param add_imp Whether to add "_imp" suffix to potential outcomes, as is 
#' needed for the calculation of causal effects
#' @return Simulated data set
#' @export
simulate_truth_scenario <- function(scenario, add_imp = FALSE) {
  dat <- simulate_scenario(n = 100000, seed = 313517734, scenario,
                           censor = FALSE, observed = FALSE)
  return(dat)
}



#' Run one replicate of a semicompeting risks scenario
#' 
#' @param n Sample size
#' @param seed Random number seed
#' @param scenario Scenario; see \code{\link{return_dgp_parameters}}
#' @param iter Number of MCMC iterations
#' @param chains Number of MCMC chains
#' @param sigma_pa Hyperparameter alpha for inverse gamma prior on sigma
#' @param sigma_pb Hyperparameter beta for inverse gamma prior on sigma. Prior 
#' mean for sigma is beta/(alpha - 1) for alpha > 1 and prior mode is 
#' beta/(alpha + 1).
#' @param init Chain initialization type ("0" or "random")
#' @param init_r Range for random starting values. Default is 0.5 for (-0.5, 0.5) 
#' initialization instead the normal (-2, 2) for fewer convergence issues.
#' @param cens_time Time to apply common administrative right-censoring 
#' (default is to have random censoring)
#' @param mc.cores Number of cores to run chains on. Default is 1.
#' @param ... Parameters to pass to Stan via \code{\link{scr_gamma_frailty_stan}}
#' @return Named list of simulated data (dat), common design matrix (xmat), and 
#' stan fit (stan_fit) objects
#' @export
run_scr_replicate <- function(n, seed, scenario, iter = 2000, chains = 4, 
                              sigma_pa = 11, sigma_pb = 2, 
                              init = "random", init_r = 0.5, 
                              cens_time = NULL,
                              mc.cores = 1,
                              ...) {
  if (!is.null(cens_time)) {
    dat <- simulate_scenario(n = n, seed = seed, scenario = scenario, 
                             censor = TRUE, cens_times = rep(cens_time, n),
                             observed = TRUE, add_imp = FALSE)
  } else {
    dat <- simulate_scenario(n = n, seed = seed, scenario = scenario, 
                             censor = TRUE, observed = TRUE, add_imp = FALSE)
  }
  xmat <- make_xmat_all_X(dat)
  stan_fit <- scr_gamma_frailty_stan(x = xmat, z = dat$z, 
                                     yr = dat$yr, yt = dat$yt, 
                                     dyr = dat$dyr, dyt = dat$dyt,
                                     use_priors = TRUE, 
                                     sigma_pa = sigma_pa, sigma_pb = sigma_pb,
                                     iter = iter, chains = chains,
                                     mc.core = mc.cores,
                                     ...)
  return(list(dat = dat, xmat = xmat, stan_fit = stan_fit))
}



#' Add the _imp suffix to potential outcome variables in a data set
#' 
#' Useful for calculating causal effects in ground truth data sets
#' 
#' @param dat Data frame containing variables named yr0, yr1, ... dyt1
#' @return Data frame with renamed columns
#' @export
add_imp_suffix <- function(dat) {
  po_names <- apply(expand.grid(c("d",""), "y", c("r","t"), c("0","1")),
                    1, paste0, collapse = "")
  po_names <- colnames(dat)[colnames(dat) %in% po_names]
  which_po <- which(colnames(dat) %in% po_names)
  po_suff <- paste0(po_names, "_imp")
  colnames(dat)[which_po] <- po_suff
  return(dat)
}



#' Simulate a large data set to be the "truth" for the scenario estimands
#' 
#' @param scenario Scenario; see \code{\link{return_dgp_parameters}}
#' @param add_imp Whether to add _imp suffix to the potential outcomes
#' @return Complete data set of potential outcomes
#' @export
simulate_scenario_truth <- function(scenario = 1, add_imp = FALSE) {
  dat <- simulate_scenario(n = 10^5, seed = 42313, scenario = scenario,
                           censor = FALSE, observed = FALSE, add_imp = add_imp)
  return(dat)
}



#' Make a table of true TV-SACE, RM-SACE, and fraction alive for simulation
#' study scenarios
#' 
#' @param scenarios Vector of simulation scenarios to summarize. Default is 1 
#' through 3.
#' @param eval_t Time points at which to evaluate effects. Default is 30 and 90.
#' @return Data frame of truths (obtained via simulation of a large data set)
#' @export
summarize_scenario_truths <- function(scenarios = 1:3, eval_t = c(30, 90)) {
  true_dat <- list()
  truths <- expand.grid(eval_t = eval_t, scenario = scenarios, 
                        rm_sace = NA, tv_sace = NA, frac_aa = NA)
  for (s in scenarios) {
    true_dat[[s]] <- simulate_scenario_truth(scenario = s, add_imp = TRUE)
    truths$tv_sace[truths$scenario == s] <- 
      v_calculate_tv_sace(eval_t = eval_t, pp = true_dat[[s]])
    truths$rm_sace[truths$scenario == s] <- 
      v_calculate_rm_sace(eval_t = eval_t, pp = true_dat[[s]])
    truths$frac_aa[truths$scenario == s] <- 
      v_calculate_frac_aa(eval_t = eval_t, pp = true_dat[[s]])
  }
  return(truths)  
}


#' Do posterior predictions from a \code{\link{run_scr_replicate}} result
#' 
#' @param rl List output from \code{\link{run_scr_replicate}}
#' @return 3d array of posterior predictions from post-warmup iterations
#' @export
pp_from_result_list <- function(rl) {
  pp <- posterior_predict_sample(stan_fit = rl$stan_fit, 
                                 yr = rl$dat$yr,
                                 yt = rl$dat$yt, 
                                 dyr = rl$dat$dyr, 
                                 dyt = rl$dat$dyt, 
                                 z = rl$dat$z,
                                 xmat = rl$xmat)
  return(pp)
}



#' Process simulation replicate result list to get estimates and credible 
#' intervals for TV-SACE, RM-SACE, and fraction always-alive at a series of
#' \code{eval_t} time points
#' 
#' @param rl Result list from \code{\link{run_scr_replicate}}
#' @param eval_t Times to evaluate the estimators. Default is time 30 and 90.
#' @param alpha Desired alpha for the credible intervals. Default is \code{0.05} 
#' for 95\% intervals
#' @return Named list of estimates (\code{estimates}) and credible intervals
#' (\code{cis})
#' @export
process_replicate_rl <- function(rl, eval_t = c(30, 90), alpha = 0.05) {
  
  stopifnot(0 < alpha, alpha < 1)
  
  # Do posterior predictions
  pp <- pp_from_result_list(rl)
  
  # Calculate quantities at every post-warmup iteration
  tv_sace <- t(apply(pp, MARGIN = 3, FUN = v_calculate_tv_sace, eval_t = eval_t))
  rm_sace <- t(apply(pp, MARGIN = 3, FUN = v_calculate_rm_sace, eval_t = eval_t))
  frac_aa <- t(apply(pp, MARGIN = 3, FUN = v_calculate_frac_aa, eval_t = eval_t))
  
  # Construct named list of estimates (means) and confidence intervals
  estimates <- list(tv_sace = colMeans(tv_sace),
                    rm_sace = colMeans(rm_sace),
                    frac_aa = colMeans(frac_aa))
  ps <- c(alpha / 2, 1 - alpha / 2)
  cis <- list(tv_sace = t(apply(tv_sace, MARGIN = 2, FUN = quantile, p = ps)),
              rm_sace = t(apply(rm_sace, MARGIN = 2, FUN = quantile, p = ps)),
              frac_aa = t(apply(frac_aa, MARGIN = 2, FUN = quantile, p = ps)))
  return(list(estimates = estimates, cis = cis))
}


#' Check whether the CI for a given quantity covers the truth at various t
#' 
#' @param r Result list from \code{\link{process_replicate_rl}}
#' @param truth Truth data set from \code{\link{summarize_scenario_truths}}
#' @param name Character name of estimand (\code{tv_sace}, \code{rm_sace}, or 
#' \code{frac_aa})
#' @export
check_ci_cover_all_t <- function(r, truth, name) {
  covers <- pmin((truth[[name]] > r$cis[[name]][, 1]),
                 (truth[[name]] < r$cis[[name]][, 2]))
  return(covers)
}



#' Calculate bias for an estimate at various t
#' 
#' @param r Result list from \code{\link{process_replicate_rl}}
#' @param truth Truth data set from \code{\link{summarize_scenario_truths}}
#' @param name Character name of estimand (\code{tv_sace}, \code{rm_sace}, or 
#' \code{frac_aa})
#' @export
calculate_bias_all_t <- function(r, truth, name) {
  bias <- r$estimates[[name]] - truth[[name]]
  return(bias)
}



#' Check whether the CI for a given quantity covers the truth at various t
#' 
#' @param r Result list from \code{\link{process_replicate_rl}}
#' @param name Character name of estimand (\code{tv_sace}, \code{rm_sace}, or 
#' \code{frac_aa})
#' @export
calculate_ci_width_all_t <- function(r, name) {
  width <- r$cis[[name]][, 2] - r$cis[[name]][, 1]
  return(width)
}



#' Process a replicate result from \code{\link{run_scr_replicate}} and calculate
#' the operating characteristics of TV-SACE, RM-SACE, and proportion always-alive
#' at a set of time points \code{eval_t}
#' 
#' @param r Result output from \code{\link{run_scr_replicate}}
#' @param truths Truth data set from \code{\link{summarize_scenario_truths}}
#' @param eval_t Vector of time points for estimator evaluation
#' @return Named list of estimates, credible intervals, bias, CI coverage, and 
#' CI width
#' @export
get_replicate_oc <- function(r, truths, scenario, eval_t = c(30, 90)) {
  r_p <- process_replicate_rl(rl = r, eval_t = eval_t)
  truth <- truths[truths$scenario == scenario, ]
  print(truths)
  shell <- data.frame(eval_t = eval_t, tv_sace = NA, rm_sace = NA, frac_aa = NA)
  oc <- list(estimates = shell, cis = shell, bias = shell, 
             ci_coverage = shell, ci_width = shell)
  oc$estimates <- as.data.frame(r_p$estimates)
  oc$cis <- as.data.frame(r_p$cis)
  for (name in names(shell)[-1]) {
    oc$ci_coverage[[name]] <- check_ci_cover_all_t(r = r_p, 
                                                   truth = truth, 
                                                   name = name)
    oc$ci_width[[name]] <- calculate_ci_width_all_t(r = r_p, name = name)
    oc$bias[[name]] <- calculate_bias_all_t(r = r_p, truth = truth, name = name)
  }
  return(oc)
}



#' Function to simulate and run a simulation study scenario as well as calulate
#' the quantities for operating characteristics
#' 
#' @param scenario Simulation scenario for \code{\link{return_dgp_parameters}}
#' @param truths Truth data set from \code{\link{summarize_scenario_truths}}
#' @param eval_t Vector of time points for estimator evaluation
#' @param ... Parameters to pass to \code{\link{run_scr_replicate}}
#' @return Named list of result (\code{r}) and operating characteristics 
#' (\code{oc})
#' @export
run_replicate_with_oc <- function(scenario, truths, eval_t = c(30, 90), ...) {

  stopifnot(scenario %in% truths$scenario)
  for (t_value in eval_t) {
    stopifnot(t_value %in% truths$eval_t)
  }
  r <- run_scr_replicate(scenario = scenario, ...) 
  oc <- get_replicate_oc(r, truths, scenario, eval_t)
  return(list(result = r, oc = oc))
}
