#' Return data generation parameters for different scenarios
#' 
#' @param scenario
#' 1 = REFERENCE scenario: moderate sharpness, moderate frailty variance
#' 2 = stratum are SHARP: covariates are highly predictive of survival, 
#' moderate frailty variance
#' 3 = stratum are BLUNT: covariates do nothing, moderate frailty variance
#' 4 = frailty STRONG: covariates moderately predictive
#' 5 = frailty (virtually) NONEXISTENT: covariates moderately predictive
#' 6 = non-terminal treatment hazard ratio CHANGES: dramatically different shape
#' parameters for non-terminal event hazard, exact same terminal event 
#' processes, moderate frailty variance
#' 7 = non-terminal treatment hazard ratio CONSTANT: exact same non-terminal 
#' event processes, but systematically higher terminal event hazards in controls
#' @param P Number of baseline covariates
#' @return Named list of data generation parameters
#' @export
return_dgp_parameters <- function(scenario, P = 4) {
  
  # Set reference rates
  sgn   <- (-1)^(1:P)
  blunt <- sgn * 0
  mod   <- sgn * (1.2)^(1 / (qnorm(0.9) - qnorm(0.1)))
  sharp <- sgn * (2.5)^(1 / (qnorm(0.9) - qnorm(0.1)))
  mod_sigma <- 0.2
  strong_sigma <- 0.7
  tiny_sigma <- 0.00001
  control <- treated <- list(beta1 = mod, beta2 = mod, beta3 = mod,
                             alpha1 = 1.1, alpha2 = 0.95, alpha3 = 1.15,
                             kappa1 = 0.003, kappa2 = 0.002, kappa3 = 0.004)
  treated$kappa1 <- control$kappa1 * exp(-0.25) # HR is 0.78, with tx having lower rate
  
  # Make modifications
  p1 <- list(treated = treated, control = control, sigma = mod_sigma)
  p2 <- list(treated = treated, control = control, sigma = mod_sigma)
  p3 <- list(treated = treated, control = control, sigma = mod_sigma)
  p4 <- list(treated = treated, control = control, sigma = strong_sigma)
  p5 <- list(treated = treated, control = control, sigma = mod_sigma)
  p6 <- list(treated = treated, control = control, sigma = mod_sigma)
  p7 <- list(treated = treated, control = control, sigma = mod_sigma)
  p2$treated$beta1 <- p2$treated$beta2 <- p2$treated$beta3 <- sharp 
  p2$control$beta1 <- p2$control$beta2 <- p2$control$beta3 <- sharp
  p3$treated$beta1 <- p3$treated$beta2 <- p3$treated$beta3 <- blunt 
  p3$control$beta1 <- p3$control$beta2 <- p3$control$beta3 <- blunt
  p5$sigma <- tiny_sigma
  p6$treated$alpha1 <- 0.8
  p6$control$alpha1 <- 1.2
  p7$treated$kappa1 <- p7$control$kappa1 # undo protective effect for nonterminal
  p7$control$kappa2 <- 2 * p7$treated$kappa2 # implement protective effect for terminal
  
  p <- switch(scenario,
              "1" = p1,
              "2" = p2,
              "3" = p3,
              "4" = p4,
              "5" = p5,
              "6" = p6,
              "7" = p7)
  return(p)
}



#' Simulate data for different scenarios
#' 
#' @param n Sample size
#' @param seed Seed
#' @param scenario Scenario; see \code{\link{return_dgp_parameters}}
#' @param censor Whether to impose censoring. Defaults to TRUE. (If FALSE, 
#' applies uniform censoring at starting at 2 * 99.999 percentile to ensure
#' essentially no censoring takes place.)
#' @param observed Whether to only return data set of observed potential outcomes
#' or to include both sets. Defaults to TRUE.
#' @param add_imp Whether to add _imp suffix to the potential outcomes
#' @return Data frame
#' @export
simulate_scenario <- function(n = 5000, seed = 123, scenario, censor = TRUE, 
                              observed = TRUE, add_imp = FALSE) {
  set.seed(seed)
  P <- 4
  params <- return_dgp_parameters(scenario = as.character(scenario), P = P)
  z <- rbinom(n, size = 1, prob = 0.5)
  x <- matrix(rnorm(n * P), ncol = P)
  colnames(x) <- paste0("X", 1:P)
  frailty <- rgamma(n, 1 / params$sigma, 1 / params$sigma)
  x1 <- x2 <- x3 <- cbind(x, log(frailty))
  if (censor) {
    cens_lb <- get_eval_t() * 0.5
    cens_ub <- cens_lb * 6  
  } else {
    cens_lb <- max(qweibull(p = 0.99999, 
                            shape = params$treated$alpha1, 
                            scale = exp(-(log(params$treated$kappa1) +
                                    x1 %*% c(params$treated$beta1, 1)) / 
                                    params$treated$alpha1)),
                   qweibull(p = 0.99999, 
                            shape = params$treated$alpha2, 
                            scale = exp(-(log(params$treated$kappa2) +
                                    x1 %*% c(params$treated$beta2, 1)) / 
                                    params$treated$alpha2)),
                   qweibull(p = 0.99999, 
                            shape = params$control$alpha1, 
                            scale = exp(-(log(params$control$kappa1) +
                                            x1 %*% c(params$control$beta1, 1)) / 
                                          params$control$alpha1)),
                   qweibull(p = 0.99999, 
                            shape = params$control$alpha2, 
                            scale = exp(-(log(params$control$kappa2) +
                                            x1 %*% c(params$control$beta2, 1)) / 
                                          params$control$alpha2))) * 2
    cens_ub <- 2 * cens_lb
  }
  treated <- SemiCompRisks::simID(id = NULL, 
                                  x1 = x1, x2 = x2, x3 = x3,
                                  beta1.true = c(params$treated$beta1, 1),
                                  beta2.true = c(params$treated$beta2, 1),
                                  beta3.true = c(params$treated$beta3, 1),
                                  alpha1.true = params$treated$alpha1, 
                                  alpha2.true = params$treated$alpha2,
                                  alpha3.true = params$treated$alpha3,
                                  kappa1.true = params$treated$kappa1, 
                                  kappa2.true = params$treated$kappa2, 
                                  kappa3.true = params$treated$kappa3,
                                  theta.true = 0,
                                  SigmaV.true = NULL,
                                  cens = c(cens_lb, cens_ub))
  control <- SemiCompRisks::simID(id = NULL, 
                                  x1 = x1, x2 = x2, x3 = x3,
                                  beta1.true = c(params$control$beta1, 1),
                                  beta2.true = c(params$control$beta2, 1),
                                  beta3.true = c(params$control$beta3, 1),
                                  alpha1.true = params$control$alpha1, 
                                  alpha2.true = params$control$alpha2,
                                  alpha3.true = params$control$alpha3,
                                  kappa1.true = params$control$kappa1, 
                                  kappa2.true = params$control$kappa2, 
                                  kappa3.true = params$control$kappa3,
                                  theta.true = 0,
                                  SigmaV.true = NULL,
                                  cens = c(cens_lb, cens_ub))
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
#' @param ... Parameters to pass to Stan via \code{\link{scr_gamma_frailty_stan}}
#' @return Named list of simulated data (dat), common design matrix (xmat), and 
#' stan fit (stan_fit) objects
#' @export
run_scr_replicate <- function(n, seed, scenario, iter = 2000, chains = 4, 
                              sigma_pa = 11, sigma_pb = 2, 
                              init = "random", init_r = 0.5, 
                              ...) {
  dat <- simulate_scenario(n = n, seed = seed, scenario = scenario, 
                           censor = TRUE, observed = TRUE, add_imp = FALSE)  
  xmat <- make_xmat_all_X(dat)
  stan_fit <- scr_gamma_frailty_stan(x = xmat, z = dat$z, 
                                     yr = dat$yr, yt = dat$yt, 
                                     dyr = dat$dyr, dyt = dat$dyt,
                                     use_priors = TRUE, 
                                     sigma_pa = sigma_pa, sigma_pb = sigma_pb,
                                     iter = iter, chains = chains,
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



