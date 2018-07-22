#' Return data generation parameters for different scenarios
#' 
#' @param scenario
#' 1 = ignoring death is MISLEADING: treatment has no effect except to kill 
#' people faster so they do not have non-terminal event
#' 2 = stratum are SHARP: covariates are highly predictive of survival
#' 3 = stratum are BLUNT: covariates do almost nothing
#' 4 = DIFFERENTIAL frailty effects: frailty is much less predictive of 
#' non-terminal event than terminal event (raised to a power < 1)
#' 5 = frailty is LOGNORMAL: misspecified frailty distribution
#' @return Named list of data generation parameters
#' @export
return_dgp_parameters <- function(scenario) {
  p <- switch(scenario,
              "1" = list(beta1.true = c(0, 0.1, 0.1, 0.1),
                         beta2.true = c(1, 0.1, 0.1, 0.1),
                         beta3.true = c(-1, 0.1, 0.1, 0.1),
                         alpha1.true = 1,
                         alpha2.true = 0.95,
                         alpha3.true = 1,
                         kappa1.true = 0.003,
                         kappa2.true = 0.002, 
                         kappa3.true = 0.004,
                         theta.true = 0.2),
              "2" = list(beta1.true = c(-0.5, 1, 1.1, 1),
                         beta2.true = c(-0.7, 1, 1, 1.2),
                         beta3.true = c(-0.5, 1.3, 1, 1),
                         alpha1.true = 1,
                         alpha2.true = 0.95,
                         alpha3.true = 1,
                         kappa1.true = 0.003,
                         kappa2.true = 0.002, 
                         kappa3.true = 0.004,
                         theta.true = 0.2),
              "3" = list(beta1.true = c(-0.5, 0.01, 0.01, 0.01),
                         beta2.true = c(-0.7, 0.01, 0.01, 0.01),
                         beta3.true = c(-0.5, 0.01, 0.01, 0.01),
                         alpha1.true = 1,
                         alpha2.true = 0.95,
                         alpha3.true = 1,
                         kappa1.true = 0.003,
                         kappa2.true = 0.002, 
                         kappa3.true = 0.004,
                         theta.true = 0.2),
              "4" = list(beta1.true = c(-0.5, 0.1, 0.1, 0.1, 0.5),
                         beta2.true = c(-0.7, 0.1, 0.1, 0.1, 1),
                         beta3.true = c(-0.5, 0.1, 0.1, 0.1, 1),
                         alpha1.true = 1,
                         alpha2.true = 0.95,
                         alpha3.true = 1,
                         kappa1.true = 0.003,
                         kappa2.true = 0.002, 
                         kappa3.true = 0.004,
                         theta.true = 0),
              "5" = list(beta1.true = c(-0.5, 0.1, 0.1, 0.1, 1),
                         beta2.true = c(-0.7, 0.1, 0.1, 0.1, 1),
                         beta3.true = c(-0.5, 0.1, 0.1, 0.1, 1),
                         alpha1.true = 1,
                         alpha2.true = 0.95,
                         alpha3.true = 1,
                         kappa1.true = 0.003,
                         kappa2.true = 0.002, 
                         kappa3.true = 0.004,
                         theta.true = 0))
  return(p)
}



#' Simulate data for different scenarios
#' 
#' @param n Sample size
#' @param seed Seed
#' @param scenario Scenario; see \code{\link{return_dgp_parameters}}
#' @return Data frame
#' @export
simulate_scenario <- function(n = 5000, seed = 123, scenario) {
  set.seed(seed)
  params <- return_dgp_parameters(scenario = as.character(scenario))
  P <- 3
  z <- rbinom(n, size = 1, prob = 0.5)
  x <- cbind(z, matrix(rnorm(n * P), ncol = 3))
  colnames(x)[2:(P+1)] <- paste0("X", 1:P)
  if (scenario %in% c(4, 5)) {
    if (scenario == 4){
      theta.true <- 0.2
      frailty <- rgamma(n, 1 / theta.true, 1 / theta.true)
      x1 <- x2 <- x3 <- cbind(x, log(frailty))
    } else if (scenario == 5) {
      # lognormal with median of 1 and variance of 0.2
      frailty <- exp(rnorm(n, mean = 0, sd = sqrt(log(0.5 + 0.3 * sqrt(5)))))
      x1 <- x2 <- x3 <- cbind(x, log(frailty))
    } 
  } else {
    x1 <- x2 <- x3 <- x 
  }
  cens_lb <- get_eval_t() * 0.5
  cens_ub <- cens_lb * 6
  dat <- SemiCompRisks::simID(id = NULL, 
                              x1 = x1, x2 = x2, x3 = x3,
                              beta1.true = params$beta1.true,
                              beta2.true = params$beta2.true,
                              beta3.true = params$beta3.true,
                              alpha1.true = params$alpha1.true, 
                              alpha2.true = params$alpha2.true,
                              alpha3.true = params$alpha3.true,
                              kappa1.true = params$kappa1.true, 
                              kappa2.true = params$kappa2.true, 
                              kappa3.true = params$kappa3.true,
                              theta.true = params$theta.true, 
                              SigmaV.true = NULL,
                              cens = c(cens_lb, cens_ub))
  return(cbind(dat, x, z = z))
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
#' @param ... Parameters to pass to Stan via \code{\link{scr_gamma_frailty_stan}}
#' @return Named list of simulated data (dat), common design matrix (xmat), and 
#' stan fit (stan_fit) objects
#' @export
run_scr_replicate <- function(n, seed, scenario, iter = 2000, chains = 4, 
                              sigma_pa = 11, sigma_pb = 2, ...) {
  dat <- simulate_scenario(n = n, seed = seed, scenario = scenario)  
  xmat <- make_xmat_all_X(dat)
  stan_fit <- scr_gamma_frailty_stan(x = xmat, z = dat$z, 
                                     yr = dat$y1, yt = dat$y2, 
                                     dyr = dat$delta1, dyt = dat$delta2,
                                     use_priors = TRUE, 
                                     sigma_pa = sigma_pa, sigma_pb = sigma_pb,
                                     iter = iter, chains = chains,
                                     ...)
  return(list(dat = dat, xmat = xmat, stan_fit = stan_fit))
}

