#' Fit semicompeting risks with gamma frailty without marginalizing over frailties
#'
#' This version allows frailty importance to vary based on destination.
#'
#' @param x N x P design matrix, no intercept
#' @param z Length-N vector of binary treatment indicators
#' @param yr Length-N vector of non-terminal event times
#' @param yt Length-N vector of terminal event times
#' @param dyr Length-N vector binary indicators for having observed the 
#' non-terminal event
#' @param dyt Length-N vector binary indicators for having observed the 
#' terminal event
#' @param use_priors Whether to use weakly informative/data-driven priors
#' @param sigma_pa Hyperparameter alpha for inverse gamma prior on sigma
#' @param sigma_pb Hyperparameter beta for inverse gamma prior on sigma. Prior 
#' mean for sigma is beta/(alpha - 1) for alpha > 1 and prior mode is 
#' beta/(alpha + 1).
#' @param ... Additional parameters to pass to `rstan::sampling`
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @export
scr_gamma_frailty_unmarginalized_variable_stan <- function(x, z, yr, yt, dyr, dyt, use_priors = TRUE, 
                                   sigma_pa = 0.7, sigma_pb = 0.7, ...) {
  if (use_priors) {
    pm <- make_prior_means(yr = yr, yt = yt, dyr = dyr, dyt = dyt)  
  } else {
    pm <- list(log_alpha_pmean = rep(0, 6), log_kappa_pmean = rep(0, 6))
  }
  out <- rstan::sampling(stanmodels$scr_gamma_frailty_unmarginalized_variable,
                         data = list(N = NROW(x), 
                                     z = z,
                                     P = NCOL(x),
                                     yr = yr,
                                     yt = yt,
                                     dyr = dyr,
                                     dyt = dyt,
                                     use_priors = use_priors * 1,
                                     log_alpha_pmean = pm$log_alpha_pmean,
                                     log_kappa_pmean = pm$log_kappa_pmean,
                                     sigma_pa = sigma_pa,
                                     sigma_pb = sigma_pb),
                         ...)
  return(out)
  
}
