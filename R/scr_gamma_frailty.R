#' Fit semicompeting risks with gamma frailty
#'
#' @param x N x P design matrix, no intercept
#' @param z Length-N vector of binary treatment indicators
#' @param yr Length-N vector of non-terminal event times
#' @param yr Length-N vector of terminal event times
#' @param dyr Length-N vector binary indicators for having observed the 
#' non-terminal event
#' @param use_priors Whether to use weakly informative/data-driven priors
#' @param dyt Length-N vector binary indicators for having observed the 
#' terminal event
#' @param v_precision Prior variance for precision parameter. Default is 1 / 0.7
#' @param ... Additional parameters to pass to `rstan::sampling`
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @examples 
#' \dontrun{
#' rstan_options(auto_write = TRUE)
#' options(mc.cores = 4)
#' library("rsemicompstan")
#' set.seed(123)
#' N <- 5000
#' x1 <-matrix(rnorm(N), ncol = 1)
#' dat <- SemiCompRisks::simID(x1 = x1, x2 = x1, x3 = x1, 
#'                             beta1.true = 0.1, beta2.true = 0.2, beta3.true = 0.3, 
#'                             alpha1.true = 1, alpha2.true = 0.95, alpha3.true = 1,
#'                             kappa1.true = 0.2, kappa2.true = 0.3, kappa3.true = 0.4,
#'                             theta.true = 0.5, SigmaV.true = NULL,
#'                             cens = c(0.5, 10))
#' z <- rbinom(N, size = 1, prob = 0.5)
#' resg <- scr_gamma_frailty_stan(x = x1, z = z, yr = dat$y1, yt = dat$y2,
#'                                dyr = dat$delta1, dyt = dat$delta2,
#'                                use_priors = TRUE,
#'                                sigma_pa = 0.6, sigma_pb = 0.6,
#'                                iter = 2000, chains = 4)
#' }
#' @export
scr_gamma_frailty_stan <- function(x, z, yr, yt, dyr, dyt, use_priors = TRUE, 
                                   sigma_pa = 0.7, sigma_pb = 0.7, ...) {
  if (use_priors) {
    pm <- make_prior_means(yr = yr, yt = yt, dyr = dyr, dyt = dyt)  
  } else {
    pm <- list(log_alpha_pmean = rep(0, 6), log_kappa_pmean = rep(0, 6))
  }
  out <- rstan::sampling(stanmodels$scr_gamma_frailty,
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
