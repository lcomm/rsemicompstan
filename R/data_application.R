#' Fit basic model in data application with pancreatic data
#' 
#' Prior for sigma is weak (~ 20 prior obs) around frailty variance in Lee 2015.
#' Applies fuzzing to handle discretized time scale. 
#' 
#' @param file File (including path) to RDS with scrData object
#' @param use_priors Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param init Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param init_r Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param sigma_pa Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param sigma_pb Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param iter Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param chains Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param ... Parameters to pass to Stan via \code{\link{scr_gamma_frailty_stan}}
#' @return List of data, design matrix, and Stan fit object
#' @export
fit_data_app <- function(file = "~/Dropbox/Semicompeting-PS/scrData.rds",
                         use_priors = TRUE, 
                         init = "random", 
                         init_r = 0.5,
                         sigma_pa = 21, 
                         sigma_pb = 7.1,
                         iter = 2000, 
                         chains = 4,
                         ...) {
  
  # Read in data
  scr_pc <- readRDS(file = file)[1:3000,]
  
  # Alias variables
  scr_pc$yr <- scr_pc$time1
  scr_pc$yt <- scr_pc$time2
  scr_pc$dyr <- scr_pc$event1
  scr_pc$dyt <- scr_pc$event2
  
  # Fuzz so that events never happen at time zero
  # Add one-half of the minimum observed time scale to sojourn times of zero
  # (so if time measured in days, same-day non-terminal and terminal events are 
  # said to be a half day apart)
  all_unique <- sort(unique(c(scr_pc$yr, scr_pc$yt,
                              (scr_pc$yt - scr_pc$yr)[scr_pc$dyr == 1])))
  fuzz       <- min(diff(all_unique)) / 2
  zero_nt    <- which(scr_pc$yr == 0)
  zero_death <- which(scr_pc$yt == 0)
  scr_pc$yr[zero_nt]    <- scr_pc$yr[zero_nt] + fuzz
  scr_pc$yt[zero_death] <- scr_pc$yt[zero_death] + fuzz
  zero_soj   <- which((scr_pc$yt == scr_pc$yr)[scr_pc$dyr == 1])
  scr_pc$yt[zero_soj]   <- scr_pc$yt[zero_soj] + fuzz
  
  # Adjustment covariates are race, standardized age, sex
  # "Treatment" is being discharged to homecare or hospice (vs. not)
  scr_pc$z <- pmax(scr_pc$disc_hospice, scr_pc$disc_homecare)
  xmat <- cbind(scr_pc$race_noWhite, scr_pc$sex_female)
  
  # Fit model
  stan_fit <- scr_gamma_frailty_stan(x = xmat, z = scr_pc$z, 
                                     yr = scr_pc$yr, yt = scr_pc$yt, 
                                     dyr = scr_pc$dyr, dyt = scr_pc$dyt,
                                     use_priors = use_priors, 
                                     init = init, init_r = init_r,
                                     sigma_pa = sigma_pa, sigma_pb = sigma_pb,
                                     iter = iter, chains = chains,
                                     ...)
  
  # Return
  return(list(dat = scr_pc, xmat = xmat, stan_fit = stan_fit))
}


