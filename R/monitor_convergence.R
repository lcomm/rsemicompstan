# Utility functions lightly adapted from 
# https://github.com/betanalpha/knitr_case_studies/blob/460dcacbcd2e1d4709515db1fdfc4bbe20bb7bae/principled_bayesian_workflow/stan_utility.R



#' Check for post-warmup divergences
#' 
#' @param stan_fit Stan fit object
#' @param quiet Whether or not to print if divergences are found
#' @return If quiet, return TRUE if problematic
#' @export
check_div <- function(stan_fit, quiet = FALSE) {
  sampler_params <- rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)
  divergent <- do.call(rbind, sampler_params)[, "divergent__"]
  diverged <- sum(divergent)
  iter <- length(divergent)
  div_pct <- diverged / iter * 100
  
  if (!quiet) {
    print(sprintf("%s of %s post-warmup iterations ended with a divergence (%s%%)",
                  diverged, iter, div_pct))
  }
  
  if (diverged == 0) {
    if (!quiet) print("No divergences after warmup")
    if (quiet) return(FALSE)
  } else {
    if (!quiet) print(">>> Divergences occurred after warmup!")
    if (quiet) return(TRUE)
  }
}



#' Check whether maximum tree depth was hit
#' 
#' @param stan_fit Stan fit object
#' @param quiet Whether or not to print if max treedepth was hit
#' @return If quiet, return TRUE if problematic
#' @export
check_treedepth <- function(stan_fit, quiet = FALSE) {
  sampler_params <- rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)
  max_depth <- attr(stan_fit@sim$samples[[1]], "args")$control$max_treedepth
  treedepths <- do.call(rbind, sampler_params)[, "treedepth__"]
  hit_max <- length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  iter <- length(treedepths)
  hit_max_pct <- hit_max / iter * 100
  
  if (!quiet) {
    print(sprintf("%s of %s iterations saturated the maximum tree depth of %s (%s%%)",
                  hit_max, iter, max_depth, hit_max_pct))
  }
  
  if (hit_max == 0) {
    if (!quiet) print("Max tree depth looks okay")
    if (quiet) return(FALSE)
  } else {
    if (!quiet) print(">>> At least one iteration saturated the maximum tree depth!")
    if (quiet) return(TRUE)
  }
}



#' Checks the energy fraction of missing information (E-FMI)
#' 
#' @param stan_fit Stan fit object
#' @param threshold Threshold below which E-FMI should not fall
#' @param quiet Whether or not to print if E-FMI is below threshold
#' @return If quiet, return TRUE if problematic
#' @export
check_energy <- function(stan_fit, threshold = 0.2, quiet = FALSE) {
  sampler_params <- rstan::get_sampler_params(stan_fit, inc_warmup = FALSE)
  no_warning <- TRUE
  for (chain_i in 1:length(sampler_params)) {
    energies <- sampler_params[chain_i][[1]][, "energy__"]
    num <- sum(diff(energies)**2) / length(energies)
    den <- var(energies)
    if (den == 0) {
      if (!quiet) {
        print(sprintf(">>> Chain %s: energy is constant!", chain_i))
      }
      no_warning <- FALSE
    } else if (num / den < threshold) {
      if (!quiet) {
        print(sprintf(">>> Chain %s: E-FMI = %s", chain_i, num / den))
      }
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print("E-FMI indicated no pathological behavior")
    if (quiet) return(FALSE)
  } else {
    if (!quiet) print(paste0(">>> E-FMI below ", threshold,
                             " suggests need for reparameterization"))
    if (quiet) return(TRUE)
  }
}



#' Checks for low effective sample size per iteration
#' 
#' @param stan_fit Stan fit object
#' @param threshold Threshold below which (n_eff / iter) should not fall
#' @param quiet Whether or not to print if n_eff / iter is too low
#' @return If quiet, return TRUE if problematic
#' @export
check_n_eff <- function(stan_fit, threshold = 0.5, quiet=FALSE) {
  fit_summary <- rstan::summary(stan_fit, probs = c(0.5))$summary
  P <- dim(fit_summary)[[1]] - 1 # (ignore lp__)
  iter <- dim(rstan::extract(stan_fit)[[1]])[[1]]
  
  no_warning <- TRUE
  for (p in 1:P) {
    ratio <- fit_summary[, 5][p] / iter
    if (ratio < threshold) {
      if (!quiet) print(sprintf(">>> n_eff / iter for parameter %s is %s!",
                                rownames(fit_summary)[p], ratio))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print("n_eff / iter looks reasonable for all parameters")
    if (quiet) return(FALSE)
  }
  else {
    if (!quiet) print(paste0(">>> n_eff / iter is less than ", threshold,"!"))
    if (quiet) return(TRUE)
  }
}

#' Checks the potential scale reduction factors
#' 
#' @param stan_fit Stan fit object
#' @param threshold Threshold above which Rhat should not fall
#' @param quiet Whether or not to print if Rhat is too high
#' @return If quiet, return TRUE if problematic
#' @export
check_rhat <- function(stan_fit, threshold = 1.02, quiet = FALSE) {
  fit_summary <- rstan::summary(stan_fit, probs = c(0.5))$summary
  P <- dim(fit_summary)[[1]]
  
  no_warning <- TRUE
  for (p in 1:P) {
    rhat <- fit_summary[, 6][p]
    if ((rhat > threshold) || is.infinite(rhat) || is.nan(rhat)) {
      if (!quiet) print(sprintf(">>> Rhat for parameter %s is %s!",
                                rownames(fit_summary)[p], rhat))
      no_warning <- FALSE
    }
  }
  if (no_warning) {
    if (!quiet) print("Rhat looks reasonable for all parameters")
    if (quiet) return(FALSE)
  } else {
    if (!quiet) print(paste0(">>> At least one Rhat above ", threshold, "!"))
    if (quiet) return(TRUE)
  }
}


#' Check all diagnostics for a stan fit object
#' 
#' @param stan_fit Stan fit object
#' @param quiet Whether or not to print if n_eff / iter is too low
#' @return If quiet, return warning code that can be parsed by 
#' \code{\link{parse_warning_code}}
#' @export
check_all_diagnostics <- function(stan_fit, quiet = FALSE) {
  if (!quiet) {
    check_n_eff(stan_fit)
    check_rhat(stan_fit)
    check_div(stan_fit)
    check_treedepth(stan_fit)
    check_energy(stan_fit)
  } else {
    warning_code <- 0
    # browser()
    if (check_n_eff(stan_fit, quiet = TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 0))
    if (check_rhat(stan_fit, quiet = TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 1))
    if (check_div(stan_fit, quiet = TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 2))
    if (check_treedepth(stan_fit, quiet = TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 3))
    if (check_energy(stan_fit, quiet = TRUE))
      warning_code <- bitwOr(warning_code, bitwShiftL(1, 4))
    
    return(warning_code)
  }
}



#' Translate sampler warning codes
#' 
#' @param warning_code Warning code returned by \code{\link{check_all_diagnostics}}
#' @return None. Warnings are printed.
#' @export
parse_warning_code <- function(warning_code) {
  if (bitwAnd(warning_code, bitwShiftL(1, 0)))
    print("n_eff / iteration warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 1)))
    print("rhat warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 2)))
    print("divergence warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 3)))
    print("treedepth warning")
  if (bitwAnd(warning_code, bitwShiftL(1, 4)))
    print("energy warning")
}



#' Check whether there were convergence problems in simulation replicate
#' 
#' @param replicate Job result from \code{\link{submit_scenario_jobs}}
#' @param n_eff_thresh n_eff / B threshold for \code{\link{check_n_eff}}
#' @param rhat_thresh R hat threshold for \code{\link{check_rhat}}
#' @param detailed Whether to return vector of each probed problem
#' type (\code{TRUE}) or just a summary of whether any problems were
#' detected (\code{FALSE}). Defaults to \code{TRUE}.
#' @return Scalar boolean for whether there were any problems or vector of
#' named problem booleans
#' @export
apply_simstudy_conv_criteria <- function(replicate, 
                                         n_eff_thresh = 0.5, 
                                         rhat_thresh = 1.1,
                                         detailed = TRUE) {
  sf <- replicate$result$stan_fit
  problems <- rep(NA, 4)
  names(problems) <- c("n_eff_bad", "div_bad", "rhat_bad", "energy_bad")
  problems[1] <- check_n_eff(sf, threshold = n_eff_thresh, quiet = TRUE)
  problems[2] <- check_div(sf, quiet = TRUE)
  problems[3] <- check_rhat(sf, threshold = rhat_thresh, quiet = TRUE)
  problems[4] <- check_rhat(sf, quiet = TRUE)
  if (!detailed) {
    problems <- any(problems)
    names(problems) <- "any problem"
  } else {
    problems <- c(problems, any(problems))
    names(problems)[5] <- "any_problem"
  }
  return(problems)
}

