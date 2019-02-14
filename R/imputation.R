#' Convert event observation indicators to censoring types
#' 
#' Vectorized -- types match my likelihood writeup
#' 
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @return Vector of types
#' @export
make_type <- function(dyr, dyt) {
  type <- 1 * (dyr == 0 & dyt == 0) +
          2 * (dyr == 1 & dyt == 0) +
          3 * (dyr == 0 & dyt == 1) +
          4 * (dyr == 1 & dyt == 1)
  return(type)
}



#' Convert linear predictor to Weibull scale
#' 
#' @param lp Length-N vector of linear predictors
#' @param alpha Weibull shape parameter
#' @return Length-N vector of Weibull scales
#' @export
make_scale <- function(lp, alpha) {
  return(exp(-lp / alpha))
}



#' Function to return the lp
#' 
#' @param xmat N x P design matrix with no intercept
#' @param beta Length-P vector of regression coefficients
#' @param kappa Scalar baseline hazard
#' @return Length-N vector of linear predictors to make scale
#' @export
make_ph_lp <- function(xmat, beta, kappa) {
  if (is.vector(xmat)) {
    xmat <- as.matrix(xmat, ncol = 1)
  }
  return(xmat %*% beta + log(kappa))
}



#' Convert parameters and covariates to a Weibull scale
#' 
#' @param xmat N x P design matrix
#' @param beta P x 3 matrix of regression coefficients
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param alpha Length-3 vector of Weibull shapes
#' @param frailty Scalar or vector of frailties. Defaults to 1 (reference value).
#' @return N x 3 matrix of Weibull scales
#' @export
make_full_scale <- function(xmat, beta, kappa, alpha, frailty = 1) {
  scale <- matrix(NA, nrow = NROW(xmat), ncol = 3)
  for (g in 1:3) {
    scale[, g] <- make_scale(make_ph_lp(xmat = matrix(xmat, ncol = NROW(beta)),
                                        beta = beta[, g], 
                                        kappa = kappa[g] * frailty), 
                             alpha = alpha[g])
  }
  return(scale)
}



#' Impute a vector of posterior frailties
#' 
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @param xmat N x P design matrix
#' @param alpha Length-3 vector of Weibull shapes
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param beta P x 3 matrix of regression coefficients
#' @param sigma Variance of frailties
#' @return Length-N vector of frailties
#' @export
impute_frailty <- function(yr, yt, dyr, dyt, xmat, alpha, kappa, beta, sigma) {
  n <- length(dyr)
  soj <- max(0, yt - yr)
  scale <- make_full_scale(xmat = xmat, 
                           beta = beta, 
                           kappa = kappa, 
                           alpha = alpha, 
                           frailty = 1)
  a1 <- 1 / sigma + dyr + dyt
  a2 <- 1 / sigma + 
        -pweibull(yr, alpha[1], scale[, 1], lower.tail = FALSE, log.p = TRUE) + 
        -pweibull(yr, alpha[2], scale[, 2], lower.tail = FALSE, log.p = TRUE) + 
        -pweibull(soj, alpha[3], scale[, 3], lower.tail = FALSE, log.p = TRUE)
  if (anyNA(c(a1, a2))) {
    browser()
  }
  return(rgamma(n, shape = a1, rate = a2))
}



#' Generate a Weibull truncated to have no mass below a lower bound
#' 
#' Vectorized!
#' 
#' @param shape Scalar or length-N vector of Weibull shape parameters
#' @param scale Scalar or length-N vector of Weibull scale parameters
#' @param lb Scalar or length-N vector of lower bounds
#' @return Length-N vector of truncated Weibulls
#' @export
rtweibull <- function(shape, scale, lb) {
  p <- pweibull(lb, shape = shape, scale = scale, lower.tail = TRUE)
  u <- runif(n = length(p), min = p, max = 1)
  if (anyNA(u)) browser()
  rv <- qweibull(p = u, shape = shape, scale = scale)
  return(rv)
}



#' Impute (if necessary) outcomes for observed arm of a single censoring type 
#' 
#' Only works for one type at a time
#' 
#' @param type Censoring type 1-4
#' @param frailty Imputed frailty
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @param xmat N x P design matrix
#' @param alpha Length-3 vector of Weibull shapes
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param beta P x 3 matrix of regression coefficients
#' @return N x 4 matrix of uncensored (yr, yt, dyr, dyt)
#' @export
impute_obs_by_type <- function(type, frailty, yr, yt, dyr, dyt, xmat, 
                               alpha, kappa, beta) {
  stopifnot(length(type) == 1)
  scale <- make_full_scale(xmat = xmat, beta = beta, kappa = kappa, alpha = alpha, 
                           frailty = frailty)
  if (type == 1) {
    R_cand     <- rtweibull(shape = alpha[1], scale[, 1], lb = yr)
    D_cand     <- rtweibull(shape = alpha[2], scale[, 2], lb = yr)
    soj_cand   <- rtweibull(shape = alpha[3], scale[, 3], lb = 0)
    dyr_uncens <- ifelse(R_cand < D_cand, 1, 0)
    yr_uncens  <- ifelse(R_cand < D_cand, R_cand, D_cand)
    yt_uncens  <- ifelse(R_cand < D_cand, R_cand + soj_cand, D_cand)
  } else if (type == 2) {
    soj_uncens <- rtweibull(shape = alpha[3], scale[, 3], lb = yt - yr)
    dyr_uncens <- dyr
    yr_uncens  <- yr
    yt_uncens  <- yr + soj_uncens
  } else if (type %in% c(3, 4)) {
    yr_uncens <- yr
    yt_uncens <- yt
    dyr_uncens <- dyr
  }
  dyt_uncens <- 1
  return(cbind(yr_uncens, yt_uncens, dyr_uncens, dyt_uncens))
}



#' Impute (if necessary) uncensored outcomes in observed arm
#' 
#' @param frailty Imputed frailty
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @param xmat N x P design matrix
#' @param alpha Length-3 vector of Weibull shapes
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param beta P x 3 matrix of regression coefficients
#' @return N x 4 matrix of observed z and imputed 
#' (frailty, yr0, yt0, dyr0, dyt0, yr1, yt1, dyr1, dyt1)
#' @export
impute_obs <- function(frailty, yr, yt, dyr, dyt, xmat, 
                       alpha, kappa, beta) {
  dat <- matrix(nrow = NROW(xmat), ncol = 4)
  colnames(dat) <- c("yr_uncens", "yt_uncens", "dyr_uncens", "dyt_uncens")
  
  type <- make_type(dyr = dyr, dyt = dyt)
  
  for (ti in 1:4) {
    if (any(type == ti)) {
      is_ti <- which(type == ti)
      dat[is_ti, ] <- impute_obs_by_type(type = ti, 
                                         frailty = frailty[is_ti], 
                                         yr = yr[is_ti], yt = yt[is_ti], 
                                         dyr = dyr[is_ti], dyt = dyt[is_ti], 
                                         xmat = xmat[is_ti, ], 
                                         alpha = alpha, 
                                         kappa = kappa, 
                                         beta = beta)
    }
  }
  return(dat)
}




#' Impute causally missing outcomes from a single parameter set
#' 
#' @param frailty Imputed frailty
#' @param xmat N x P design matrix
#' @param alpha Length-3 vector of Weibull shapes
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param beta P x 3 matrix of regression coefficients
#' @return N x 4 matrix of uncensored (yr, yt, dyr, dyt) for counter-to-fact arm
#' @export
impute_mis <- function(frailty, xmat, alpha, kappa, beta) {
  dat <- impute_obs_by_type(type = 1, frailty = frailty, 
                            yr = 0, yt = 0, dyr = 0, dyt = 0, 
                            xmat = xmat, 
                            alpha = alpha, kappa = kappa, beta = beta)
  colnames(dat) <- c("yr_mis", "yt_mis", "dyr_mis", "dyt_mis")
  return(dat)
}



#' Do outcome posterior prediction for a single parameter draw
#' 
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @param z Assigned treatment vector
#' @param xmat N x P design matrix
#' @param alpha Length-6 vector of Weibull shapes
#' @param kappa Length-6 vector of Weibull baseline hazards
#' @param beta P x 3 or P x 6 matrix of regression coefficients
#' @param sigma Scalar frailty variance
#' @param frailty_type Whether to impute that all frailties are exactly 
#' one ("reference"), impute from data ("impute"), use provided frailties 
#' ("given"), or take a quantile of the distribution implied by sigma 
#' ("quantile"). Reference is useful for ranking covariate patterns by risk. 
#' Defaults to "impute".
#' @param frailty Provided length-N frailty vector if frailty_type = "given"
#' @param frailty_q Quantile (from 0 to 1) if frailty_type = "quantile"
#' @param data_type Whether to impute for sample ("sample") or for new people
#' ("new")
#' @return N-row data frame of uncensored (z, frailty, yr, yt, dyr, dyt)
#' @export 
posterior_predict_draw <- function(yr, yt, dyr, dyt, z, xmat,
                                   alpha, kappa, beta, sigma,
                                   frailty_type = "impute", 
                                   frailty = NULL, 
                                   frailty_q = NULL,
                                   data_type = "sample") {
  if (!(data_type %in% c("sample", "new"))) {
    stop("Need to specify data_type as 'sample' or 'new'")
  } else if (data_type == "sample" && frailty_type == "quantile") {
    stop("Quantile frailty specification only makes sense for new observations")
  }
  if (frailty_type == "given" && is.null(frailty)) {
    stop("Need to provide vector of frailties if frailty_type == given!")
  } else if (frailty_type == "quantile" && is.null(frailty_q)) {
    stop("Need to provide quantile for frailty if frailty_type == quantile!")
  }
  n  <- length(z)
  if (frailty_type != "given") {
    frailty <- rep(NA, n)
  }
  shared_beta <- (ncol(beta) == 3) * 1
  xmat <- as.matrix(xmat)
  obs <- mis <- matrix(NA, nrow = n, ncol = 4)
  for (zval in 0:1) {
    z_which <- which(z == zval)
    nz <- length(z_which)
    if (zval == 0) {
      obs_beta_indices <- obs_indices <- 1:3
      mis_indices <- 4:6
      if (shared_beta == 1) {
        mis_beta_indices <- 1:3
      } else {
        mis_beta_indices <- 4:6
      }
    } else if (zval == 1) {
      obs_indices <- 4:6
      mis_beta_indices <- mis_indices <- 1:3
      if (shared_beta == 1) {
        obs_beta_indices <- 1:3
      } else {
        obs_beta_indices <- 4:6
      }
    }
    if (frailty_type == "impute") {
      frailty[z_which] <- impute_frailty(yr = yr[z_which], yt = yt[z_which], 
                                         dyr = dyr[z_which], dyt = dyt[z_which], 
                                         xmat = as.matrix(xmat[z_which, ], 
                                                          nrow = nz), 
                                         alpha = alpha[obs_indices], 
                                         kappa = kappa[obs_indices], 
                                         beta = beta[ , obs_beta_indices], 
                                         sigma = sigma)
    } else if (frailty_type == "reference") {
      frailty[z_which] <- rep(1, length(z_which))  
    } else if (frailty_type == "quantile") {
      frailty[z_which] <- rep(qgamma(p = frailty_q, 1 / sigma, 1 / sigma), 
                              length(z_which))
    } 
    
    if (data_type == "sample") {
      obs[z_which, ]   <- impute_obs(frailty = frailty[z_which], 
                                     yr = yr[z_which], yt = yt[z_which], 
                                     dyr = dyr[z_which], dyt = dyt[z_which], 
                                     xmat = as.matrix(xmat[z_which, ], 
                                                      nrow = nz), 
                                     alpha = alpha[obs_indices], 
                                     kappa = kappa[obs_indices], 
                                     beta = beta[ , obs_beta_indices])
    } else if (data_type == "new") {
      # New people = simulate from "observed" treatment arm process but do not 
      # truncate to agree with existing outcomes
      obs[z_which, ]   <- impute_mis(frailty = frailty[z_which], 
                                     xmat = matrix(xmat[z_which, ], nrow = nz), 
                                     alpha = alpha[obs_indices],
                                     kappa = kappa[obs_indices], 
                                     beta = beta[ , obs_beta_indices])
    }
    mis[z_which, ]   <- impute_mis(frailty = frailty[z_which], 
                                   xmat = matrix(xmat[z_which, ], nrow = nz), 
                                   alpha = alpha[mis_indices],
                                   kappa = kappa[mis_indices], 
                                   beta = beta[ , mis_beta_indices])
  }
  out0 <- out1 <- obs * NA
  out0[z == 0, ] <- obs[z == 0, ]
  out1[z == 0, ] <- mis[z == 0, ]
  out0[z == 1, ] <- mis[z == 1, ]
  out1[z == 1, ] <- obs[z == 1, ]
  
  colnames(out0) <- c("yr0_imp", "yt0_imp", "dyr0_imp", "dyt0_imp")
  colnames(out1) <- c("yr1_imp", "yt1_imp", "dyr1_imp", "dyt1_imp")
  
  return(cbind(z, frailty, out0, out1))
}



#' Do outcome posterior prediction for all parameter draws from Stan fit
#' 
#' @param stan_fit Stan fit object
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @param z Assigned treatment vector
#' @param xmat N x P design matrix
#' @param frailty_type Whether to impute that all frailties are exactly 
#' one ("reference"), impute from data ("impute"), use provided frailties 
#' ("given"). Reference is useful for ranking covariate patterns by risk. 
#' Defaults to "impute".
#' @param frailty Provided length-N frailty vector or N x R matrix if 
#' frailty_type = "given"
#' @return N-row data frame of uncensored (z, yr, yt, dyr, dyt)
#' @export 
posterior_predict_sample <- function(stan_fit, yr, yt, dyr, dyt, z, xmat,
                                     frailty_type = "impute", frailty = NULL) {
  
  aalpha <- t(as.array(extract(stan_fit, par = "alpha")[["alpha"]]))
  akappa <- t(as.array(extract(stan_fit, par = "kappa")[["kappa"]]))
  if (frailty_type != "reference") {
    asigma <- as.array(extract(stan_fit, par = "sigma")[["sigma"]])  
  }
  abeta <- aperm(as.array(extract(stan_fit, par = "beta")[["beta"]]), 
                 c(2, 3, 1))
  afrailty <- frailty
  
  R <- NCOL(aalpha)
  n <- length(z)
  pps <- array(NA, dim = c(n, 10, R))
  for (r in 1:R) {
    alpha <- aalpha[, r]
    kappa <- akappa[, r]
    if (frailty_type == "reference") {
      sigma = 0
    } else {
      sigma <- asigma[r]
    } 
    
    if (!is.null(afrailty) && NROW(afrailty) == n && NCOL(afrailty) == R) {
      frailty <- afrailty[, r]
    }
    
    beta <- matrix(abeta[ , , r], nrow = NCOL(xmat))
    a <- posterior_predict_draw(yr = yr, yt = yt, dyr = dyr, dyt = dyt, 
                                z = z, xmat = xmat,
                                alpha = alpha, kappa = kappa, 
                                beta = beta, sigma = sigma,
                                frailty_type = frailty_type,
                                frailty = frailty)
    pps[ , , r] <- a
  }
  dimnames(pps)[[2]] <- colnames(a)
  
  return(pps)
}



#' Calculate principal states at a scalar t
#' 
#' @param eval_t Time at which to evaluate survival
#' @param pp Posterior predictive data set containing yt0_imp and yt1_imp.
#' (this is one slice of the array from output of 
#' \code{\link{posterior_predict_sample}}, i.e., \code{res[, , 1]})
#' @return Character vector of principal states
#' @export
make_pstates <- function(eval_t, pp) {
  pp <- as.data.frame(pp)
  pstate <- rep(NA, NROW(pp))
  pstate[(pp$yt0_imp > eval_t) & (pp$yt1_imp > eval_t)] <- "AA"
  pstate[(pp$yt0_imp > eval_t) & (pp$yt1_imp <= eval_t)] <- "TK"
  pstate[(pp$yt0_imp <= eval_t) & (pp$yt1_imp > eval_t)] <- "TS"
  pstate[(pp$yt0_imp <= eval_t) & (pp$yt1_imp <= eval_t)] <- "DD"
  return(pstate)
}



#' Calculate SACE at a certain time point
#' 
#' @param eval_t Time at which to evaluate principal state and cumulative 
#' incidence of the non-terminal event
#' @param pp Posterior predictive data set containing yt0_imp and yt1_imp.
#' (this is one slice of the array from output of 
#' \code{\link{posterior_predict_sample}}, i.e., \code{res[, , 1]})
#' @return Scalar estimate of SACE at that eval_t
#' @export
calculate_tv_sace <- function(eval_t, pp) {
  pp <- as.data.frame(pp)
  pstate <- make_pstates(eval_t, pp)
  r1_by_t <- r0_by_t <- rep(0, NROW(pp))
  r0_by_t[pp$yr0 < eval_t] <- 1
  r1_by_t[pp$yr1 < eval_t] <- 1
  diff_by_t <- r1_by_t - r0_by_t
  tv_sace <- mean(diff_by_t[pstate == "AA"])
  return(tv_sace)
}



#' Calculation of TV-SACE for multiple time points from one set of P.O.
#' 
#' Vectorized version of \code{\link{calculate_tv_sace}}
#' @param eval_t Length-T vector of times at which to evaluate principal state 
#' and cumulative incidence of the non-terminal event
#' @param pp Posterior predictive data set containing yt0_imp and yt1_imp
#' @return Length-T vector of TV-SACE for each eval_t
#' @export
v_calculate_tv_sace <- Vectorize(FUN = calculate_tv_sace, 
                                 vectorize.args = "eval_t")



#' Calculated the restricted mean SACE given a set of posterior predictive draws
#' 
#' @param eval_t Time at which to evaluate principal state and accumulated 
#' benefit with respect to the non-terminal event
#' @param pp Posterior predictive data set containing yt0_imp and yt1_imp.
#' (this is one slice of the array from output of 
#' \code{\link{posterior_predict_sample}}, i.e., \code{res[, , 1]})
#' @return Scalar estimate of RM-SACE at that eval_t
#' @export
calculate_rm_sace <- function(eval_t, pp) {
  pp <- as.data.frame(pp)
  pstate <- make_pstates(eval_t, pp)
  r1_by_t <- r0_by_t <- rep(0, NROW(pp))
  r0_by_t <- pmin(eval_t, pp$yr0)
  r1_by_t <- pmin(eval_t, pp$yr1)
  diff_by_t <- r1_by_t - r0_by_t
  rm_sace <- mean(diff_by_t[pstate == "AA"])
  return(rm_sace)
}



#' Calculation of RM-SACE for multiple time points from one set of P.O.
#' 
#' Vectorized version of \code{\link{calculate_rm_sace}}
#' @param eval_t Length-T vector of times at which to evaluate principal state 
#' and cumulative incidence of the non-terminal event
#' @param pp Posterior predictive data set containing yt0_imp and yt1_imp
#' @return Length-T vector of RM-SACE for each eval_t
#' @export
v_calculate_rm_sace <- Vectorize(FUN = calculate_rm_sace, 
                                 vectorize.args = "eval_t")



#' Give reasonable middle-ish times for evaluation of causal effects
#' 
#' First time is median non-terminal event time in control (if no semicompeting risk)
#' Second time is double that.
#' 
#' @return Vector of two times for evaluating causal effects
#' @export
get_eval_t <- function() {
  params <- return_dgp_parameters(scenario = "1")
  kappa1 <- params$control$kappa1
  alpha1 <- params$control$alpha1
  t1 <- exp(-log(kappa1)/alpha1) * log(2)^(1 / alpha1)
  return(c(t1 = t1, t2 = 2 * t1))
}



#' Make covariate risk scores for the terminal event 
#' 
#' Create risk scores for each covariate pattern by calculating the the posterior
#' probability of the individual being always-alive at eval_t, ignoring frailties. 
#' (The exact values will be different for different eval_t, but the ordering will 
#' not.) Individuals with the same covariate pattern will have the same risk score.
#' 
#' @param stan_fit Stan fit object
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @param z Assigned treatment vector
#' @param xmat N x P design matrix
#' @param eval_t Time at which to evaluate P(always-alive)
#' @return Length-N vector of covariate risk scores
#' @export
make_terminal_risk_scores <- function(stan_fit, yr, yt, dyr, dyt, z, xmat, 
                                      eval_t) {
  pp_ref <- posterior_predict_sample(stan_fit, yr, yt, dyr, dyt, z, xmat,
                                     frailty_type = "reference")
  is_aa_t1_ref <- (apply(X = pp_ref, MARGIN = 3, FUN = make_pstates, 
                         eval_t = eval_t) == "AA")
  risk_score <- rowMeans(is_aa_t1_ref)
  return(risk_score)
}



#' Function to do posterior prediction for a new design matrix based on
#' parameter draws from an existing Stan fit object
#' 
#' @param xnew Design matrix with J rows (observations) and P columns 
#' (covariates)
#' @param stan_fit Stan fit object like that returned as an element of the list
#' output from \code{\link{run_scr_replicate}}
#' @param frailty Frailty, if provided. Useful for passing in a 1 to get reference.
#' @param frailty_q Frailty quantile, if desiring to draw from frailty distribution
#' implied by posterior for sigma
#' @return Array of posterior predictions
#' @export
posterior_predict_xnew <- function(xnew, stan_fit, frailty = NULL, frailty_q = 0.5) {
  aalpha <- t(as.array(extract(stan_fit, par = "alpha")[["alpha"]]))
  akappa <- t(as.array(extract(stan_fit, par = "kappa")[["kappa"]]))
  asigma <- as.array(extract(stan_fit, par = "sigma")[["sigma"]])  
  abeta  <- aperm(as.array(extract(stan_fit, par = "beta")[["beta"]]), 
                  c(2, 3, 1))
  R <- NCOL(aalpha)
  nnew <- NROW(xnew)
  znew <- round(seq(0, 1, length.out = nnew))
  zeros <- rep(0, nnew)
  pps <- array(NA, dim = c(nnew, 10, R))
  for (r in 1:R) {
    alpha <- aalpha[, r]
    kappa <- akappa[, r]
    sigma <- asigma[r]
    beta <- matrix(abeta[ , , r], nrow = NCOL(xnew))
    if (!(is.null(frailty))) {
      a <- posterior_predict_draw(yr = zeros, 
                                  yt = zeros, 
                                  dyr = zeros, 
                                  dyt = zeros, 
                                  z = znew, 
                                  xmat = xnew,
                                  alpha, kappa, beta, sigma,
                                  frailty_type = "given", 
                                  frailty = frailty,
                                  data_type = "new")
    } else {
      a <- posterior_predict_draw(yr = zeros, 
                                  yt = zeros, 
                                  dyr = zeros, 
                                  dyt = zeros, 
                                  z = znew, 
                                  xmat = xnew,
                                  alpha, kappa, beta, sigma,
                                  frailty_type = "quantile", 
                                  frailty_q = frailty_q,
                                  data_type = "new")
    }
    pps[ , , r] <- a
  }
  dimnames(pps)[[2]] <- colnames(a)
  
  return(pps)
}



