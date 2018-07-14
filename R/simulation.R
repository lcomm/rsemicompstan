#' Assign a binary treatment Z 
#' 
#' @param n Number of observations to randomize
#' @param treated_fraction Fraction of sample to have = 1
#' @return Length-n vector of 0s (if control) or 1s (if treated)
#' @export
assign_Z <- function(n, treated_fraction = 0.5) {
  return((1:n %in% sample(1:n, round(n * treated_fraction))) * 1)
}



#' Function to simulate basic covariate data
#' 
#' @param n Number of observations
#' @param p Number of covariates to simulate
#' @return Data frame of 3 observations and 3 columns
#' @export
simulate_covariates <- function(n = 100, p = 3) {
  
  # Very basic parameter check
  stopifnot(p > 0)
  
  # Simulate binary covariates of moderate prevalences
  # Covariates centered around true prevalence
  dat <- data.frame(row.names = 1:n)
  if (p > 1) {
    prevs <- seq(0.4, 0.2, length.out = p)
    for (j in 1:p) {
      dat[[paste0("X", j)]] <- rbinom(n, size = 1, prob = prevs[j]) - prevs[j]
    }
  } else if (p == 1) {
    dat$X1 <- rbinom(n, size = 1, prob = 0.5) - 0.5
  }
  
  # Return
  return(dat)

}



#' Function to generate frailties
#' 
#' Draw n i.i.d. frailties from a Gamma($\sigma^{-1}$, $\sigma^{-1}$)\
#'
#' @param n Number of observations
#' @param sigma Variance for gamma frailty
#' @return Length-n numeric matrix of frailties
#' @export
simulate_frailty <- function(n, sigma = 1) {
  
  stopifnot(is.numeric(sigma), length(sigma) == 1, sigma > 0)
  
  # Variance for rgamma parameterization is shape * scale^2
  return(rgamma(n, shape = 1 / sigma, scale = sigma))
  
}


#' Fit a maximum likelihood model to center log(alpha) and log(kappa)
#' at reasonable values
#' 
#' @param time Time of event or censoring
#' @param event Event indicator (if not provided, time is assumed to have no
#' censoring)
#' @return Named vector of prior means for log(alpha) and log(kappa)
#' @export
get_prior_mean_from_mle <- function(time, event = NULL) {
  
  require("survival")
  
  if (is.null(event)) {
    fit <- survreg(Surv(time = time) ~ 1, dist = "weibull")
  } else {
    stopifnot(is.numeric(event), all(event %in% c(0, 1)))
    fit <- survreg(Surv(time = time, event = event) ~ 1, dist = "weibull")
  }
  alpha <- 1 / fit$scale
  log_alpha <- log(alpha)
  rwei_scale <- exp(unname(coef(fit)[1]))
  kappa <- rwei_scale^(-alpha)
  log_kappa <- log(kappa)
  
  return(c(log_alpha = log_alpha, log_kappa = log_kappa))
  
}



#' Make a full frailty or regression coefficient matrix from vector or smaller 
#' matrix
#' 
#' @param x Length-N vector of common values to be used for all 6 models
#' or an N x 3 matrix of transition-specific values
#' @return N x 6 matrix
#' @export
expand_to_full <- function(x) {
  
  stopifnot(is.numeric(x), anyNA(x) == FALSE)
  n <- NROW(x)
  
  if (is.matrix(x)) {
    stopifnot(NCOL(x) %in% c(3, 6))
    x <- matrix(x, nrow = n, ncol = 6)
  } else if (is.null(dim(x))){
    x <- matrix(x, nrow = n, ncol = 6)
  } else {
    stop("Cannot parse expansion to 6 model structure")
  }
  
  return(x)
  
}



#' Simulate Weibull potential outcomes from scratch, making a science table
#' 
#' @param xmat N x P design matrix
#' @param beta P x 3 or P x 6 matrix for regression coefficient i in model j
#' @param alpha Length-6 vector of Weibull shape parameters
#' @param kappa Length-6 vector of baseline hazards for PH model
#' @param frailty Length-N vector, N x 3, or N x 6 matrix of frailties
#' @return N x 8 dataframe with potential outcomes under control and treated
#' @export
simulate_science <- function(xmat, beta, alpha, kappa, frailty) {
  
  # Parameter checks
  stopifnot(all(alpha > 0), all(kappa > 0), 
            NCOL(xmat) == NROW(beta), NROW(xmat) == NROW(frailty))
  stopifnot(length(kappa) == 6, length(alpha) == 6)
  
  # Initialize data frame for potential outcomes
  n <- NROW(xmat)
  dat <- data.frame(row.names = 1:n)
  
  # Expand if not already in 6-column format
  beta <- expand_to_full(beta)
  frailty <- matrix(frailty, nrow = n, ncol = 6)
  kappa <- matrix(rep(kappa, each = n), ncol = 6)
  alpha <- matrix(rep(alpha, each = n), ncol = 6)

  # Pre-calculate linear predictor and Weibull scale parameters
  lp <- xmat %*% beta + log(frailty) + log(kappa)
  scales <- exp(-lp/alpha)
  
  # Pre-calculate shared censoring time
  # Censor uniformly between 50th and 90th percentiles of process 2
  # under the control condition
  lb <- qweibull(0.5, shape = alpha[1,2], scale = mean(scales[1, 2]))
  ub <- qweibull(0.9, shape = alpha[1,2], scale = mean(scales[1, 2]))
  C <- runif(n, min = lb, max = ub)
  
  for (start_i in c(1, 4)) {
    z <- ifelse(start_i == 1, 0, 1)
    
    # Generate data
    R <- rweibull(n, shape = alpha[, start_i], scale = scales[, start_i])
    D <- rweibull(n, shape = alpha[, start_i + 1], scale = scales[, start_i + 1])
    soj <- rweibull(n, shape = alpha[, start_i + 2], scale = scales[, start_i + 2])
    yesR <- (R < D)
    if (!any(yesR)) {
      stop("simulate_poutcomes: did not simulate any non-terminal events!")
    } else if (all(yesR)) {
      stop("simulate_poutcomes: simulated all terminal events of same type!")
    }
    D[yesR] <- R[yesR] + soj[yesR]
    deltaR <- rep(NA, n)
    deltaD <- rep(NA, n)
    yr <- R
    yt <- D
    ind01 <- which((D < R) & (D < C))
    yr[ind01] <- D[ind01]
    deltaR[ind01] <- 0
    deltaD[ind01] <- 1
    ind10 <- which((R < D) & (R < C) & (D >= C))
    yt[ind10] <- C[ind10]
    deltaR[ind10] <- 1
    deltaD[ind10] <- 0
    ind00 <- which((R >= C) & (D >= C))
    yr[ind00] <- C[ind00]
    yt[ind00] <- C[ind00]
    deltaR[ind00] <- 0
    deltaD[ind00] <- 0
    ind11 <- which((R < C) & (D < C) & (R < D))
    deltaR[ind11] <- 1
    deltaD[ind11] <- 1
    dyr <- deltaR
    dyt <- deltaD
  
    # Save in data set
    to_name <- c("R", "D", "deltaR", "deltaD", "yr", "yt", "dyr", "dyt")
    arm_names <- paste0(to_name, z)
    dat[[arm_names[1]]] <- R
    dat[[arm_names[2]]] <- D
    dat[[arm_names[3]]] <- deltaR
    dat[[arm_names[4]]] <- deltaD
    dat[[arm_names[5]]] <- yr
    dat[[arm_names[6]]] <- yt
    dat[[arm_names[7]]] <- dyr
    dat[[arm_names[8]]] <- dyt
    
  }
  
  # Return
  return(dat)
  
}



#' Convert a full science table to observed data and counterfactuals
#' 
#' @param dat Data frame containing assigned treatment Z, 
#' event times R0, R1, D0, D1, and observation indicators deltaR0, deltaR1,
#' deltaD0, and deltaD1
#' @return Data frame with observed Y and counter-to-fact outcomes (_mis)
#' @export
convert_science_to_obs <- function(dat) {
  
  dat$yr  <- ifelse(dat$Z, dat$yr1, dat$yr0)
  dat$yt  <- ifelse(dat$Z, dat$yt1, dat$yt0)
  dat$dyr <- ifelse(dat$Z, dat$dyr1, dat$dyr0)
  dat$dyt <- ifelse(dat$Z, dat$dyt1, dat$dyt0)
  
  dat$yr_mis  <- ifelse(dat$Z == 0, dat$yr1, dat$yr0)
  dat$yt_mis  <- ifelse(dat$Z == 0, dat$yt1, dat$yt0)
  dat$dyr_mis <- ifelse(dat$Z == 0, dat$dyr1, dat$dyr0)
  dat$dyt_mis <- ifelse(dat$Z == 0, dat$dyt1, dat$dyt0)
  
  return(dat)
}


#' Simulate data
#' 
#' @param n Number of observations
#' @param alpha Weibull shape parameters
#' @param beta Regression coefficient matrix
#' @param kappa Baseline hazard vector
#' @param sigma Variance for frailties
#' @param p Number of covariates
#' @return List containg: dat, a data frame of N observations
#' @examples
#' \dontrun{
#' set.seed(42)
#' dat <- simulate_data(n = 5000, alpha = c(1, 1.2, 1, 1, 1.2, 1),
#'                      beta = matrix((1:3)/10, nrow = 3, ncol = 6),
#'                      kappa = 1 + (1:6)/10, sigma = 0.0001, p = 3)
#' bfit <- scr_no_frailty_stan(x = cbind(dat$X1, dat$X2, dat$X3), 
#'                             z = as.numeric(dat$Z), yr = dat$yr, yt = dat$yt,
#'                             dyr = dat$dyr, dyt = dat$dyt)
#' }
#' @export
simulate_data <- function(n, alpha, beta, kappa, sigma, p = 3) {
  
  stopifnot(NROW(beta) == p, NCOL(beta) %in% c(1, 3, 6))
  
  cov_dat <- simulate_covariates(n = n, p = p)
  frailty <- simulate_frailty(n = n, sigma = sigma)
  xmat <- model.matrix(~ -1 + ., data = cov_dat)
  out_dat <- simulate_science(xmat = xmat, beta = beta, alpha = alpha, 
                              kappa = kappa, frailty = frailty)
  out_dat <- cbind(cov_dat, out_dat)
  out_dat$Z <- assign_Z(n, treated_fraction = 0.5)
  out_dat <- convert_science_to_obs(out_dat)
  return(out_dat)
  
}



#' Make a design matrix of any covariate with X + a number
#' 
#' Will not include an intercept
#' 
#' @param dat Data frame of N observations and P X covariates
#' @return N x P model matrix
#' @export
make_xmat_all_X <- function(dat){
  is_X <- grep("X[:digit:]*", x = colnames(dat), fixed = FALSE)
  return(as.matrix(model.matrix(~ -1 + ., data = dat[, is_X])))
}


#' Make list containing prior means for alpha and kappa
#' 
#' Fits frequentist intercept-only models to center log-alpha and log-kappa
#' around reasonable values
#' 
#' @param yr Nonterminal event time
#' @param yt Terminal event time
#' @param dyr Nonterminal event observation indicator
#' @param dyt Terminal event observation indicator
#' @return Named list of length-6 vectors of prior means for log(alpha) and 
#' log(kappa)
make_prior_means <- function(yr, yt, dyr, dyt) {
  log_kappa_pmean <- log_alpha_pmean <- rep(NA, 6)
  soj <- yt - yr
  time_h1  <- yr
  event_h1 <- dyr
  time_h2  <- yt[dyr == 0]
  event_h2 <- dyt[dyr == 0]
  time_h3  <- soj[dyr == 1]
  event_h3 <- dyt[dyr == 1]
  pmeans   <- rbind(get_prior_mean_from_mle(time_h1, event_h1),
                    get_prior_mean_from_mle(time_h2, event_h2),
                    get_prior_mean_from_mle(time_h3, event_h3))
  return(list(log_alpha_pmean = rep(pmeans[, 1], 2),
              log_kappa_pmean = rep(pmeans[, 2], 2)))
}



