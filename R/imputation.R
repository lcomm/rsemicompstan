#' Convert event observation indicators to censoring types
#' 
#' Vectorized -- matches my likelihood writeup
#' 
#' @param dyr Whether nonterminal event was observed
#' @param dyt Whether terminal event was observed
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
    scale[, g] <- make_scale(make_ph_lp(xmat = xmat, 
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
#' @param dyr Non-terminal event observation indicator
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
        -pweibull(yr, alpha[1], scale[, 1], lower = FALSE, log = TRUE) + 
        -pweibull(yr, alpha[2], scale[, 2], lower = FALSE, log = TRUE) + 
        -pweibull(soj, alpha[3], scale[, 3], lower = FALSE, log = TRUE)
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
  return(qnorm(p = p, shape = shape, scale = scale))
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
#' @param dyr Non-terminal event observation indicator
#' @param xmat N x P design matrix
#' @param alpha Length-3 vector of Weibull shapes
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param beta P x 3 matrix of regression coefficients
#' @return N x 4 matrix of uncensored (yr, yt, dyr, dyt)
#' @export
impute_obs_by_type <- function(type, frailty, yr, yt, dyr, dyt, xmat, 
                               alpha, kappa, beta) {
  stopifnot(length(type) == 1)
  scale <- make_full_scale(xmat = x, beta = beta, kappa = kappa, alpha = alpha, 
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



#' Impute (if necessary) outcomes for observed arm
#' 
#' @param type Length-N vector of types
#' @param frailty Imputed frailty
#' @param yr Last observed non-terminal time
#' @param yt Last observed terminal time
#' @param dyr Non-terminal event observation indicator
#' @param dyr Non-terminal event observation indicator
#' @param xmat N x P design matrix
#' @param alpha Length-3 vector of Weibull shapes
#' @param kappa Length-3 vector of Weibull baseline hazards
#' @param beta P x 3 matrix of regression coefficients
#' @return N x 4 matrix of uncensored (yr, yt, dyr, dyt)
#' @export
impute_obs_by_type <- function(type, frailty, yr, yt, dyr, dyt, xmat, 
                               alpha, kappa, beta) {
  
  dat <- matrix(nrow = NROW(xmat), col = 4)
  colnames(dat) <- c("yr_uncens", "yt_uncens", "dyr_uncens", "dyt_uncens")
  
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



#TODO(LCOMM): sample splitting into treated and control
#TODO(LCOMM): full imputation within an arm
#TODO(LCOMM): synthesize output with previous code for SACE and RM-SACE
#TODO(LCOMM): unit tests

# summary(resg)
# 
# aalpha <- t(as.array(extract(resg, par = "alpha")[["alpha"]]))
# akappa <- t(as.array(extract(resg, par = "kappa")[["kappa"]]))
# asigma <- as.array(extract(resg, par = "sigma")[["sigma"]])
# abeta <- aperm(as.array(extract(resg, par = "beta")[["beta"]]), c(2, 3, 1))
# alpha <- aalpha[, 1]
# kappa <- akappa[, 1]
# sigma <- asigma[1]
# beta <- matrix(abeta[ , , 1], ncol = 3)
# 
# x = x1
# z = z
# yr = dat$y1
# yt = dat$y2
# dyr = dat$delta1
# dyt = dat$delta2