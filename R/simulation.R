# Functions for simulating illness-death data


#' Function to simulate basic covariate data
#' 
#' @param n Number of observations
#' @param p Number of covariates to simulate
#' @return Data frame of 3 observations and 3 columns
simulate_covariate_data <- function(n = 100, p = 3) {
  
  # Very basic parameter check
  stopifnot(p > 1)
  
  # Simulate binary covariates of moderate prevalences
  dat <- data.frame(row.names = 1:n)
  if (p > 1) {
    prevs <- seq(0.4, 0.2, length.out = p)
    for (j in 1:p) {
      dat[[paste0("X", j)]] <- rbinom(n, size = 1, prob = prevs[j])
    }
  } else if (p == 1) {
    dat$X1 <- rbinom(n, size = 1, prob = 0.5)
  }
  
  # Return
  return(dat)
  
}

  
#' Make design matrices from data (or simulate all data)
#' 
#' @param n Number of observations
#' @param formulas List of 6 regression formulas; if none, all covariate main 
#' effects will be included. If 1, same formula will be used for all 6
#' @param covariate_data Data frame containing covariates; if none, then will
#' simulate data frame with 3 binary covariates of size n
#' @return List of 6 model.matrix outputs
make_Xmats <- function(n, formulas = NULL, covariate_data = NULL) {
  
  # Basic parameter checks
  #TODO(LCOMM): add ability to have unequal numbers of covariates?
  #TODO(LCOMM): add formula checking?

  # Make covariate data if none exists
  if (is.null(covariate_data)) {
    covariate_data <- simulate_covariate_data(n, p = 3)
  } else {
    stopifnot(is.data.frame(covariate_data))
    stopifnot(anyNA(covariate_data) == FALSE)
    stopifnot(NROW(covariate_data) == n)
  }
  
  # Assume all covariates in all formulas if none provided
  if (is.null(formulas)) {
    formulas <- replicate(6, formula(~ 1 + .))
  } else {
    stopifnot(is.list(formulas), length(formulas) == 6)
  }
  
  # Make Xmats list
  Xmats <- as.list(rep(NA, 6))
  for (i in 1:6) {
    Xmats[[i]] <- model.matrix(formulas[[i]], data = covariate_data)
  }
  
  # Return
  return(Xmats)
  
}


#' Function to generate frailty matrix
#' 
#' Columns correspond to hazard models, rows are individuals
#' Individuals are assumed independent
#' Frailties within a person are multivariate normal
#' @param n Number of observations
#' @param ftype Type of frailty pattern 
#' 0 = All 0
#' 1 = single, 
#' 2 = Tx-specific but shared within treatment arms
#' 3 = transition-specific but shared across treatment arms
#' 4 = Completely 
#' @param distn Distribution; default = "norm"
#' (required to be multivariate normal for now)
#' @param fSig User-provided variance-covariance matrix for frailty vector (optional); 
#' default is to set variances to 1 and any correlations to 0.5
#' @return n x 6 numeric matrix of frailties
simulate_frailties <- function(n, ftype = 1, distn = "norm", fSig = NULL) {
  
  # Require mvtnorm package
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("mvtnorm needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # Very basic parameter checking
  stopifnot(NROW(ftype) == 1, ftype %in% 0:4, distn == "norm")
  
  # Return all zeros if requested
  if (ftype == 0) {
    return(matrix(0, nrow = n, ncol = 6))
  }
  
  if (is.null(fSig)) {
    
    # Type 1: completely shared frailty
    fSig <- matrix(1, nrow = 6, ncol = 6)
    
    # Type 2: 2 treatment-specific frailties, potentially correlated
    # Default is a correlation of 0.5
    fSig <- matrix(0.5, nrow = 6, ncol = 6)
    fSig[1:3, 1:3] <- 1
    fSig[4:6, 4:6] <- 1
  
    # Type 3: 3 transition-specific frailties, potentially correlated
    # Default block is a compound symmetric with off-diagonal 0.5
    b <- matrix(0.5, nrow = 3, ncol = 3)
    diag(b) <- 1
    r <- cbind(b, b)
    fSig <- rbind(r, r)
    
    # Type 4: 6 hazard-specific frailties, all potentially correlated
    # Default block is a compound symmetric with off-diagonal 0.5
    fSig <- matrix(0.5, nrow = 6, ncol = 6)
    diag(fSig) <- 1
    
  } else {
    
    # Need matrixcalc package
    if (!requireNamespace("matrixcalc", quietly = TRUE)) {
      stop("matrixcalc needed for this function to work. Please install it.",
           call. = FALSE)
    }
    
    # Basic parameter checking on provided covariance matrix
    stopifnot(is.numeric(fSig), is.matrix(fSig), dim(fSig) == c(6,6),
              matrixcalc::is.positive.semi.definite(fSig))
    
  }
  
  # Return frailty matrix
  return(mvtnorm::rmvnorm(n, sigma = fSig))
  
}


#' Make vector of frailty coefficients
#'
#' Position is the model the coefficient belongs to
#' 1 = healthy-ill transition in control
#' 2 = healthy-dead transition in control
#' 3 = ill-dead transition in control
#' 4 = healthy-ill transition in treated
#' 5 = healthy-dead transition in treated
#' 6 = ill-dead transition in treated
#' @param fctype Frailty coefficient type 0-6; default is 1
#' Type 0: Underlying frailty is irrelevant
#' Type 1: All coeffients are identical, not 0
#' Type 2: Frailty coefficients are treatment-specific
#' Type 3: Frailty coefficients are transition-specific
#' Type 4: All frailty coefficients are different
#' @return Length 6 vector of frailty coefficients
give_frailty_coefs <- function(fctype = 1) {
  
  # Set frailty coefficients based on type
  betas <- switch(as.character(fctype),
                  "0" = rep(0, 6),
                  "1" = rep(0.2, 6),
                  "2" = c(rep(0.15, 3), rep(0.05, 3)),
                  "3" = c(c(0.1, 0.2, 0.3), c(0.1, 0.2, 0.3)),
                  "4" = seq(0.30, 0.05, length.out = 6))
  
  # Return
  return(betas)

}

#' Create list of regression and frailty coefficients 
#' 
#' @param Xmats List of 6 design matrices, of potentially varying column #s
#' @param fctype Type of the frailty coefficients; default is 1
#' see \code{\link{give_frailty_coefs}} for details
#' @param exp Whether
#' @return Named list of regression coefficients (rcoefs) and frailty 
#' coefficients (fcoefs) by arm ("control" and "treated")
give_coefs <- function(Xmats, fctype = 1, exp = TRUE) {
  
  # Extract number of columns for the design matrices
  ps <- sapply(Xmats, dim, simplify = TRUE)[2, ]
  
  # Regression coefficients, made so that treatment effect on logHR is -0.2
  # Corresponds to HR of ~0.81
  # Same for both non-terminal and terminal hazards
  # (Treatment effect is difference in intercepts since everything else is same)
  rcoefs <- sapply(ps, function(x) { c(-1, seq(0.1, 0.2, length.out = x - 1)) })
  rcoefs[1, 4:6] <- -1.2
  
  # Shape parameters for Weibull
  # Moderate deviation from Exponential if not forced to be Exponential
  shapes <- if (exp) { rep(1, 6) } else { c(0.95, 1.0, 1.05, 0.95, 1.0, 1.05) }
  
  # Frailty coefficients
  fcoefs <- give_frailty_coefs(fctype)
  
  # Package nicely
  coefs <- list()
  coefs[["control"]] <- list(rcoefs = rcoefs[, 1:3], fcoefs = fcoefs[1:3],
                             shapes = shapes[1:3])
  coefs[["treated"]] <- list(rcoefs = rcoefs[, 4:6], fcoefs = fcoefs[4:6],
                             shapes = shapes[4:6])
  
  # Return
  return(coefs)
  
}

#' Simulate Weibull potential outcomes in both exposure arms
#' 
#' @param Xmats List of 6 design matrices (first 3 control, 2nd 3 are tx)
#' @param coefs List of control and treated regression coefs and Weibull shapes
#' see output of \code{\link{give_coefs}} for details on input type
#' @param frailties n x 6 matrix of frailties
#' see output of \code{\link{simulate_frailties}} for details on input type
#' @return n x 8 dataframe with potential outcomes under control and treated
simulate_poutcomes <- function(Xmats, coefs, frailties) {

  # Very basic parameter checks
  stopifnot(is.list(Xmats), length(Xmats) == 6)
  stopifnot(is.list(coefs), names(coefs) == c("control", "treated"))
  stopifnot(is.matrix(frailties), NCOL(frailties) == 6)
  
  # Initialize data frame for potential outcomes
  n <- NROW(Xmats[[1]])
  dat <- data.frame(row.names = 1:n)
  
  for (arm_i in 1:2) {
    z <- arm_i - 1
    arm <- c("control", "treated")[arm_i]
    start_i <- c(1, 4)[arm_i]
    end_i <- start_i + 2
    arm_rcoefs <- coefs[[arm]][["rcoefs"]]
    arm_fcoefs <- coefs[[arm]][["fcoefs"]]
    arm_shapes <- coefs[[arm]][["shapes"]]
    arm_frailties <- frailties[, start_i:end_i]
    arm_Xmats <- Xmats[start_i:end_i]
    arm_lps <- matrix(NA, nrow = NROW(arm_Xmats[[1]]), ncol = 3)
    for (j in 1:3) {
      arm_lps[, j] <- arm_Xmats[[j]] %*% arm_rcoefs[, j] + arm_frailties[, j]
    }
    arm_scales <- exp(-arm_lps/arm_shapes)
    
    # Generate data
    R <- rweibull(n, shape = arm_shapes[1], scale = arm_scales[1])
    D <- rweibull(n, shape = arm_shapes[2], scale = arm_scales[2])
    soj <- rweibull(n, shape = arm_shapes[3], scale = arm_scales[3])
    yesR <- (R < D)
    D[yesR] <- R[yesR] + soj[yesR]
    deltaR <- yesR * 1
    deltaD <- rep(1, n) # assuming no censoring
    R[deltaR == 0] <- D[deltaR == 0]
    
    # TODO(LCOMM): Add censoring later
    
    # Save in data set
    arm_names <- paste0(c("R", "D", "deltaR", "deltaD"), z)
    dat[[arm_names[1]]] <- R
    dat[[arm_names[2]]] <- D
    dat[[arm_names[3]]] <- deltaR
    dat[[arm_names[4]]] <- deltaD
    
  }
  
  # Return
  return(dat)
  
}


#' Function to simulate data set with entire truth (frailties, po's, etc)
#' 
#' Generate a truth table (with both potential outcome sets) for an 
#' illness-death model
#' 
#' @param n Number of observations
#' @param ftype Frailty type
#' @param fctype Frailty coefficient type
#' @param exp Whether exponential (TRUE) or Weibull
#' @param p Number of covariates to simulate
#' @param formulas List of formulas for models
#' @return Data frame of n observations
simulate_entire_truth <- function(n = 100, ftype = 1, fctype = 1, exp = TRUE, 
                                  p = 3, formulas = NULL, ...) {
  
  # Make data
  cov_dat <- simulate_covariate_data(n = n, p = p)
  Xmats <- make_Xmats(n, formulas = formulas, covariate_data = cov_dat)
  coefs <- give_coefs(Xmats = Xmats, fctype = fctype, exp = exp)
  frailties <- simulate_frailties(n, ftype, ...)
  colnames(frailties) <- paste0("f", 1:6)
  dfrailties <- as.data.frame(frailties)
  dat <- simulate_poutcomes(Xmats, coefs, frailties)
  
  
  # Combine into single outcome data set
  dfrailties$id <- cov_dat$id <- dat$id <- 1:n
  dat <- merge(cov_dat, dat, by = "id")
  dat <- merge(dat, dfrailties, by = "id")
  
  # Return entire data set
  return(dat)
  
}




