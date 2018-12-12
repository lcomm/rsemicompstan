#' Make data application design matrix
#' 
#' @param df Data frame containing adjustment covariates
#' @return Design matrix for data application
#' @export
make_da_xmat <- function(df) {
  
  binary_vars <- c("race_noWhite", "sex_female", "deyo2_", "adm")
  continuous_vars <- c("age_std", "los_std")
  stopifnot(c(binary_vars, continuous_vars) %in% colnames(df))
  
  # Mean-center binary variables and center/scale continuous ones
  rescaled_df <- df
  for (bvar in binary_vars) {
    rescaled_df[[bvar]] <- scale(rescaled_df[[bvar]], center = TRUE, scale = FALSE)
  }
  for (cvar in continuous_vars) {
    rescaled_df[[cvar]] <- scale(rescaled_df[[cvar]], center = TRUE, scale = TRUE)
  }
  
  xmat <- rescaled_df[, c(binary_vars, continuous_vars)]
  return(xmat)
}



#' Make a new "scaled/centered" data 
#' 
#' Useful for making posterior predictive draws for a covariate pattern
#' 
#' @param reference_df Data frame for mean
#' @param age Age (default is 85)
#' @param los Index hospitalization length of stay (default is 10 days)
#' @param race_noWhite 0/1 indicator of non-White race
#' @param sex_female 0/1 indicator of being female
#' @param deyo2_ 0/1 indicator of comorbidity score > 2
#' @param adm 0/1 indicator of admission route
#' @return Scaled design matrix for new "data" with scaling done by as it
#' would have been in original data set
#' @export
make_new_da_xmat <- function(reference_df, 
                             age = c(85, 65), los = 0, 
                             race_noWhite = c(1, 0), sex_female = c(0, 1), 
                             deyo2_ = 0, adm = 0) {
  
  # Apply scaling done originally
  age_std <- ifelse(is.na(age), 0, (age - 85) / 5)
  los_std <- ifelse(is.na(los), 0, (los - 10) / 7)
  rescaled_df <- data.frame(age_std = age_std, los_std = los_std, 
                            race_noWhite = race_noWhite, sex_female = sex_female, 
                            deyo2_ = deyo2_, adm = adm)
  
  # Reverse the scaling and centering
  binary_vars <- c("race_noWhite", "sex_female", "deyo2_", "adm")
  continuous_vars <- c("age_std", "los_std")
  for (bvar in binary_vars) {
    # Add back on mean
    rescaled_df[[bvar]] <- rescaled_df[[bvar]] + mean(reference_df[[bvar]])
  }
  for (cvar in continuous_vars) {
    # Rescale by SD and add back on mean
    rescaled_df[[cvar]] <- rescaled_df[[cvar]] * sd(reference_df[[cvar]]) + 
                           mean(reference_df[[cvar]])
  }
  xmat <- rescaled_df[, c(binary_vars, continuous_vars)]
  return(xmat)
}



#' Fit the model for the data application propensity scores
#' 
#' @param dat Data frame containing z and adjustment covariates
#' @return GLM fit object
#' @export
fit_da_ps_model <- function(dat) {
  fit <- glm(z ~ race_noWhite + sex_female + deyo2_ + adm + 
               age_std + los_std + 
               age_std:age_std + I(age_std^2) + I(los_std^2),
             data = dat)  
  return(fit)
}



#' Fit basic model in data application with pancreatic data
#' 
#' Prior for sigma is weak (~ 20 prior obs) around frailty variance in Lee 2015.
#' Applies fuzzing to handle discretized time scale. 
#' 
#' @param file File (including path) to RDS with scrData object
#' @param use_priors Passed to \code{\link{scr_gamma_frailty_stan}}
#' @param ps_use Whether to trim controls from the sample based on 
#' propensity score overlap problems using \code{link{trim_da_with_ps}} ("trim"),
#' only use a sample with matching to treated ("match"), or do nothing at all ("none"). 
#' Defaults to "match"
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
                         ps_use = "match",
                         init = "random", 
                         init_r = 0.5,
                         sigma_pa = 21, 
                         sigma_pb = 7.1,
                         iter = 2000, 
                         chains = 4,
                         ...) {
  
  # Read in data
  scr_pc <- readRDS(file = file)
  
  # Alias variables
  scr_pc$yr <- scr_pc$time1
  scr_pc$yt <- scr_pc$time2
  scr_pc$dyr <- scr_pc$event1
  scr_pc$dyt <- scr_pc$event2
  
  # Fuzz so that events never happen at time zero
  # Add one-half day to sojourn times of zero
  not_zero_soj   <- which(!((scr_pc$yt == scr_pc$yr) & (scr_pc$dyr == 1)))
  zero_soj <- (1:NROW(scr_pc))[-not_zero_soj]
  scr_pc$yt[zero_soj] <- scr_pc$time2[zero_soj] + 0.5
  
  # Exclude anyone discharged to something besides home or home with care
  scr_pc <- scr_pc[pmax(scr_pc$disc_snf_icf,
                        scr_pc$disc_hospice,
                        scr_pc$disc_other) == 0, ]
  
  # Adjustment covariates are race, standardized age, sex, comorbidity score, 
  # admission route, and hospital length of (initial) stay
  # "Treatment" is being discharged to home with care (vs. home with no care)
  # This is reversal of home vs. not-home from before!
  scr_pc$z <- scr_pc$disc_homecare
  if (ps_use == "trim") {
    # Make design matrix for propensity score calculations
    ps_fit <- fit_da_ps_model(dat = scr_pc)
    scr_pc$ps <- fitted(ps_fit)
    min_treated_ps <- min(scr_pc$ps[scr_pc$z == 1])
    scr_pc$exclude <- ifelse((scr_pc$z == 0) & (scr_pc$ps < min_treated_ps),
                             1, 0)
    # Apply exclusion
    excluded <- scr_pc[scr_pc$exclude == 1, ]
    scr_pc   <- scr_pc[scr_pc$exclude == 0, ]
  } else if (ps_use == "match") {
    m_scr <- MatchIt::matchit(z ~ race_noWhite + sex_female + deyo2_ + adm + 
                                age_std + los_std + age_std:age_std + 
                                I(age_std^2) + I(los_std^2), 
                              data = scr_pc,
                              method = "nearest", 
                              discard = "control",
                              ratio = 1,
                              replace = FALSE)
    # Apply exclusion
    matched <- as.numeric(c(row.names(m_scr$match.matrix), 
                            m_scr$match.matrix[ , 1]))
    scr_pc$exclude <- 1
    scr_pc$exclude[row.names(scr_pc) %in% matched] <- 0
    excluded <- scr_pc[scr_pc$exclude == 1, ]
    scr_pc   <- scr_pc[scr_pc$exclude == 0, ]
  } else if (ps_use == "none") {
    excluded <- NULL
  }
  xmat <- make_da_xmat(df = scr_pc)
  
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
  return(list(dat = scr_pc, xmat = xmat, stan_fit = stan_fit, 
              excluded = excluded))
}



#' Make a block for the posterior prediction table for a single covariate
#' pattern and frailty type at a series of times
#' 
#' @param xnew Design matrix for new observations
#' @param stan_fit Stan fit object for making predictions
#' @param eval_t Length-K vector of times for evaluating predictions
#' @param frailty_q Quantile for frailty (Special case: 0 = mean = frailty of 1)
#' @param thin_out Number of MCMC iterations to base calculations on
#' @return K-row data frame of principal state probabilities and causal effects
ppred_block <- function(xnew, stan_fit, eval_t, frailty_q = 0.5, thin_out = 1000) {  
  if (frailty_q == 0) {
    ppnew <- posterior_predict_xnew(xnew = xnew, stan_fit = stan_fit, 
                                    frailty = rep(1, NROW(xnew)))
  } else {
    ppnew <- posterior_predict_xnew(xnew = xnew, stan_fit = stan_fit, 
                                    frailty_q = frailty_q)  
  }
  
  R <- dim(ppnew)[3]
  stopifnot(R >= thin_out)
  which_r <- round(seq(1, R, length.out = thin_out))
  ppnew <- ppnew[ , , which_r]
  
  f_dd <- f_ck <- f_tk <- f_aa <- eval_t * 0 
  for (r in 1:thin_out) {
    pstates <- v_make_pstates(eval_t = eval_t, pp = ppnew[ , , r])
    f_aa <- f_aa + colMeans(pstates == "AA") / thin_out
    f_ck <- f_ck + colMeans(pstates == "TS") / thin_out
    f_tk <- f_tk + colMeans(pstates == "TK") / thin_out
    f_dd <- f_dd + colMeans(pstates == "DD") / thin_out  
  }
  
  # Data frame containing proportion in each state at each t
  state_names <- c("AA", "CK", "TK", "DD")
  comp_dat <- data.frame(Time = rep(eval_t, times = length(state_names)),
                         State = rep(state_names, each = length(eval_t)))
  comp_dat$Proportion <- c(f_aa, f_ck, f_tk, f_dd)
  comp_dat$State <- factor(comp_dat$State, levels = rev(state_names), 
                           ordered = TRUE)
  res <- reshape2::dcast(data = comp_dat, Time ~ State, 
                         value.var = "Proportion")
  res$Time <- eval_t
  res <- res[, c("Time", "AA", "CK", "TK")]
  res$TVSACE <- rowMeans(apply(ppnew, 3, FUN = v_calculate_tv_sace, 
                               eval_t = eval_t), na.rm = TRUE)
  res$RMSACE <- rowMeans(apply(ppnew, 3, FUN = v_calculate_rm_sace, 
                               eval_t = eval_t), na.rm = TRUE)
  return(res)
}



#' Make posterior predictive table for data application
#' 
#' @param stan_res Stan fit object from \code{\link{fit_data_app}}
#' @param eval_t Times for predictions/causal effect calculations
#' @param frailty_q Quantiles for frailty used in predictions; 0 = average frailty of 1
#' @param thin_out Thinning of MCMC iterations for calculation across posterior draws. If 
#' \code{thin_out == 0}, all MCMC iterations are used.
#' @param nrep Number of Monte Carlo replicates for each covariate/frailty pattern
#' @return Knitr kable table that can be passed to \code{\link{make_ppred_table_pretty}}
#' @export
ppred_table <- function(stan_res,
                        eval_t = c(30, 90), 
                        frailty_q = c(0.9, 0, 0.1),
                        thin_out = 0, nrep = 10000) {
  # Alias
  reference_df <- stan_res$dat
  stan_fit     <- stan_res$stan_fit
  if (thin_out == 0) {
    thin_out <- length(rstan::extract(stan_fit, par = "lp__")[[1]])
  }
  
  # Make list of data frame from desired covariate patterns
  list_xnew <- list(make_new_da_xmat(reference_df = reference_df, 
                                     age = rep(85, nrep),
                                     race_noWhite = 1, sex_female = 0,
                                     los = 0, deyo2_ = 0, adm = 0),
                    make_new_da_xmat(reference_df = reference_df,
                                     age = rep(65, nrep),
                                     race_noWhite = 0, sex_female = 1,
                                     los = 0, deyo2_ = 0, adm = 0))
  
  # Make results shell
  ncovs <- length(list_xnew)
  nfrails <- length(frailty_q)
  nt <- length(eval_t)
  res <- rev(expand.grid(RMSACE = NA, TVSACE = NA, TK = NA, CK = NA, AA = NA,
                         Time = eval_t, Frailty = 1:nfrails, Pattern = 1:ncovs))
  
  # Loop over blocks
  nblocks <- nfrails * ncovs
  for (block_i in 1:nblocks) {
    row_start <- (block_i - 1) * 2 + 1
    row_end <- row_start + length(eval_t) - 1
    res_block <- ppred_block(xnew = list_xnew[[res$Pattern[row_start]]],
                             stan_fit, 
                             eval_t, 
                             frailty_q = frailty_q[res$Frailty[row_start]],
                             thin_out = thin_out) 
    res[row_start:row_end, 3:ncol(res)] <- res_block
  }
  
  # Return
  return(res)
}



#' Pretty up the posterior predictive table for data application
#' 
#' @param tab Kable output from \code{\link{ppred_table}}
#' @param digits Number of digits to pass to knitr::kable
#' @param caption Table caption
#' @return Kabled table with prettier labels and nicer column headings
#' @export
make_ppred_table_pretty <- function(tab, digits = 3, 
                                    caption = "Posterior  predictive  means  for  
                                    principal  state  probabilities  and  
                                    principal  stratum  causal effects for new 
                                    patients of two covariate patterns") {
  
  # Extract and do basic checks that table seems right dimensions
  nfrails <- length(unique(tab$Frailty))
  stopifnot((nfrails == 3), (length(unique(tab$Pattern)) == 2))
  nt <- length(unique(tab$Time))
  
  # Remove some of the duplicate row contents
  tab$Frailty[seq(2, NROW(tab), by = nt)] <- nfrails + 1
  lh_desc <- c("Frail", "Average", "Healthy")
  slh <- paste0("\\multirow[t]{", nt, "}{*}{")
  elh <- "}"
  full_lh <- paste0(slh, lh_desc, elh)
  tab$Frailty <- factor(tab$Frailty,
                        levels = c(1:4),
                        labels = c(full_lh, ""))
  xnew_desc <- c("Nonwhite male aged 85, \\\\ average comorbidity score \\\\ and hospital length of stay",
                 "White female aged 65, \\\\ average comorbidity score \\\\ and hospital length of stay")
  spc <- paste0("\\multirow[t]{", nfrails * nt, "}{*}{\\begin{tabular}[t]{@{}l@{}}")
  epc <- paste0("\\end{tabular}}")
  full_pc <- paste0(spc, xnew_desc, epc)
  tab$Pattern[-seq(1, NROW(tab), by = nt * nfrails)] <- 3
  tab$Pattern <- factor(tab$Pattern,
                        levels = c(1:3),
                        labels = c(full_pc, ""))
  tv <- "\\begin{tabular}[c]{@{}c@{}}Difference in readmission \\\\ incidence by $t$ \\end{tabular}"
  rm <- "\\begin{tabular}[c]{@{}c@{}}Additional readmission-free \\\\ days accumulated by $t$ \\end{tabular}"
  tabcolnames <- c("Patient characteristics", "Latent health\\tnote{2}",
                   "Day $t$", "AA", "CK", "TK", tv, rm)
  
  # Do basic kabling and some minor tweaks to add footnotes, etc
  ktab <- knitr::kable(tab, format = "latex", booktabs = TRUE, 
                       digits = digits, escape = FALSE, col.names = tabcolnames, 
                       caption = caption)
  ktab <- gsub(pattern = "\\addlinespace", x = ktab, replacement = "", 
               fixed = TRUE)
  ktab <- gsub(pattern = "\\begin{tabular}{llrrrrrr}\n\\toprule", x = ktab,
               replacement = "\\scriptsize \\begin{threeparttable}[t] \\begin{tabular}{llcccccc}\\toprule
               \\multicolumn{1}{c}{} & \\multicolumn{1}{c}{} &  & 
               \\multicolumn{3}{c}{\\begin{tabular}[c]{@{}c@{}}Principal State \\\\ 
               Probabilities at $t$\\tnote{1}\\end{tabular}} & 
               \\multicolumn{2}{c}{\\begin{tabular}[c]{@{}c@{}}If always-alive, \\\\
               causal effect of being \\\\ 
               discharged to home with support \\\\ (vs. without) \\end{tabular}} 
               \\\\ \\midrule",
               fixed = TRUE)
  ktab <- gsub(pattern = "\\bottomrule\n\\end{tabular}", x = ktab,
               replacement = "\\bottomrule\\end{tabular}\\begin{tablenotes}
                \\item[1] Always-alive (AA), dead only under control (CK), and dead only under treatment (TK)
                \\item[2] Frail and healthy correspond to the $90^{th}$ and
               $10^{th}$ percentiles of
               $\\gamma$, while average health corresponds to $\\gamma = 1$
               \\end{tablenotes}
               \\end{threeparttable}",
               fixed = TRUE)
  ktab <- gsub(pattern = "\\label{tab:}", x = ktab,
               replacement = "\\label{tab:pps}", fixed = TRUE)
  return(ktab)
}



#' Plot propensity score densities overlapping one another
#' 
#' @param dat Data frame with propensity scores in \code{ps}
#' @return ggplot object
#' @export
plot_ps_overlap <- function(dat) {
  p <- ggplot(data = dat, 
              aes_string(x = "ps", group = "as.factor(z)", 
                         fill = "as.factor(z)")) + 
    geom_density(alpha = 0.5) + 
    xlim(c(0,1)) + 
    labs(x = "Propensity score", y = "Density", 
         fill = "Treated")
  return(p)
}



#' Make data for graphing discrepancy measures for data application
#' 
#' @param stan_res Object resulting from \code{\link{fit_data_app}}
#' @param seed Random number see for posterior predictive draws
#' @param eval_t Length-K vector of times to calculate statistic at
#' @param subsamp Whether to subsample the MCMC iterations for calculating
#' discrepancy metrics. If \code{FALSE}, all MCMC iterations are used. 
#' @return Data frame that can be reshaped/plotted
#' @export
prepare_discrep_plot_dat <- function(stan_res, seed, 
                                     eval_t = c(15, 30, 45, 60, 90),
                                     subsamp = FALSE) {
  kms <- calculate_km_disc(eval_t = eval_t, stan_res = stan_res, 
                           cens_times = rep(90, NROW(stan_res$dat)), 
                           seed = seed, subsamp = subsamp)
  aas <- calculate_frac_aa_disc(eval_t = eval_t, stan_res = stan_res, 
                                cens_times = rep(90, NROW(stan_res$dat)), 
                                seed = seed, subsamp = subsamp)
  disc_plot_dat <- as.data.frame(cbind(kms, AA = aas, eval_t = eval_t))
  disc_plot_dat <- reshape2::melt(data = disc_plot_dat, 
                                  id.vars = "eval_t", 
                                  measure.vars = c("KM0", "KM1", "AA"),
                                  value.name = "disc_test_stat",
                                  variable.name = "measure_type")
  return(disc_plot_dat)
}



#' Function to plot discrepancy statistics
#' 
#' @param stan_res Object resulting from \code{\link{fit_data_app}}
#' @param seed Random number see for posterior predictive draws
#' @param eval_t Length-K vector of times to calculate statistic at
#' @param subsamp Whether to subsample the MCMC iterations for calculating
#' discrepancy metrics. If \code{FALSE}, all MCMC iterations are used. 
#' @param color_vals Color values for plot (order: KM0, KM1, AA)
#' @return ggplot object
#' @export
make_da_discrepancy_plot <- function(stan_res, seed, 
                                     eval_t = c(15, 30, 45, 60, 90), 
                                     subsamp = 100,
                                     color_vals = c("#0B353B",
                                                    "#EC3F19",
                                                    "#2E655D")) {
  plot_dat <- prepare_discrep_plot_dat(stan_res = stan_res, seed = seed, 
                                       eval_t = eval_t,
                                       subsamp = subsamp)
  measure_labels <- c(bquote(T[KM*",0"]), bquote(T[KM*",1"]), bquote(T[AA]))
  p <- ggplot(plot_dat,
              aes_string(y = "disc_test_stat", 
                         x = "eval_t", 
                         group = "measure_type", 
                         color = "measure_type",
                         linetype = "measure_type")) +
    
    scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) + 
    scale_x_continuous(breaks = c(0, eval_t), limits = c(0, 90)) + 
    scale_linetype_manual(values = c("dotted", "dashed", "twodash"),
                          labels = measure_labels) + 
    scale_color_manual(values = color_vals, 
                       labels = measure_labels) + 
    geom_line() +
    geom_point(size = 2) + 
    labs(color = "Measure", linetype = "Measure", 
         x = "Day", 
         y = "Test statistic", 
         title = "Discrepancy measure test statistics")
  return(p)
}



#' Compare imputed in-sample frailty density to theoretical density implied by 
#' frailty variance sigma
#' 
#' Used as diagnostic for AA discrepancy metric
#' 
#' @param stan_fit Stan fit from data application
#' @param pp Posterior prediction array from \code{\link{posterior_predict_sample}}
#' @return ggplot2 plot object
#' @export
make_da_frailty_density <- function(stan_fit, pp, 
                                    color_vals = c("#800026", "#0B353B")) {
  
  # Extract imputed frailties in observed data set
  frailties <- pp[ , "frailty", ]
  
  # Sample equivalent number of replicate frailties from the theoretical 
  # distribution based on posterior for sigma
  sigmas <- unlist(rstan::extract(stan_fit, par = "sigma"))
  R <- length(sigmas)
  length_out <- dim(pp)[3]
  thinning_r <- round(seq(1, R, length.out = length_out))
  nrep <- dim(pp)[1]
  gamma_rep <- sapply(sigmas[thinning_r], 
                      FUN = function(x) rgamma(n = nrep, 1 / x, 1 / x))
  dat <- data.frame(type = rep(c("In-sample", "Theoretical"), 
                               each = nrep * length_out),
                    frailties = c(frailties, gamma_rep))
  
  # Make plot
  p <- ggplot(dat, aes_string(x = "frailties", color = "type", fill = "type")) + 
    geom_density(alpha = 0.3) + 
    scale_color_manual(values = color_vals) + 
    scale_fill_manual(values = color_vals) + 
    xlim(c(0, 4)) + 
    labs(x = "Frailty", y = "Density",
         fill = "Frailty Type", 
         color = "Frailty Type") 
  
  return(p)
}

