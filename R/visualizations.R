#' Get the median time for coloring ridges based on P(AA)
#' 
#' Calculates time at which half of referent-hazard (i.e., frailty = 1) 
#' individuals will have experienced at least one of the event times, assuming 
#' all covariates are centered.
#' 
#' @param scenario Data generation scenario
#' @param P Number of covariates for DGP. Defaults to 4.
#' @param max_t Maximum time (used only for upper bounding rootfinding)
#' @return Scalar time at which to evaluate the probability of being an 
#' always-survivor
#' @export
get_color_t <- function(scenario, P = 4, max_t = get_eval_t()[2]) {
  
  f_median <- function(x, a1, b1, a2, b2) {
    log(0.5) + (x / b1)^a1 + (x / b2)^a2
  }
  
  params <- return_dgp_parameters(scenario, P = P)
  a1 <- params$alpha1
  b1 <- make_scale(lp = log(params$control$kappa1), 
                   alpha = params$control$alpha1)
  a2 <- params$control$alpha2
  b2 <- make_scale(lp = log(params$control$kappa2), 
                   alpha = params$control$alpha2)
  
  color_t <- uniroot(f_median, lower = 0, upper = max_t * 2, 
                     tol = 1e-10, 
                     a1 = a1, b1 = b1, 
                     a2 = a2, b2 = b2)$root
  
  return(color_t)
}



#' Calculate fraction of population that is always-alive
#' 
#' @param eval_t Time at which to evaluate survival
#' @param pp Posterior predictive draws (for a single MCMC iteration)
#' @return Scalar value between 0 and 1 representing proportion always-alive
#' @export
calculate_frac_aa <- function(eval_t, pp) {
  pp <- as.data.frame(pp)
  frac_aa <- mean(make_pstates(eval_t, pp) == "AA")
  return(frac_aa)
}



#' Vectorized calculate fraction of population that is always-alive
#' 
#' Evaluates for a vector of times
#' 
#' @param eval_t Times at which to evaluate survival
#' @param pp Posterior predictive draws (for a single MCMC iteration)
#' @return Scalar value between 0 and 1 representing proportion always-alive
#' @export
v_calculate_frac_aa <- Vectorize(FUN = calculate_frac_aa, 
                                 vectorize.args = "eval_t")



#' Calculate fraction of population that has event time past some time
#' 
#' Assumes no censoring
#' 
#' @param eval_t Time at which to evaluate survival
#' @param pp Posterior predictive draws (for a single MCMC iteration)
#' @return Scalar value between 0 and 1 representing proportion alive
#' @export
calculate_frac_alive <- function(eval_t, pp, eventvar) {
  pp <- as.data.frame(pp)
  frac_alive <- mean(pp[ , eventvar] > eval_t)
  return(frac_alive)
}



#' Vectorized calculate fraction of population with event time > t
#' 
#' Evaluates for a vector of times
#' 
#' @param eval_t Times at which to evaluate survival
#' @param pp Posterior predictive draws (for a single MCMC iteration)
#' @param eventvar Name of event time variable
#' @return Scalar value between 0 and 1 representing proportion always-alive
#' @export
v_calculate_frac_alive <- Vectorize(FUN = calculate_frac_alive, 
                                    vectorize.args = "eval_t")



#' Vectorized calculation of principal states
#' 
#' Evaluates for a vector of times
#' 
#' @param eval_t Vector of T times at which to evaluate survival
#' @param pp Posterior predictive draws (for a single MCMC iteration)
#' @return N x T matrix of principal state characters ("AA"/"TK"/"CK"/"DD")
#' @export
v_make_pstates <- Vectorize(FUN = make_pstates, 
                            vectorize.args = "eval_t")



#' Turn posterior predictive data into a long format amenable to graphing the 
#' causal effect
#' 
#' @param pp Posterior prediction array
#' @param max_t Maximum time for calculating effects
#' @param length_out Number of time points at which to calculate effects. More points
#' means smoother effect lines
#' @return Data frame with MCMC iteration (r), evaluation time point (eval_t),
#' time-varying survivor average causal effect (tv_sace), and restricted mean survivor
#' average causal effect (rm_sace), proportion always-alive (frac_aa), fraction alive under
#' treatment (frac_a_t), and alive under control (frac_a_c)
#' @export
prepare_graph_data <- function(pp, max_t, length_out = 10) {
  R <- dim(pp)[3]
  xt <- seq(0, max_t, length.out = length_out)
  res <- as.data.frame(expand.grid(r = 1:R, eval_t = xt, 
                                   frac_aa = NA, frac_a_t = NA, frac_a_c = NA,
                                   tv_sace = NA, rm_sace = NA))
  
  for (r in 1:R) {
    pp_mcmc <- as.data.frame(pp[, , r])
    for (eval_t in xt) {
      res[(res$r == r) & (res$eval_t == eval_t), "frac_a_t"] <- 
        calculate_frac_alive(eval_t = eval_t, pp = pp_mcmc, eventvar = "yt1_imp")
      res[(res$r == r) & (res$eval_t == eval_t), "frac_a_c"] <- 
        calculate_frac_alive(eval_t = eval_t, pp = pp_mcmc, eventvar = "yt0_imp")
      res[(res$r == r) & (res$eval_t == eval_t), "frac_aa"] <- 
        calculate_frac_aa(eval_t = eval_t, pp = pp_mcmc)
      res[(res$r == r) & (res$eval_t == eval_t), "tv_sace"] <- 
        calculate_tv_sace(eval_t = eval_t, pp = pp_mcmc)
      res[(res$r == r) & (res$eval_t == eval_t), "rm_sace"] <- 
        calculate_rm_sace(eval_t = eval_t, pp = pp_mcmc)
    }
  }
  
  return(res)
}



#' Calculate TV-SACE(r, t) among a group already subsetted to always-alive at t
#' 
#' Note: by assuming the data have already been subsetted to be always-alive
#' based on some t > max(eval_t), it is assumed than any times for r0 or r1 are 
#' event times. This may be problematic for r = t.
#' 
#' @param r1 Potential outcomes under treatment (assumed not censoring!)
#' @param r0 Potential outcomes under control (assumed not censoring!)
#' @param eval_t Scalar value for r argument in TV-SACE(r, t)
#' @export
tv_sace_aa_only <- function(r1, r0, eval_t) {
  r1_by_t <- r0_by_t <- rep(0, length(r1))
  r0_by_t[r0 < eval_t] <- 1
  r1_by_t[r1 < eval_t] <- 1
  diff_by_t <- r1_by_t - r0_by_t
  tv_sace <- mean(diff_by_t)
  return(tv_sace)
}



#' Vectorized calculation of TV-SACE when passed always-alive subset
#' 
#' Evaluates for a vector of times
#' 
#' @param r1 Potential outcomes under treatment (assumed not censoring!)
#' @param r0 Potential outcomes under control (assumed not censoring!)
#' @param eval_t Vector values for r argument in TV-SACE(r, t)
#' @export
v_tv_sace_aa_only <- Vectorize(tv_sace_aa_only, 
                               vectorize.args = "eval_t")



#' Calculate RM-SACE(r, t) among a group already subsetted to always-alive at t
#' 
#' Note: by assuming the data have already been subsetted to be always-alive
#' based on some t > max(eval_t), it is assumed than any times for r0 or r1 are 
#' event times. This may be problematic for r = t.
#' 
#' @param r1 Potential outcomes under treatment (assumed not censoring!)
#' @param r0 Potential outcomes under control (assumed not censoring!)
#' @param eval_t Scalar value for r argument in RM-SACE(r, t)
#' @export
rm_sace_aa_only <- function(r1, r0, eval_t) {
  diff_by_t <- pmin(eval_t, r1) - pmin(eval_t, r0)
  rm_sace <- mean(diff_by_t)
  return(rm_sace)
}



#' Vectorized calculation of RM-SACE when passed always-alive subset
#' 
#' Evaluates for a vector of times
#' 
#' @param r1 Potential outcomes under treatment (assumed not censoring!)
#' @param r0 Potential outcomes under control (assumed not censoring!)
#' @param eval_t Vector values for r argument in RM-SACE(r, t)
#' @export
v_rm_sace_aa_only <- Vectorize(rm_sace_aa_only, 
                               vectorize.args = "eval_t")



#' Turn posterior predictive data into a long format amenable to graphing the 
#' causal effects separately by always-alive cohorts
#' 
#' @param pp Posterior prediction array
#' @param cohort Time(s) for defining the always-alive cohort
#' @param length_out Number of time points at which to calculate effects. More points
#' means smoother effect lines
#' @return Data frame with MCMC iteration (r), evaluation time point (eval_t),
#' time-varying survivor average causal effect (tv_sace), and restricted mean survivor
#' average causal effect (rm_sace)
#' @export
prepare_cohort_graph_data <- function(pp, cohort, by = 1) {
  R <- dim(pp)[3]
  res <- as.data.frame(expand.grid(eval_t = seq(0, max(cohort), by = by), 
                                   cohort = cohort, 
                                   tv_sace = NA,
                                   rm_sace = NA,
                                   frac_aa = NA))
  res <- res[res$eval_t <= res$cohort, ]
  for (cohort_t in cohort) {
    xt <- seq(0, cohort_t, by = by)
    rm_sace <- tv_sace <- xt * 0
    frac_aa <- 0
    for (r in 1:R) {
      pp_mcmc <- as.data.frame(pp[, , r])
      pstates <- make_pstates(eval_t = cohort_t, pp = pp_mcmc)
      if (any(pstates == "AA")) {
        aa <- pp_mcmc[pstates == "AA", ]
        tv_sace <- tv_sace + v_tv_sace_aa_only(r1 = aa$yr1_imp, 
                                               r0 = aa$yr0_imp, 
                                               eval_t = xt) / R 
        rm_sace <- rm_sace + v_rm_sace_aa_only(r1 = aa$yr1_imp, 
                                               r0 = aa$yr0_imp, 
                                               eval_t = xt) / R 
        frac_aa <- frac_aa + (NROW(aa) / NROW(pp_mcmc)) / R
      } else {
        stop("Division by R not warranted!")
      }
    }
    res$eval_t[res$cohort == cohort_t]  <- xt
    res$tv_sace[res$cohort == cohort_t] <- tv_sace
    res$rm_sace[res$cohort == cohort_t] <- rm_sace
    res$frac_aa[res$cohort == cohort_t] <- frac_aa
  }
  res$cohort_id <- as.factor(res$cohort)
  return(res)
}



#' Make always-alive survival plot
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot object
#' @export
make_aa_kmplot <- function(plot_dat, time_unit = "Time") {

  # Set alpha to lower if many replicates
  alpha_val <- ifelse(length(unique(plot_dat$r)) < 100, 0.15, 0.0085)

  # Calculate overall mean 
  f_aa_mean_dat <- aggregate(frac_aa ~ eval_t, data = plot_dat, FUN = mean)
  f_aa_mean_dat$r <- 1
  
  # Make plot
  p <- ggplot(data = plot_dat,
              aes_string(y = "frac_aa", x = "eval_t", group = "r")) + 
    geom_line(alpha = alpha_val) + 
    ylim(0, 1) + 
    geom_line(data = f_aa_mean_dat,
              aes_string(y = "frac_aa", x = "eval_t")) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    labs(x = paste(tools::toTitleCase(time_unit), "r"), 
         y = "P(Always-alive at t)")
  
  return(p)
}



#' Make always-alive survival plot
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param frac_var Fraction variable name (frac_a_t or frac_a_c)
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot object
#' @export
make_pp_z_kmplot <- function(plot_dat, frac_var, time_unit = "Time") {
  
  # Alias
  plot_dat$frac_a <- plot_dat[[frac_var]]
  
  # Set alpha to lower if many replicates
  alpha_val <- ifelse(length(unique(plot_dat$r)) < 100, 0.15, 0.0085)
  
  # Calculate overall mean 
  f_a_mean_dat <- aggregate(frac_a ~ eval_t, data = plot_dat, FUN = mean)
  f_a_mean_dat$r <- 1
  
  # Make plot
  p <- ggplot(data = plot_dat,
              aes_string(y = "frac_a", x = "eval_t", group = "r")) + 
    geom_line(alpha = alpha_val) + 
    ylim(0, 1) + 
    geom_line(data = f_a_mean_dat,
              aes_string(y = "frac_a", x = "eval_t")) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    labs(x = paste(tools::toTitleCase(time_unit), "t"), 
         y = "P(Alive at t)")
  
  return(p)
}



#' Make combined Kaplan-Meiers for AA, z = 0, and z = 1
#' 
#' @param plot_dat Plot data object
#' @param color_vals Colors for lines (order: AA, z=0, z=1)
#' @param legend_position Position for legend. Defaults to "none"
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot2 object
#' @export
make_kmplot_combined <- function(plot_dat, 
                                 color_vals = c("#2E655D", "#0B353B", "#EC3F19"),
                                 legend_position = "none",
                                 time_unit = "Time") {

  # Set alpha to lower if many replicates
  alpha_val <- ifelse(length(unique(plot_dat$r)) < 100, 0.15, 0.0085)
  
  # Calculate overall mean for each (currently not used)
  # f_mean_dat <- reshape2::melt(aggregate(cbind(frac_aa, frac_a_c, frac_a_t) ~
  #                                          eval_t,
  #                                        data = plot_dat, FUN = mean),
  #                              id.vars = "eval_t",
  #                              variable.name = "surv_type",
  #                              value.name = "survival")
  # f_mean_dat$r <- as.numeric(f_mean_dat$surv_type)
  plot_dat <- plot_dat[, colnames(plot_dat) %in% c("eval_t", "r", "frac_aa",
                                                   "frac_a_c", "frac_a_t")]
  plot_dat_long <- reshape2::melt(plot_dat,
                                  id.vars = c("eval_t", "r"),
                                  variable.name = "surv_type",
                                  value.name = "survival")
  plot_dat_long$r <- as.factor(plot_dat_long$r)
  
  # Make plot
  p <- ggplot(data = plot_dat_long,
              aes_string(y = "survival", x = "eval_t", 
                         color = "surv_type", 
                         group = "interaction(r, surv_type)")) + 
    guides(group = FALSE) + 
    geom_line(alpha = alpha_val) + 
    ylim(0, 1) + 
    scale_color_manual("Survival Type", 
                       values = color_vals,
                       labels = c("Always-alive", 
                                  "Under z = 1", 
                                  "Under z = 0")) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) + 
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = paste(tools::toTitleCase(time_unit), "t"), y = "S(t)") + 
    theme(legend.position = legend_position)
  return(p)
}



#' Make TV-SACE curve plot (where r = t)
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param legend_position Position for legend. Defaults to "none"
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot object
#' @export
make_tvsace_plot <- function(plot_dat, 
                             legend_position = "none", 
                             time_unit = "Time") {
  tv_mean_dat <- aggregate(tv_sace ~ eval_t, data = plot_dat, FUN = mean)
  
  p <- ggplot(plot_dat, aes_string(x = "eval_t", y = "tv_sace", group = "r",
                                   color = "frac_aa")) + 
    geom_line(alpha = 0.4) + 
    geom_line(data = tv_mean_dat, 
              aes_string(x = "eval_t", y = "tv_sace", group = NULL), 
              color = "black",
              size = 0.7) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits = c(0, 1),
                          breaks = c(0, 0.5, 1),
                          labels = c("0%", "50%", "100%")) +
    guides(color = guide_colourbar(title = "Percent of \npopulation",
                                   title.hjust = 0.5,
                                   label.position = "left")) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          legend.text.align = 0.5) +
    labs(x = paste(tools::toTitleCase(time_unit), "t"), 
         y = "TV-SACE(t, t)")
  return(p)
}



#' Make RM-SACE curve plot (where r = t)
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param legend_position Position for legend. Defaults to "none"
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot object
#' @export
make_rmsace_plot <- function(plot_dat, 
                             legend_position = "none", 
                             time_unit = "Time") {
  rm_mean_dat <- aggregate(rm_sace ~ eval_t, data = plot_dat, FUN = mean)
  
  p <- ggplot(plot_dat, aes_string(x = "eval_t", y = "rm_sace", group = "r",
                                   color = "frac_aa")) + 
    geom_line(alpha = 0.2) +
    geom_line(data = rm_mean_dat, 
              aes_string(x = "eval_t", y = "rm_sace", group = NULL), 
              color = "black",
              size = 0.7) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits = c(0, 1),
                          breaks = c(0, 0.5, 1),
                          labels = c("0%", "50%", "100%")) +
    guides(color = guide_colourbar(title = "Percent of \npopulation",
                                   title.hjust = 0.5,
                                   label.position = "left"),
           alpha = FALSE) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          legend.text.align = 0.5) +
    labs(x = paste(tools::toTitleCase(time_unit), "t"),
         y = "RM-SACE(t, t)")
  return(p)
}



#' Make the TV-SACE plot with always-alive survival curve under it
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param aspect_ratio Shared aspect ratio for plots
#' @return ggplot object
#' @export
make_tvsace_plotpair <- function(plot_dat, aspect_ratio = 0.7) {
  p_top <- make_tvsace_plot(plot_dat) +
             theme(axis.title.x = element_blank(),
                   aspect.ratio = aspect_ratio)
  p_bot <- make_aa_kmplot(plot_dat) + 
             theme(aspect.ratio = aspect_ratio)
  p <- gridExtra::grid.arrange(p_top, p_bot, ncol = 1)
  return(p)
}



#' Make the RM-SACE plot with always-alive survival curve under it
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param aspect_ratio Shared aspect ratio for plots
#' @return ggplot object
#' @export
make_rmsace_plotpair <- function(plot_dat, aspect_ratio = 0.7) {
  p_top <- make_rmsace_plot(plot_dat) +
             theme(axis.title.x = element_blank(),
                   aspect.ratio = aspect_ratio)
  p_bot <- make_aa_kmplot(plot_dat) + 
             theme(aspect.ratio = aspect_ratio)
  p <- gridExtra::grid.arrange(p_top, p_bot, ncol = 1)
  return(p)
}



#' Make both TV-SACE and RM-SACE plots with always-alive survival curve under 
#' them
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @param aspect_ratio Shared aspect ratio for plots
#' @return ggplot object
#' @export
make_ce_plottrio <- function(plot_dat, aspect_ratio = 0.7) {
  p1 <- make_tvsace_plot(plot_dat) +
    theme(axis.title.x = element_blank(),
          aspect.ratio = aspect_ratio,
          legend.position = "right")
  p2 <- make_rmsace_plot(plot_dat) +
    theme(axis.title.x = element_blank(),
          aspect.ratio = aspect_ratio,
          legend.position = "none")
  p3 <- make_aa_kmplot(plot_dat) + 
    theme(aspect.ratio = aspect_ratio) + 
    ggtitle("Size of always-alive principal stratum")
  p <- gridExtra::grid.arrange(p1, p2, p3, ncol = 1)
  return(p)
}



#' Make ridgeplot colored by probability of stratum membership
#' 
#' Useful to compare sharp and blunt cases
#' 
#' @param res Result object from scenario run
#' @param risk_t Time at which to order individuals by probability of 
#' double-survival
#' @param color_t Time for coloring ridges by probability of double-survival
#' @return Combined ggplot object
#' @export
make_frailty_ridgeplot <- function(res, risk_t = get_eval_t()[1], 
                                   color_t = get_eval_t()[1]) {
    
    # Extract components of results object
    stan_fit <- res$stan_fit
    yr <- res$dat$yr
    yt <- res$dat$yt
    dyr <- res$dat$dyr
    dyt <- res$dat$dyt
    z <- res$dat$z
    xmat <- res$xmat
    
    # Make risk scores
    risks <- make_terminal_risk_scores(stan_fit, yr, yt, dyr, dyt, z, xmat, 
                                       eval_t = risk_t)
    n_sub <- 20
    subsamp <- order(risks)[floor(seq(1, NROW(xmat), length.out = n_sub))]
    pp <- posterior_predict_sample(stan_fit, yr, yt, dyr, dyt, z, xmat)
    
    df1 <- reshape2::melt(pp[ ,"frailty", ], varnames = c("id", "iter"),
                          value.name = "frailty")
    df_sub <- df1[df1$id %in% c(subsamp), ]
    df_sub$id <- factor(df_sub$id, levels = subsamp)
    
    # Get P(V = AA) for each individual averaged across MCMC
    is_aa_color_t <- (apply(pp, MARGIN = 3, FUN = make_pstates, 
                            eval_t = color_t) == "AA")
    sub_tx <- data.frame(id = as.factor(subsamp),
                         z = res$dat[subsamp, "z"])
    sub_tx$paa_t1 <- rowMeans(is_aa_color_t)[subsamp]
    df_sub <- merge(df_sub, sub_tx, by = "id")
    df_sub$risk <- factor(df_sub$id, levels = rev(levels(df_sub$id)))
    
    # Make ridgeplot
    p <- ggplot(df_sub,
                aes_string(x = "frailty", y = "risk", fill = "paa_t1")) + 
           ggridges::geom_density_ridges(rel_min_height = 0.05) + 
           scale_fill_gradientn("P(V = AA)",
                                colors = RColorBrewer::brewer.pal(11, "Spectral"),
                                limits = c(0, 1),
                                breaks = c(0, 1)) +
           ggridges::theme_ridges(center_axis_labels = TRUE, grid = FALSE) +
           scale_x_continuous(limits = c(0, 3.5)) +
           scale_y_discrete(expand = expand_scale(mult = c(0.01, 0.1)),
                            labels = c("Lowest", rep("", n_sub - 2), "Highest")) +
           labs(title = "Posterior distributions for the frailty for
                subsample of 20 individuals",
                x = "Frailty",
                y = "Covariate risk score")
    
  return(p)
}



#' Make side-by-side ridgeplots of sharp and blunt cases
#' 
#' @param res_sharp Result object from a sharp scenario replicate
#' @param res_blunt Result object from a blunt scenario replicate
#' @param titles Title labels for the two plots
#' @param color_t Vector of times at which to evaluate probabilities
#' of always-alive for coloring ridges
#' @return Combined ggplot object
#' @export
make_sharp_blunt_plotpair <- function(res_sharp, res_blunt,
                                      titles = c("Sharp", "Blunt"),
                                      color_t = c(get_color_t(scenario = 2),
                                                  get_color_t(scenario = 3))) {
  p_s <- make_frailty_ridgeplot(res_sharp, color_t = color_t[1])
  p_b <- make_frailty_ridgeplot(res_blunt, color_t = color_t[2])
  
  p <- ggarrange(p_s + theme(legend.position = "none") + ggtitle(titles[1]),
                 p_b + theme(legend.position = "none",
                             # axis.title.y = element_blank(),
                             # axis.text.y = element_blank()
                             ) + 
                   ggtitle(titles[2]),
                 ncol = 2,
                 common.legend = TRUE, legend = "bottom")
  
  return(p)
}



#' Make plot of TV-SACE(r,t) lines for cohorts defined by t
#' 
#' @param cohort_plot_dat Data set of cohorts and effects made by 
#' \link{\code{prepare_cohort_graph_data}}
#' @param legend_position Position for legend. Defaults to "none"
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot2 object
#' @export
make_cohort_tvsace_plot <- function(cohort_plot_dat, 
                                    legend_position = "none",
                                    time_unit = "Time") {
  
  p <- ggplot(data = cohort_plot_dat,
              aes_string(x = "eval_t", y = "tv_sace", 
                         group = "cohort_id", color = "frac_aa")) +
    geom_smooth(se = FALSE, alpha = 0.7) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits = c(0, 1),
                          breaks = c(0, 0.5, 1),
                          labels = c("0%", "50%", "100%")) +
    guides(color = guide_colourbar(title = "Percent of \npopulation",
                                   title.hjust = 0.5,
                                   label.position = "left")) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          legend.text.align = 0.5) +
    labs(x = paste(tools::toTitleCase(time_unit), "r"), 
         y = "TV-SACE(r, t)") +
    directlabels::geom_dl(aes(label = cohort_id), 
                          method = list(directlabels::dl.trans(x = x + 0.25),
                                        "last.points", 
                                        cex = 1))
  return(p)
}



#' Make plot of RM-SACE(r,t) lines for cohorts defined by t
#' 
#' @param cohort_plot_dat Data set of cohorts and effects made by 
#' \link{\code{prepare_cohort_graph_data}}
#' @param legend_position Position for legend. Defaults to "none"
#' @param time_unit Label for time units. Defaults to "Time"
#' @return ggplot2 object
#' @export
make_cohort_rmsace_plot <- function(cohort_plot_dat, 
                                    legend_position = "none",
                                    time_unit = "Time") {
  
  p <- ggplot(data = cohort_plot_dat,
              aes_string(x = "eval_t", y = "rm_sace", 
                         group = "cohort_id", color = "frac_aa")) +
    geom_smooth(se = FALSE, alpha = 0.7) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits = c(0, 1),
                          breaks = c(0, 0.5, 1),
                          labels = c("0%", "50%", "100%")) +
    guides(color = guide_colourbar(title = "Percent of \npopulation",
                                   title.hjust = 0.5,
                                   label.position = "left")) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          legend.text.align = 0.5) +
    labs(x = paste(tools::toTitleCase(time_unit), "r"), 
         y = "RM-SACE(r, t)") +
    directlabels::geom_dl(aes(label = cohort_id), 
                          method = list(directlabels::dl.trans(x = x + 0.25),
                                        "last.points", 
                                        cex = 1))
  return(p)
}



#' Function to take series of posterior predictive draws and plot mean proportions
#' of each principal state
#' 
#' @param pp Array of posterior predictive draws
#' @param maxt Maximum t for the plot
#' @param length_out Number of time points at which to evaluate curves (more = 
#' smoother)
#' @param color_vals Colors for principal states (order: DD, TK, CK, AA)
#' @return ggplot of state composition over time
#' @export
make_state_composition_plot <- function(pp, maxt = 90, length_out = 10,
                                        color_vals = c("#0B353B", "#800026",
                                                       "#ECB433", "#2E655D"),
                                        time_unit = "Time") {
  
  # Sequence of points to evaluate
  xt <- seq(0, maxt, length.out = length_out)
  
  # Get state matrix
  R <- dim(pp)[3]
  f_dd <- f_ck <- f_tk <- f_aa <- xt * 0 
  for (r in 1:R) {
    pstates <- v_make_pstates(eval_t = xt, pp = pp[, , r])
    f_aa <- f_aa + colMeans(pstates == "AA") / R
    f_ck <- f_ck + colMeans(pstates == "TS") / R
    f_tk <- f_tk + colMeans(pstates == "TK") / R
    f_dd <- f_dd + colMeans(pstates == "DD") / R  
  }
  
  # Data frame containing proportion in each state at each t
  state_names <- c("AA", "CK", "TK", "DD")
  comp_dat <- data.frame(Time = rep(xt, times = length(state_names)),
                         State = rep(state_names, each = length_out))
  comp_dat$Proportion <- c(f_aa, f_ck, f_tk, f_dd)
  comp_dat$State <- factor(comp_dat$State, levels = rev(state_names), 
                           ordered = TRUE)
  
  # Make the plot
  cplot <- ggplot(comp_dat, 
                  aes_string(x = "Time", y = "Proportion", fill = "State")) + 
    xlab(tools::toTitleCase(time_unit)) + 
    ylab("Proportion in state") +
    geom_area(color = "black") +
    scale_fill_manual(values = color_vals) +
    scale_x_continuous(breaks = seq(0, maxt, by = 30)) +
    ggtitle("Principal state composition") + 
    theme_minimal() + 
    theme(legend.position = "bottom")
  
  # Return
  return(cplot)
}


#TODO(LCOMM): add ROCs
