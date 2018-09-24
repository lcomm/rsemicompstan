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



#' Turn posterior predictive data into a long format amenable to graphing the 
#' causal effect
#' 
#' @param pp Posterior prediction array
#' @param max_t Maximum time for calculating effects
#' @param length_out Number of time points at which to calculate effects. More points
#' means smoother effect lines
#' @return Data frame with MCMC iteration (r), evaluation time point (eval_t),
#' time-varying survivor average causal effect (tv_sace), and restricted mean survivor
#' average causal effect (rm_sace)
#' @export
prepare_graph_data <- function(pp, max_t, length_out = 10) {
  R <- dim(pp)[3]
  xt <- seq(0, max_t, length.out = length_out)
  res <- as.data.frame(expand.grid(r = 1:R, eval_t = xt, 
                                   frac_aa = NA, tv_sace = NA, rm_sace = NA))
  
  for (r in 1:R) {
    pp_mcmc <- as.data.frame(pp[, , r])
    for (eval_t in xt) {
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



#' Make always-alive survival plot
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @return ggplot object
#' @export
make_aa_kmplot <- function(plot_dat) {
  # Calculate overall mean 
  f_aa_mean_dat <- aggregate(frac_aa ~ eval_t, data = plot_dat, FUN = mean)
  f_aa_mean_dat$r <- 1
  
  # Make plot
  p <- ggplot(data = plot_dat,
              aes_string(y = "frac_aa", x = "eval_t", group = "r")) + 
    geom_line(alpha = 0.15) + 
    ylim(0, 1) + 
    geom_line(data = f_aa_mean_dat,
              aes_string(y = "frac_aa", x = "eval_t")) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + 
    labs(x = "Time t", y = "P(Always-alive at t)")
  
  return(p)
}



#' Make TV-SACE curve plot (where r = t)
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @return ggplot object
#' @export
make_tvsace_plot <- function(plot_dat) {
  tv_mean_dat <- aggregate(tv_sace ~ eval_t, data = plot_dat, FUN = mean)
  
  p <- ggplot(plot_dat, aes_string(x = "eval_t", y = "tv_sace", group = "r",
                                   color = "frac_aa")) + 
    geom_line(alpha = 0.4) + 
    geom_line(data = tv_mean_dat, 
              aes_string(x = "eval_t", y = "tv_sace", group = NULL), 
              color = "black",
              size = 0.7) +
    scale_color_gradientn("Proportion Always-Alive",
                          colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                          limits = c(0, 1),
                          breaks = c(0, 1),
                          guide = TRUE) + 
    guides(alpha = FALSE) + 
    theme_minimal() + 
    theme(legend.position = "right") + 
    labs(x = "Time t", y = "TV-SACE(t)") + 
    ggtitle("Time-varying survivor average causal effect")
  
  return(p)
}



#' Make RM-SACE curve plot (where r = t)
#' 
#' @param plot_dat Plot data set from \code{\link{prepare_graph_data}}
#' @return ggplot object
#' @export
make_rmsace_plot <- function(plot_dat) {
  rm_mean_dat <- aggregate(rm_sace ~ eval_t, data = plot_dat, FUN = mean)
  
  p <- ggplot(plot_dat, aes_string(x = "eval_t", y = "rm_sace", group = "r",
                                   color = "frac_aa")) + 
    geom_line(alpha = 0.4) +
    geom_line(data = rm_mean_dat, 
              aes_string(x = "eval_t", y = "rm_sace", group = NULL), 
              color = "black",
              size = 0.7) +
    scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd"),
                         limits = c(0, 1),
                         breaks = c(0, 1)) + 
    guides(alpha = FALSE, color = FALSE) + 
    theme_minimal() + 
    labs(x = "Time t", y = "RM-SACE(t)") + 
    ggtitle("Restricted mean survivor average causal effect")
  
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


#TODO(LCOMM): add ROCs
