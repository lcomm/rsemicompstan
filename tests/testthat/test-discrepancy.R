library("rsemicompstan")

context("Discrepancy measures")
test_that("T_AA of a model fit to simulated data is not extreme", {
  
  # set.seed(123)
  # obs_dat <- simulate_from_param(n = 5000,
  #                                params = return_dgp_parameters(scenario = 4),
  #                                censor = TRUE)
  # 
  # obs_fit <- run_scr_replicate(n = 2000, seed = 55, scenario = 2, iter = 4000, 
  #                          chains = 1, init = "0", warmup = 3500,
  #                          # shared_beta = 1,
  #                          control = list(adapt_delta = 0.85))
  # obs_fit <- fit$stan_fit
  # ds <- prepare_discrep_plot_dat(stan_res = obs_fit, seed = 2254, eval_t = 75,
  #                          subsamp = 100)
  # calculate_frac_aa_disc(eval_t = 75, stan_res = fit$stan_fit, 
  #                        cens_times = )
  #                          
  # obs <- scr_gamma_frailty_stan()
  # rstan_options(mc.cores = 4)
  # obs_fit <- run_scr_replicate(n = 5000, seed = 10, scenario = 4, chains = 2)
  
  #   frailties <- simulate_frailty(n = 10^6, sigma = sigma)
  #   expect_equal(round(var(frailties), 1), sigma)  
  # }
  
})
