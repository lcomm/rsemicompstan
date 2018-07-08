set.seed(12345)

MODELS_HOME <- "src/stan_files"
INCLUDE_DIR <- "include"

# context("setup")
# test_that("Stan programs are available", {
#   expect_true(file.exists(MODELS_HOME))
# })

library("rstan")
Sys.unsetenv("R_TESTS")

functions <- sapply(dir(MODELS_HOME, pattern = "stan$", full.names = TRUE), 
                    function(f) {
  mc <- readLines(f)
  mc <- grep("^#include", mc, invert = TRUE, value = TRUE)
  start <- grep("^functions[[:blank:]]*\\{[[:blank:]]*$", mc)
  if (length(start) == 1) {
    end <- grep("^}[[:blank:]]*$", mc)[1]
    if (end == (start + 1L)) return(as.character(NULL))
    return(mc[(start + 1L):(end - 1L)])
  } else return(as.character(NULL))
})
names(functions) <- basename(names(functions))
functions <- lapply(functions, FUN = function(x) {
  if (length(x) == 0) { x <- " "; }; return(x); })
model_code <- paste(c("functions {", unlist(functions), "}", collapse = "\n"))
model_code <- paste(model_code, collapse = "\n")
stanc_ret <- rstan::stanc(model_code = model_code, model_name = "Stan Functions",
                   allow_undefined = TRUE)
rstan::expose_stan_functions(stanc_ret)



context("Hazard functions")
# test_that("Cumulative Weibull hazard function matches R", {
#   log_alpha = log(1.145)
#   log_kappa = log(1/100)
#   alpha <- exp(log_alpha)
#   kappa <- exp(log_kappa)
#   eval_t = 0.5
#   
#   # Cumulative hazard c_haz H(t) = - log(1 - F(t))
#   expect_equal(-pweibull(eval_t, scale = exp(-log(kappa)/alpha), shape = alpha, 
#                          lower = FALSE, log = TRUE),
#                c_haz(eval_t, log_kappa, log_alpha, lp = 0))
#   
# })

# test_that("Log-instantaneous Weibull hazard function matches R", {
#   log_alpha = log(1.145)
#   log_kappa = log(1/100)
#   alpha <- exp(log_alpha)
#   kappa <- exp(log_kappa)
#   eval_t = 0.5
#   
#   # Log-instantaneous hazard log_i_haz
#   # h(t) = f(t) / S(t)
#   expect_equal(log(dweibull(eval_t, scale = exp(-log(kappa)/alpha), 
#                             shape = alpha) /
#                      (pweibull(eval_t, scale = exp(-log(kappa)/alpha), 
#                                shape = alpha, lower = FALSE))),
#                log_i_haz(eval_t, log_kappa, log_alpha, lp = 0))
#   
#   
# })



