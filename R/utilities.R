#' Return utility DGP
#' 
#' @param u_scenario Utility scenario
#' 0 = 0 utility
#' 1 = full utility
#' 2 = flat over (0, 1)
#' 3 = (p = 0.2; prior sample size = 12)
#' 4 = (p = 0.7; prior sample size = 12)
# hist(rbeta(n = 10000, p * pss, (1 - p) * pss), xlim = c(0,1))
#' @param P Number of columns of the design matrix
#' @return Named list of coefficients
#' @export
# return_utility_dgp <- function(u_scenario) {
#   switch("0" = rbeta(0, ),
#          )
# }