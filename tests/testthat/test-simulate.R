library("rsemicompstan")

# Verify whether all elements are equal
all_elements_equal <- function(x, tol = .Machine$double.eps^0.5) {
  if (length(x) == 1) return(TRUE)
  x <- range(x) / mean(x)
  isTRUE(all.equal(x[1], x[2], tolerance = tol))
}


context("Data simulation")

test_that("Simulated gamma frailties have variance sigma", {
  
  set.seed(123)
  for (sigma in c(0.4, 2, 5)) {
    frailties <- simulate_frailty(n = 10^6, sigma = sigma)
    expect_equal(round(var(frailties), 1), sigma)  
  }
  
})

test_that("Frailty coefficients match descriptions", {
  
  # Type 0: All-zero case
  expect_equal(give_frailty_coefs(fctype = 0), rep(0, 6))
  
  # Type 1: All equal, and not zero
  x <- give_frailty_coefs(fctype = 1)
  expect_gt(min(abs(x)), 0)
  expect_true(all_elements_equal(x))
  
  # Type 2: same within tx and varies across tx
  x <- give_frailty_coefs(fctype = 2)
  expect_true(all_elements_equal(x[1:3]))
  expect_true(all_elements_equal(x[4:6]))
  expect_false(isTRUE(all.equal(x[1:3], x[4:6])))
  
  # Type 3: same within transition type and varies across type
  x <- give_frailty_coefs(fctype = 3)
  expect_equal(x[1], x[4])
  expect_equal(x[2], x[5])
  expect_equal(x[3], x[6])
  expect_false(all_elements_equal(x))
  
  # Type 4: all are different
  x <- give_frailty_coefs(fctype = 4)
  expect_gt(min(abs(diff(sort(x)))), .Machine$double.eps^0.5)
  
})
  
test_that("Coefficient output has right form", {
  
  Xmats <- replicate(6, matrix(1:40, nrow = 10, ncol = 4), simplify = FALSE)
  coefs <- give_coefs(Xmats, fctype = 1, exp = TRUE)
  
  expect_named(coefs, expected = c("control", "treated"))
  expect_named(coefs[[1]], expected = c("rcoefs", "fcoefs", "shapes"))
  expect_named(coefs[[2]], expected = c("rcoefs", "fcoefs", "shapes"))
  expect_equal(coefs[[1]][["shapes"]], c(1, 1, 1))
  expect_equal(coefs[[2]][["shapes"]], c(1, 1, 1))
  
  coefs <- give_coefs(Xmats, fctype = 1, exp = FALSE)
  expect_gt(max(abs(coefs[[1]][["shapes"]] - 1)), .Machine$double.eps^0.5)
  
})

test_that("Simulated covariates have right form", {
  
  cov_dat <- simulate_covariate_data(n = 25, p = 5)
  expect_is(cov_dat, "data.frame")
  expect_true(NROW(cov_dat) == 25)
  expect_true(NCOL(cov_dat) == 5)
  expect_false(anyNA(cov_dat))
  expect_named(cov_dat, expected = paste0("X",1:5))
  
})

test_that("Design matrices have right form", {
  
  # When no data provided
  Xmats <- make_Xmats(n = 20, formulas = NULL, covariate_data = NULL)
  expect_is(Xmats, "list")
  expect_length(Xmats, 6)
  expect_is(Xmats[[1]], "matrix")
  expect_false(anyNA(Xmats))
  
  # When data already provided
  cov_dat <- simulate_covariate_data(n = 45, p = 5)
  Xmats <- make_Xmats(n = 45, formulas = NULL, covariate_data = cov_dat)
  expect_error(make_Xmats(n = 20, formulas = NULL, covariate_data = cov_dat))
  expect_is(Xmats, "list")
  expect_length(Xmats, 6)
  expect_is(Xmats[[1]], "matrix")
  expect_false(anyNA(Xmats))
  expect_true(NROW(Xmats[[1]]) == 45)
  expect_true(NCOL(Xmats[[1]]) == 5 + 1) # +1 since includes intercept
  
})


test_that("Potential outcomes output has right form", {
  
  n <- 15
  Xmats <- make_Xmats(n)
  coefs <- give_coefs(Xmats)
  frailties <- simulate_frailties(n)
  po_dat <- simulate_poutcomes(Xmats, coefs, frailties)
  expect_true(is.data.frame(po_dat))
  expect_true(NROW(po_dat) == n & NCOL(po_dat) == 8)
  vars <- c("R", "D", "deltaR", "deltaD")
  expect_named(po_dat, c(paste0(vars, "0"), paste0(vars, "1")))
  expect_false(anyNA(po_dat))
  expect_true(all(po_dat$R0 <= po_dat$D0))
  expect_true(all(po_dat$R1 <= po_dat$D1))
  expect_true(all(po_dat$deltaR0 %in% c(0, 1)))
  expect_true(all(po_dat$deltaR1 %in% c(0, 1)))
  expect_true(all(po_dat$deltaD0 %in% c(0, 1)))
  expect_true(all(po_dat$deltaD1 %in% c(0, 1)))
  
})


