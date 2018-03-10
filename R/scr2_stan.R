#' Semicompeting risks frailty case 2
#'
#' @param tx Binary vector of length N for observed treatment assignment
#' @param x_1 N x P_1 design matrix for hazard model 1 (nonterminal)
#' @param x_2 N x P_2 design matrix for hazard model 3 (terminal, without nt)
#' @param x_3 N x P_3 design matrix for hazard model 3 (sojourn)
#' @param yr Vector of observed nonterminal event/censoring/truncation times
#' @param dyr Event indicator of observing nonterminal event (0=no, 1=yes)
#' @param yt Vector of observed terminal event/censoring times
#' @param dyt Event indicator of observing terminal event (0=no, 1=yes)
#' @param \dots Additional arguments to pass to to `rstan::sampling()`
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @export
scr2_stan <- function(tx, x_1, x_2, x_3, yr, dyr, yt, dyt, ...) {
  out <- rstan::sampling(stanmodels$scr2, 
                         data = list(N = length(tx), tx = tx, 
                                     x_1 = x_1, x2 = x_2, x_3 = x_3,
                                     P_1 = ncol(x_1), 
                                     P_2 = ncol(x_2), 
                                     P_3 = ncol(x_3),
                                     yr = yr, dyr = dyr, yt = yt, dyt = dyt),
                         ...)
  return(out)
}
