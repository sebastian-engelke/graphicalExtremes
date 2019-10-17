#' Simulates HR extremal functions
#'
#' Simulate the Huessler-Reiss extremal functions
#'
#' @param no.simu Integer. Number of simulations, by default equal to 1
#' @param idx Integer.
simu_px_HR <- function(no.simu=1, idx, trend, chol.mat) {
  stopifnot(length(idx)==1)
  d <- nrow(chol.mat)
  proc <- t(chol.mat)%*%matrix(rnorm(d*no.simu), ncol=no.simu) - trend
  proc <- exp(t(proc) - proc[idx,])
  return(proc)
}
