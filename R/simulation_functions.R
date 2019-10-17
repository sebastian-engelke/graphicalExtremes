#' Simulate samples of multivariate Pareto distribution
#'
#' Simulates exact samples of multivariate Pareto distributions
#' @param model String. The parametric model type. Is one of:
#' \itemize{
#' \item \code{HR}
#' \item \code{logistic}
#' \item \code{neglogistic}
#' \item \code{dirichlet}.
#' }
#' @param d Positive integer. Dimension of the multivariate Pareto
#' distribution.
#' @param no.simu Positive integer. Number of simulations
#' (by default equal to 1).
#' @param par Is the respective parameter for the given \code{model}.
#' Is one of:
#' \itemize{
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, if \code{model = logistic}
#' \item \eqn{\theta > 0}, if \code{model = neglogistic}
#' \item \eqn{\alpha}, numeric vector of size \code{d},
#' if \code{model = dirichlet}
#' \item \eqn{\Gamma}, numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram, if \code{model = HR}.
#' }
#' @return List. The list is made of:
#' \itemize{
#' \item \code{res} Numeric matrix of size \eqn{no.simu \times d}{no.simu x d}.
#' The simulated multivariate Pareto data.
#' \item \code{counter} Positive integer. The number of times needed to sweep
#' over the \code{d} variables to simulate \code{no.simul} multivariate
#' observations.
#' }
rmpareto <- function(model, d, no.simu=1, par) {

  stopifnot((d==round(d)) & (d>=1))
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))
  stopifnot(model %in% c("HR", "logistic", "neglogistic", "dirichlet"))

  if (model=="HR") {
    stopifnot(is.matrix(par))
    Gamma = par
    stopifnot(nrow(Gamma) == d & ncol(Gamma) == d)
    cov.mat <- sapply(1:d, function(i) sapply(1:d, function(j)
      (Gamma[i,1] + Gamma[j,1] - Gamma[i,j])/2))
    cov.mat <- cov.mat + 1e-3 ##add constant random effect to avoid numerical problems
    chol.mat <- chol(cov.mat)
  } else if (model=="logistic") {
    stopifnot(length(par) == 1 & 1e-12 < par & par < 1 - 1e-12)
    theta = par
  } else if (model=="neglogistic") {
    stopifnot(par > 1e-12)
    theta = par
  } else if (model=="dirichlet") {
    alpha = par
    stopifnot(length(alpha) == d)
    stopifnot(all(alpha>1e-12))
  }

  counter <- 0
  res <- numeric(0)
  n.total <- 0
  while (n.total < no.simu){
    counter <- counter + 1
    shift <- sample(1:d, no.simu, replace=TRUE)
    for(k in 1:d){
      if (model == "HR") {
        trend <- sapply(1:d, function(j) Gamma[j,k]/2)
      }
      n.k <- sum(shift==k)

      if(n.k>0){
        proc <- switch(model,
                       "HR"           = simu_px_HR(no.simu=n.k, idx=k, trend=trend, chol.mat=chol.mat),
                       "logistic"     = simu_px_logistic(no.simu=n.k, idx=k, N=d, theta=theta),
                       "neglogistic"  = simu_px_neglogistic(no.simu=n.k, idx=k, N=d, theta=theta),
                       "dirichlet"    = simu_px_dirichlet(no.simu=n.k, idx=k, N=d, alpha=alpha)
        )
        stopifnot(dim(proc)==c(n.k, d))
        proc <- proc/rowSums(proc) / (1-runif(nrow(proc)))
        idx.sim <- which(apply(proc,1,max) > 1)
        res <- rbind(res, proc[idx.sim,])
        n.total <- nrow(res)
      }
    }
  }
  return(list(res=res[sample(1:nrow(res), no.simu, replace=FALSE),], counter=counter))
}
