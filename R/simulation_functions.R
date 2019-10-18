#' Simulate samples of multivariate Pareto distribution
#'
#' Simulates exact samples of multivariate Pareto distributions
#' @param no.simu Positive integer. Number of simulations.
#' @param model String. The parametric model type. Is one of:
#' \itemize{
#' \item \code{HR}
#' \item \code{logistic}
#' \item \code{neglogistic}
#' \item \code{dirichlet}.
#' }
#' @param d Positive integer. Dimension of the multivariate Pareto
#' distribution.
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
#'
#' ## !!! add examples (define params and call function)
#' }
rmpareto <- function(no.simu, model, d, par) {

  stopifnot((d==round(d)) & (d>=1))
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))
  stopifnot(model %in% c("HR", "logistic", "neglogistic", "dirichlet"))

  if (model=="HR") {
    stopifnot(is.matrix(par))
    Gamma = par
    stopifnot(nrow(Gamma) == d & ncol(Gamma) == d)
    cov.mat <- Gamma2Sigma(Gamma, k=1, full=FALSE)
    chol.mat <- matrix(0,d,d)
    chol.mat[-1,-1] <- chol(cov.mat)
    # !!! add error if cannot chol()
    # !!! put here the trend and save it as a matrix (each row is one variable)
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
  } else if (model == "dirichlet_mix"){ # !!!
    # !!!
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
                       "logistic"     = simu_px_logistic(no.simu=n.k, idx=k, d=d, theta=theta),
                       "neglogistic"  = simu_px_neglogistic(no.simu=n.k, idx=k, d=d, theta=theta),
                       "dirichlet"    = simu_px_dirichlet(no.simu=n.k, idx=k, d=d, alpha=alpha),
                       "dirichlet_mix" = simu_px_dirichlet_mix(no.simu = n.k,
                                                               idx = k, d = d,
                                                               weights=..., alpha=..., norm.alpha=...)
        )
        stopifnot(dim(proc)==c(n.k, d))
        proc <- proc/rowSums(proc) / (1-runif(nrow(proc)))
        idx.sim <- which(apply(proc,1,max) > 1)
        res <- rbind(res, proc[idx.sim,])
        n.total <- nrow(res)
      }
    }
  }
  return(list(res=res[sample(1:nrow(res), no.simu, replace=FALSE),],
              counter=counter))
}
