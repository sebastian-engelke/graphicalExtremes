#' Simulates HR extremal functions
#'
#' Simulate the Huessler-Reiss extremal functions
#'
#' @param no.simu Integer. Number of simulations, by default equal to 1
#' @param idx Integer. Related index corresponding to the variable  over which variable we use  extremal function is simulated.
#' @param trend
simu_px_HR <- function(no.simu=1, idx, trend, chol.mat, d) {
  stopifnot(length(idx)==1)
  d <- nrow(chol.mat)
  proc <- t(chol.mat)%*%matrix(rnorm(d*no.simu), ncol=no.simu) - trend
  proc <- exp(t(proc) - proc[idx,])
  return(proc)
}



### Internal: simulates logistic extremal functions
simu_px_logistic <- function(no.simu=1, idx, N, theta) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- matrix(1/gamma(1-theta)*(-log(runif(no.simu*N)))^(-theta), nrow=no.simu, ncol=N)
  res[cbind(1:no.simu,idx)] <- 1/gamma(1-theta)*rgamma(no.simu,shape=1-theta)^(-theta)
  return(res/res[cbind(1:no.simu,idx)])
}


### Internal: simulates negative logistic extremal functions
simu_px_neglogistic <- function(no.simu=1, idx, N, theta) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- matrix(rweibull(no.simu*N, shape=theta, scale=1/gamma(1+1/theta)), nrow=no.simu, ncol=N)
  res[cbind(1:no.simu,idx)] <- 1/gamma(1+1/theta)*rgamma(no.simu,shape=1+1/theta)^(1/theta)
  return(res/res[cbind(1:no.simu,idx)])
}


### Internal: simulates Dirichlet extremal functions
simu_px_dirichlet <- function(no.simu, idx, N, alpha) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  shape <- alpha
  shape[idx] <- alpha[idx] + 1
  shape.mat <- matrix(shape, nrow=N, ncol=no.simu)
  rate.mat <- matrix(alpha, nrow=N, ncol=no.simu)
  res <- t(matrix(rgamma(N*no.simu, shape=shape.mat, rate=rate.mat), nrow=N, ncol=no.simu))
  return(res/res[cbind(1:no.simu,idx)])
}


### Internal: simulates Dirichlet mixture extremal functions
simu_px_dirichlet_mix <- function(no.simu, idx, N, weights, alpha, norm.alpha) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  if (length(idx)==1) {
    k <- sample(1:length(weights), no.simu, replace=TRUE, prob=N*weights*norm.alpha[idx,])
  } else {
    k <- sapply(1:no.simu, function(i) sample(1:length(weights), 1, prob=N*weights*norm.alpha[idx[i],]))
  }
  shape.mat <- alpha[,k,drop=FALSE]
  shape.mat[cbind(idx,1:no.simu)] <- shape.mat[cbind(idx,1:no.simu)]+1
  res <- t(matrix(rgamma(N*no.simu, shape=shape.mat), nrow=N, ncol=no.simu))
  return(res/res[cbind(1:no.simu,idx)])
}


### Internal: simulates HR extremal functions on a tree
simu_px_tree_HR <- function(no.simu=1, G.vec, A) {
  res <- exp(A %*% matrix(rnorm(length(G.vec)*no.simu, mean= -G.vec/2, sd=sqrt(G.vec)), ncol=no.simu))
  return(t(res))
}

### Internal: simulates logistic extremal functions on a tree
simu_px_tree_logistic <- function(no.simu=1, idx, nb.edges, theta, A) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- exp(A[[idx]] %*% log(matrix(1/gamma(1-theta)*(-log(runif(no.simu*nb.edges)))^(-theta) /
                                             (1/gamma(1-theta)*rgamma(no.simu*nb.edges,shape=1-theta)^(-theta)), ncol=no.simu)))
  return(t(res))
}

### Internal: simulates Dirichlet extremal functions on a tree
simu_px_tree_dirichlet <- function(no.simu=1, alpha.start, alpha.end, A) {
  e = length(alpha.start)
  shape.start = matrix(alpha.start + 1, nrow=e, ncol=no.simu)
  rate.start = matrix(alpha.start, nrow=e, ncol=no.simu)
  shape.end = matrix(alpha.end, nrow=e, ncol=no.simu)
  rate.end = matrix(alpha.end, nrow=e, ncol=no.simu)
  sim.start = matrix(rgamma(e*no.simu, shape=shape.start, rate=rate.start), nrow=e, ncol=no.simu)
  sim.end = matrix(rgamma(e*no.simu, shape=shape.end, rate=rate.end), nrow=e, ncol=no.simu)
  res       <- exp(A %*% log(sim.end / sim.start))
  return(t(res))
}



