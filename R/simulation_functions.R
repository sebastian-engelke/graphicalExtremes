#' Simulate samples of multivariate Pareto distribution
#'
#' Simulates exact samples of multivariate Pareto distributions
#' @param n Positive integer. Number of simulations.
#' @param model String. The parametric model type. Is one of:
#' \itemize{
#' \item \code{HR} (default),
#' \item \code{logistic},
#' \item \code{neglogistic},
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
#' \item list, made of ???, if \code{model = dirichlet_mix}
#' }
#' @return List. The list is made of:
#' \itemize{
#' \item \code{res} Numeric matrix of size \eqn{n \times d}{n x d}.
#' The simulated multivariate Pareto data.
#' \item \code{counter} Positive integer. The number of times needed to sweep
#' over the \code{d} variables to simulate \code{n} multivariate
#' observations.
#'
#' ## !!! add examples (define params and call function)
#' }
rmpareto <- function(n,
                     model = c("HR", "logistic", "neglogistic", "dirichlet")[1],
                     d, par) {

  # methods
  methods_nms <- c("HR", "logistic", "neglogistic", "dirichlet")

  # check arguments ####
  if (d != round(d) | d < 1){
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) & n < 1){
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% methods_nms)){
    stop(paste("The model must be one of", methods_nms))
  }

  if (model == "HR") {
    if (!is.matrix(par)){
      stop("The argument par must be a matrix, when model = HR.")
    }

    if (nrow(par) != d | ncol(par) != d){
      stop("The argument par must be a d x d matrix, when model = HR.")
    }

  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12){
      stop("The argument par must be scalar between 1e-12 and 1 - 1e-12,
           when model = logistic.")
    }

  } else if (model == "neglogistic") {
    if (par <= 1e-12){
      stop("The argument par must be scalar greater than 1e-12,
           when model = neglogistic.")
    }

  } else if (model == "dirichlet") {
    if (length(par) != d){
      stop("par must be a vector with d elements,
           when model = dirichlet.")
    }

    if (any(par <= 1e-12)){
      stop("The elements of par must be greater than 1e-12,
           when model = dirichlet.")
    }

  }

  # prepare arguments ####
  if (model == "HR") {
    Gamma <- par

    # compute cholesky decomposition
    cov.mat <- Gamma2Sigma(Gamma, k = 1, full = FALSE)
    chol_mat <- matrix(0, d, d)

    result <- tryCatch({
      chol_mat[-1, -1] <- chol(cov.mat)
    },
    error = function(e) {
      stop("The covariance matrix associated to Gamma cannot be factorized
           with Cholesky.")
    })

    # compute trend (matrix where each row is one variable)
    trend <- t(sapply(1:d, function(k){
      sapply(1:d, function(j){
        Gamma[j, k] / 2
        }
      )}))

  } else if (model == "logistic") {
    theta <- par

  } else if (model == "neglogistic") {
    theta <- par

  } else if (model == "dirichlet") {
    alpha <- par

  } else if (model == "dirichlet_mix"){
    # ???
  }

  # function body ####
  counter <- 0
  res <- numeric(0)
  n.total <- 0
  while (n.total < n){
    counter <- counter + 1
    shift <- sample(1:d, n, replace = TRUE)
    for (k in 1:d){

      n.k <- sum(shift == k)

      if (n.k > 0){
        proc <-
          switch(model,
                 "HR" =
                   simu_px_HR(n = n.k, idx = k, d = d, trend = trend[k, ],
                              chol_mat = chol_mat),
                 "logistic" =
                   simu_px_logistic(n = n.k, idx = k, d = d, theta = theta),
                 "neglogistic" =
                   simu_px_neglogistic(n = n.k, idx = k, d = d, theta = theta),
                 "dirichlet" =
                   simu_px_dirichlet(n = n.k, idx = k, d = d, alpha = alpha)
        )

        if (dim(proc) != c(n.k, d)) {
          stop("The generated sample has wrong size.")
        }

        proc <- proc / rowSums(proc) / (1 - runif(nrow(proc)))
        idx.sim <- which(apply(proc, 1, max) > 1)
        res <- rbind(res, proc[idx.sim, ])
        n.total <- nrow(res)

      }
    }
  }

  return(list(res = res[sample(1:nrow(res), n, replace=FALSE), ],
              counter = counter))
}


### !!! This function simulates tree graphical models, either multivariate
#Pareto or max-stable distributions
#tree: a graph object that must be a tree
#model: either "HR" or "logistic" for HR or logistics tree model, respectively
#method: either "mpareto" or "maxstable"
#loc, scale, shape: if method="maxstable", output is transformed to
#general GEV margins
#n: number of simulations
#Gamma: parameter matrix if model="HR"
#theta: parameter if model="logsitic"
simu_tree <- function(tree, model, method, n=1, Gamma=NULL, theta=NULL,
                      alpha.mat=NULL, loc=1, scale=1, shape=1) {
  require("igraph")
  adj =  as.matrix(as_adj(tree))
  d <- nrow(adj)
  e <- ecount(tree)
  ends.mat = ends(tree, E(tree))

  stopifnot(model %in% c("logistic", "HR", "dirichlet"))
  stopifnot((d==round(d)) & (d>=1))
  stopifnot((n==round(n)) & (n>=1))

  if (length(loc)  ==1) loc   <- rep(loc  , times=d)
  if (length(scale)==1) scale <- rep(scale, times=d)
  if (length(shape)==1) shape <- rep(shape, times=d)
  stopifnot(all(scale>1e-12))

  if (model=="logistic") {
    stopifnot(1e-12 < theta & theta < 1 - 1e-12)
  } else if (model=="HR") {
    par.vec = Gamma[ends.mat]
  } else if (model=="dirichlet") {
    stopifnot(nrow(alpha.mat) == d-1 & ncol(alpha.mat) == 2)
  }

  ## Define a matrix A[[k]] choosing the paths from k to other vertices
  idx.e <- matrix(0, nrow=d, ncol=d)
  idx.e[ends.mat] = 1:e
  idx.e = idx.e + t(idx.e)

  A <- e.start <- e.end <- list() #e.start[[k]][h] gives the index (1 or 2) of the starting node in the h edge in the tree rooted at k

  for (k in 1:d) {
    A[[k]] <- matrix(0, nrow=d, ncol=e)
    e.start[[k]] = e.end[[k]] = numeric(e)
    short.paths <- shortest_paths(tree, from = k, to=1:d)
    for(h in 1:d){
      path = short.paths$vpath[[h]]
      idx.tmp = idx.e[cbind(path[-length(path)], path[-1])]
      A[[k]][h,idx.tmp] <- 1
      e.start[[k]][idx.tmp] = apply(ends.mat[idx.tmp,] == matrix(path[-length(path)], nrow = length(idx.tmp), ncol=2), MARGIN=1, FUN = function(x) which(x==TRUE)) #path[-length(path)]
      e.end[[k]][idx.tmp] = apply(ends.mat[idx.tmp,] == matrix(path[-1], nrow = length(idx.tmp), ncol=2), MARGIN=1, FUN = function(x) which(x==TRUE))  #path[-1]
    }
  }

  if(method=="mpareto")
  {
    counter <- 0
    res <- numeric(0)
    n.total <- 0
    while (n.total < n) {
      counter <- counter + 1
      shift <- sample(1:d, n, replace=TRUE)
      for(k in 1:d){
        n.k <- sum(shift==k)
        if(n.k>0){
          proc <- switch(model,
                         "HR" = simu_px_tree_HR(n=n.k, G.vec=par.vec, A = A[[k]]),
                         "logistic"     = simu_px_tree_logistic(n=n.k, idx=k, nb.edges=e, theta=theta, A=A),
                         "dirichlet"     = simu_px_tree_dirichlet(n=n.k, alpha.start = alpha.mat[cbind(1:e, e.start[[k]])],
                                                                  alpha.end = alpha.mat[cbind(1:e, e.end[[k]])], A=A[[k]])
          )
          stopifnot(dim(proc)==c(n.k, d))
          proc <- proc/rowSums(proc) / (1-runif(nrow(proc)))
          idx.sim <- which(apply(proc,1,max) > 1)
          res <- rbind(res, proc[idx.sim,])
          n.total <- nrow(res)
        }
      }
    }
  }else if(method=="maxstable"){
    counter <- rep(0, times=n)
    res <- matrix(0, nrow=n, ncol=d)
    for (k in 1:d) {
      poisson <- rexp(n)

      while (any(1/poisson > res[,k])) {
        ind <- (1/poisson > res[,k])
        n.ind <- sum(ind)
        idx <- (1:n)[ind]
        counter[ind] <- counter[ind] + 1
        proc <- switch(model,
                       "HR" = simu_px_tree_HR(n=n.ind, G.vec=par.vec, A = A[[k]]),
                       "logistic"     = simu_px_tree_logistic(n=n.ind, idx=k, nb.edges=e, theta=theta, A=A),
                       "dirichlet"     = simu_px_tree_dirichlet(n=n.ind, alpha.start = alpha.mat[cbind(1:e, e.start[[k]])],
                                                                alpha.end = alpha.mat[cbind(1:e, e.end[[k]])], A=A[[k]])
        )
        stopifnot(dim(proc)==c(n.ind, d))
        if (k==1) {
          ind.upd <- rep(TRUE, times=n.ind)
        } else {
          ind.upd <- sapply(1:n.ind, function(i)
            all(1/poisson[idx[i]]*proc[i,1:(k-1)] <= res[idx[i],1:(k-1)]))
        }
        if (any(ind.upd)) {
          idx.upd <- idx[ind.upd]
          res[idx.upd,] <- pmax(res[idx.upd,], 1/poisson[idx.upd]*proc[ind.upd,])
        }
        poisson[ind] <- poisson[ind] + rexp(n.ind)
      }
    }
    res <- sapply(1:d, function(i) {
      if (abs(shape[i]<1e-12)) {
        return(log(res[,i])*scale[i] + loc[i])
      } else {
        return(1/shape[i]*(res[,i]^shape[i]-1)*scale[i] + loc[i])
      }
    })

  }
  return(list(res=res[sample(1:nrow(res), n, replace=FALSE),], counter=counter))
}


