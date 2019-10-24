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
#' \item \eqn{\Gamma}, numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram, if \code{model = HR}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, if \code{model = logistic}.
#' \item \eqn{\theta > 0}, if \code{model = neglogistic}.
#' \item \eqn{\alpha > 0}, numeric vector of size \code{d},
#' if \code{model = dirichlet}.
#' }
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{res} Numeric matrix of size \eqn{n \times d}{n x d}.
#' The simulated multivariate Pareto data.
#' \item \code{counter} Positive integer. The number of times needed to sweep
#' over the \code{d} variables to simulate \code{n} multivariate
#' observations.
#' }
#' ## !!! add examples (define params and call function)
#'
rmpareto <- function(n,
                     model = c("HR", "logistic", "neglogistic", "dirichlet")[1],
                     d, par) {

  # methods
  model_nms <- c("HR", "logistic", "neglogistic", "dirichlet")

  # check arguments ####
  if (d != round(d) | d < 1){
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1){
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)){
    stop(paste("The model must be one of", model_nms))
  }

  if (model == "HR") {
    if (!is.matrix(par)){
      stop("The argument par must be a matrix, when model = HR.")
    }

    if (NROW(par) != d | NCOL(par) != d){
      stop("The argument par must be a d x d matrix, when model = HR.")
    }

  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12){
      stop(paste("The argument par must be scalar between 1e-12 and 1 - 1e-12,",
                 "when model = logistic."))
    }

  } else if (model == "neglogistic") {
    if (par <= 1e-12){
      stop(paste("The argument par must be scalar greater than 1e-12,",
                 "when model = neglogistic."))
    }

  } else if (model == "dirichlet") {
    if (length(par) != d){
      stop(paste("The argument par must be a vector with d elements,",
                 "when model = dirichlet."))
    }

    if (any(par <= 1e-12)){
      stop(paste("The elements of par must be greater than 1e-12,",
                 "when model = dirichlet."))
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
      stop(paste("The covariance matrix associated to Gamma",
                 "cannot be factorized with Cholesky."))
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

        if (any(dim(proc) != c(n.k, d))) {
          stop("The generated sample has wrong size.")
        }

        proc <- proc / rowSums(proc) / (1 - runif(NROW(proc)))
        idx.sim <- which(apply(proc, 1, max) > 1)
        res <- rbind(res, proc[idx.sim, ])
        n.total <- NROW(res)

      }
    }
  }

  return(list(res = res[sample(1:NROW(res), n, replace=FALSE), ],
              counter = counter))
}



#' Simulate samples of multivariate Pareto distribution from a tree
#'
#' Simulates a tree graphical model, following a multivariate Pareto
#' distribution.
#'
#' @param n Positive integer. Number of simulations.
#' @param model String. The parametric model type. Is one of:
#' \itemize{
#' \item \code{HR} (default),
#' \item \code{logistic},
#' \item \code{dirichlet}.
#' }
#' @param tree igraph object. An igraph object representing a tree, i.e., an
#' undirected graph that is connected and has no cycles.
#' @param par Is the respective parameter for the given \code{model}.
#' Is one of:
#' \itemize{
#' \item \eqn{\Gamma}, numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram, if \code{model = HR}. Alternatively, can be a vector of
#' length \eqn{d - 1}, containing the entries of the variogram corresponding
#' to the edges of the given \code{tree}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, if \code{model = logistic}.
#' \item a matrix of size \eqn{(d - 1) \times 2}{(d - 1) x 2} containing
#' the relative \eqn{\alpha > 0} coefficients, if \code{model = dirichlet}.
#' }
#'
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{res} Numeric matrix of size \eqn{n \times d}{n x d}.
#' The simulated multivariate Pareto data.
#' \item \code{counter} Positive integer. The number of times needed to sweep
#' over the \code{d} variables to simulate \code{n} multivariate
#' observations.
#' }
#' ## !!! add examples (define params and call function)
#'
rmpareto_tree <- function(n, model = c("HR", "logistic", "dirichlet")[1],
                          tree, par) {

  # methods
  model_nms <- c("HR", "logistic", "dirichlet")

  # graph theory objects ####
  # check if it is directed
  if (igraph::is_directed(tree)){
    warning("The given tree is directed. Converted to undirected.")
    tree <- igraph::as.undirected(tree)
  }

  # set graph theory objects
  adj =  as.matrix(igraph::as_adj(tree))
  d <- NROW(adj)
  e <- igraph::ecount(tree)
  ends.mat = igraph::ends(tree, igraph::E(tree))

  # check if it is tree
  is_connected <- igraph::is_connected(tree)
  is_tree <- is_connected & (e == d - 1)

  if (!is_tree){
    stop("The given graph is not a tree.")
  }


  # check arguments ####
  if (d != round(d) | d < 1){
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1){
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)){
    stop(paste("The model must be one of", model_nms))
  }

  if (model == "HR") {
    if (!is.matrix(par)){
      if (length(par) != d){
        stop(paste("The argument par must be a d x d matrix,",
                   "or a vector with d elements, when model = HR."))
      }
    } else {
      if (NROW(par) != d | NCOL(par) != d){
        stop(paste("The argument par must be a d x d matrix,",
                   "or a vector with d elements, when model = HR."))
      }
    }
  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12){
      stop(paste("The argument par must be scalar between 1e-12 and 1 - 1e-12,",
                 "when model = logistic."))
    }

  } else if (model == "dirichlet") {
    if (NROW(par) != d-1 | NCOL(par) != 2){
      stop(paste("The argument par must be a (d-1) x 2 ,",
                 "when model = dirichlet."))
    }
    if (any(par <= 1e-12)){
      stop(paste("The elements of par must be greater than 1e-12,",
                 "when model = dirichlet."))
    }
  }

  # prepare arguments ####
  if (model == "HR"){
    if (is.matrix(par)){
      par.vec <- par[ends.mat]
    } else {
      par.vec <- par
    }
  } else if(model == "logistic"){
    theta <- par
  } else if(model == "dirichlet"){
    alpha.mat <- par
  }

  # function body ####
  ## Define a matrix A[[k]] choosing the paths from k to other vertices
  idx.e <- matrix(0, nrow=d, ncol=d)
  idx.e[ends.mat] = 1:e
  idx.e = idx.e + t(idx.e)

  # e.start[[k]][h] gives the index (1 or 2) of the starting node in the h edge
  # in the tree rooted at k
  A <- e.start <- e.end <- list()
  for (k in 1:d) {
    A[[k]] <- matrix(0, nrow=d, ncol=e)
    e.start[[k]] = e.end[[k]] = numeric(e)
    short.paths <- igraph::shortest_paths(tree, from = k, to=1:d)
    for(h in 1:d){
      path = short.paths$vpath[[h]]
      idx.tmp = idx.e[cbind(path[-length(path)], path[-1])]
      A[[k]][h,idx.tmp] <- 1
      e.start[[k]][idx.tmp] =
        apply(ends.mat[idx.tmp,] ==
                matrix(path[-length(path)], nrow = length(idx.tmp), ncol=2),
              MARGIN=1, FUN = function(x) which(x==TRUE))
      e.end[[k]][idx.tmp] =
        apply(ends.mat[idx.tmp,] ==
                matrix(path[-1], nrow = length(idx.tmp), ncol=2),
              MARGIN=1, FUN = function(x) which(x==TRUE))
    }
  }

  counter <- 0
  res <- numeric(0)
  n.total <- 0
  while (n.total < n) {
    counter <- counter + 1
    shift <- sample(1:d, n, replace=TRUE)
    for(k in 1:d){
      n.k <- sum(shift==k)
      if(n.k>0){
        proc <-
          switch(model,
                 "HR" =
                   simu_px_tree_HR(n=n.k, G.vec=par.vec, A = A[[k]]),
                 "logistic" =
                   simu_px_tree_logistic(n=n.k, idx=k, nb.edges=e,
                                         theta=theta, A=A),
                 "dirichlet" =
                   simu_px_tree_dirichlet(n = n.k,
                                          alpha.start =
                                            alpha.mat[cbind(1:e, e.start[[k]])],
                                          alpha.end =
                                            alpha.mat[cbind(1:e, e.end[[k]])],
                                          A=A[[k]])
        )

        if (any(dim(proc) != c(n.k, d))) {
          stop("The generated sample has wrong size.")
        }

        proc <- proc/rowSums(proc) / (1-runif(NROW(proc)))
        idx.sim <- which(apply(proc,1,max) > 1)
        res <- rbind(res, proc[idx.sim,])
        n.total <- NROW(res)
      }
    }
  }

  return(list(res=res[sample(1:NROW(res), n, replace=FALSE),], counter=counter))
}



#' Simulate samples of max-stable process
#'
#' Simulates exact samples of max-stable process
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
#' \item \eqn{\Gamma}, numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram, if \code{model = HR}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, if \code{model = logistic}.
#' \item \eqn{\theta > 0}, if \code{model = neglogistic}.
#' \item \eqn{\alpha > 0}, numeric vector of size \code{d},
#' if \code{model = dirichlet}.
#' }
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{res} Numeric matrix of size \eqn{n \times d}{n x d}.
#' The simulated multivariate Pareto data.
#' \item \code{counter} Positive integer. The number of times needed to sweep
#' over the \code{d} variables to simulate \code{n} multivariate
#' observations.
#' }
#' ## !!! add examples (define params and call function)
#'
rmstable <- function(n,
                     model = c("HR", "logistic", "neglogistic", "dirichlet")[1],
                     d, par) {

  # methods
  model_nms <- c("HR", "logistic", "neglogistic", "dirichlet")

  # check arguments ####
  if (d != round(d) | d < 1){
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1){
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)){
    stop(paste("The model must be one of", model_nms))
  }

  if (model == "HR") {
    if (!is.matrix(par)){
      stop("The argument par must be a matrix, when model = HR.")
    }

    if (NROW(par) != d | NCOL(par) != d){
      stop("The argument par must be a d x d matrix, when model = HR.")
    }

  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12){
      stop(paste("The argument par must be scalar between 1e-12 and 1 - 1e-12,",
                 "when model = logistic."))
    }

  } else if (model == "neglogistic") {
    if (par <= 1e-12){
      stop(paste("The argument par must be scalar greater than 1e-12,",
                 "when model = neglogistic."))
    }

  } else if (model == "dirichlet") {
    if (length(par) != d){
      stop(paste("The argument par must be a vector with d elements,",
                 "when model = dirichlet."))
    }

    if (any(par <= 1e-12)){
      stop(paste("The elements of par must be greater than 1e-12,",
                 "when model = dirichlet."))
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
      stop(paste("The covariance matrix associated to Gamma",
                 "cannot be factorized with Cholesky."))
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

  }

  # function body ####
  counter <- rep(0, times=n)
  res <- matrix(0, nrow=n, ncol=d)
  for (k in 1:d) {
    poisson <- rexp(n)

    while (any(1/poisson > res[,k])) {
      ind <- (1/poisson > res[,k])
      n.ind <- sum(ind)
      idx <- (1:n)[ind]
      counter[ind] <- counter[ind] + 1
      proc <-
        switch(model,
               "HR" =
                 simu_px_HR(n = n.ind, idx = k, d = d, trend = trend[k, ],
                            chol_mat = chol_mat),
               "logistic" =
                 simu_px_logistic(n = n.ind, idx = k, d = d, theta = theta),
               "neglogistic" =
                 simu_px_neglogistic(n = n.ind, idx = k, d = d, theta = theta),
               "dirichlet" =
                 simu_px_dirichlet(n = n.ind, idx = k, d = d, alpha = alpha))

      if (any(dim(proc) != c(n.ind, d))) {
        stop("The generated sample has wrong size.")
      }

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

  return(list(res = res[sample(1:NROW(res), n, replace=FALSE), ],
              counter = counter))
}



#' Simulate samples of max-stable process from a tree
#'
#' Simulates a tree graphical model, following a max-stable process.
#'
#' @inheritParams rmpareto_tree
#'
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{res} Numeric matrix of size \eqn{n \times d}{n x d}.
#' The simulated multivariate Pareto data.
#' \item \code{counter} Positive integer. The number of times needed to sweep
#' over the \code{d} variables to simulate \code{n} multivariate
#' observations.
#' }
#' ## !!! add examples (define params and call function)
#'
rmstable_tree <- function(n, model = c("HR", "logistic", "dirichlet")[1],
                          tree, par) {

  # methods
  model_nms <- c("HR", "logistic", "dirichlet")

  # graph theory objects ####
  # check if it is directed
  if (igraph::is_directed(tree)){
    warning("The given tree is directed. Converted to undirected.")
    tree <- igraph::as.undirected(tree)
  }

  # set graph theory objects
  adj =  as.matrix(igraph::as_adj(tree))
  d <- NROW(adj)
  e <- igraph::ecount(tree)
  ends.mat = igraph::ends(tree, igraph::E(tree))

  # check if it is tree
  is_connected <- igraph::is_connected(tree)
  is_tree <- is_connected & (e == d - 1)

  if (!is_tree){
    stop("The given graph is not a tree.")
  }


  # check arguments ####
  if (d != round(d) | d < 1){
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1){
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)){
    stop(paste("The model must be one of", model_nms))
  }

  if (model == "HR") {
    if (!is.matrix(par)){
      if (length(par) != d){
        stop(paste("The argument par must be a d x d matrix,",
                   "or a vector with d elements, when model = HR."))
      }
    } else {
      if (NROW(par) != d | NCOL(par) != d){
        stop(paste("The argument par must be a d x d matrix,",
                   "or a vector with d elements, when model = HR."))
      }
    }
  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12){
      stop(paste("The argument par must be scalar between 1e-12 and 1 - 1e-12,",
                 "when model = logistic."))
    }

  } else if (model == "dirichlet") {
    if (NROW(par) != d-1 | NCOL(par) != 2){
      stop(paste("The argument par must be a (d-1) x 2 ,",
                 "when model = dirichlet."))
    }
    if (any(par <= 1e-12)){
      stop(paste("The elements of par must be greater than 1e-12,",
                 "when model = dirichlet."))
    }
  }

  # prepare arguments ####
  if (model == "HR"){
    if (is.matrix(par)){
      par.vec <- par[ends.mat]
    } else {
      par.vec <- par
    }
  } else if(model == "logistic"){
    theta <- par
  } else if(model == "dirichlet"){
    alpha.mat <- par
  }

  # function body ####
  ## Define a matrix A[[k]] choosing the paths from k to other vertices
  idx.e <- matrix(0, nrow=d, ncol=d)
  idx.e[ends.mat] = 1:e
  idx.e = idx.e + t(idx.e)

  # e.start[[k]][h] gives the index (1 or 2) of the starting node in the h edge
  # in the tree rooted at k
  A <- e.start <- e.end <- list()
  for (k in 1:d) {
    A[[k]] <- matrix(0, nrow=d, ncol=e)
    e.start[[k]] = e.end[[k]] = numeric(e)
    short.paths <- igraph::shortest_paths(tree, from = k, to=1:d)
    for(h in 1:d){
      path = short.paths$vpath[[h]]
      idx.tmp = idx.e[cbind(path[-length(path)], path[-1])]
      A[[k]][h,idx.tmp] <- 1
      e.start[[k]][idx.tmp] =
        apply(ends.mat[idx.tmp,] ==
                matrix(path[-length(path)], nrow = length(idx.tmp), ncol=2),
              MARGIN=1, FUN = function(x) which(x==TRUE))
      e.end[[k]][idx.tmp] =
        apply(ends.mat[idx.tmp,] ==
                matrix(path[-1], nrow = length(idx.tmp), ncol=2),
              MARGIN=1, FUN = function(x) which(x==TRUE))
    }
  }


  counter <- rep(0, times=n)
  res <- matrix(0, nrow=n, ncol=d)
  for (k in 1:d) {
    poisson <- rexp(n)

    while (any(1/poisson > res[,k])) {
      ind <- (1/poisson > res[,k])
      n.ind <- sum(ind)
      idx <- (1:n)[ind]
      counter[ind] <- counter[ind] + 1
      proc <-
        switch(model,
               "HR" =
                 simu_px_tree_HR(n=n.ind, G.vec=par.vec, A = A[[k]]),
               "logistic" =
                 simu_px_tree_logistic(n=n.ind, idx=k, nb.edges=e,
                                       theta=theta, A=A),
               "dirichlet" =
                 simu_px_tree_dirichlet(n=n.ind,
                                        alpha.start =
                                          alpha.mat[cbind(1:e, e.start[[k]])],
                                        alpha.end =
                                          alpha.mat[cbind(1:e, e.end[[k]])],
                                        A=A[[k]])
      )

      if (any(dim(proc) != c(n.ind, d))) {
        stop("The generated sample has wrong size.")
      }

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

  return(list(res=res[sample(1:NROW(res), n, replace=FALSE),], counter=counter))
}
