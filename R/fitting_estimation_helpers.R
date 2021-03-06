#' Compute exponent measure
#'
#' Computes the exponent measure of HR distribution.
#'
#' @param x Numeric vector with \eqn{d} positive elements
#' where the exponent measure is to be evaluated.
#' @param par Numeric vector with
#' \eqn{\frac{d(d - 1)}{2}}{d x (d - 1) / 2} elements.
#' It represents the upper triangular portion of a
#' variogram matrix \eqn{\Gamma}.
#'
#' @return Numeric. The exponent measure of the HR distribution.
#'
#' @keywords internal
V_HR <- function(x, par) {
  # helper function ####
  f1 <- function(i, x) {
    S <- Gamma2Sigma(G, k = i)
    return(1 / x[i] * mvtnorm::pmvnorm(
      upper = (log(x / x[i]) + G[, i] / 2)[-i],
      mean = rep(0, d - 1), sigma = S
    )[1])
  }

  # function body ####
  if (any(is_leq(x, 0))) {
    stop("The elements of x must be positive.")
  }

  d <- length(x)
  G <- par2Gamma(par)

  if (NROW(G) != d) {
    stop("The length of par must be d * (d - 1) / 2.")
  }

  return(sum(apply(cbind(1:d), 1, f1, x = x)))
}


#' Compute the exponent measure density of HR distribution
#'
#' Computes the exponent measure density of HR distribution.
#'
#' @param x Numeric matrix \eqn{n\times d}{n x d} or vector with \eqn{d}
#' elements.
#' @inheritParams V_HR
#'
#' @return Numeric. The censored exponent measure of the HR distribution.
#'
#' @keywords internal
logdV_HR <- function(x, par) {
  if (any(is_leq(x, 0))) {
    stop("The elements of x must be positive.")
  }

  if (is.vector(x)) {
    d <- length(x)
  }
  if (is.matrix(x)) {
    d <- ncol(x)
  }

  i <- 1
  G <- par2Gamma(par)

  if (NROW(G) != d) {
    stop("The length of par must be d * (d - 1) / 2.")
  }

  S <- Gamma2Sigma(G, k = i)
  cholS <- chol(S)
  Sm1 <- chol2inv(cholS)
  logdetS <- 2 * sum(log(diag(cholS)))
  if (is.vector(x)) {
    y <- (log(x / x[i]) + G[, i] / 2)[-i]
    logdv <- -sum(log(x)) - log(x[i]) - ((d - 1) / 2) * log(2 * pi) -
      1 / 2 * logdetS - 1 / 2 * t(y) %*% Sm1 %*% y
  }
  if (is.matrix(x)) {
    y <- (t(t(log(x / x[, i])) + G[, i] / 2))[, -i, drop = FALSE]
    logdv <- -apply(log(x), 1, sum) - log(x[, i]) - ((d - 1) / 2) * log(2 * pi) -
      1 / 2 * logdetS - 1 / 2 * fast_diag(y, Sm1)
  }
  return(logdv)
}





#' Compute censored exponent measure
#'
#' Computes the censored exponent measure density of HR distribution.
#'
#' @param x Numeric vector with \eqn{d} positive elements
#' where the censored exponent measure is to be evaluated.
#' @param K Integer vector, subset of \eqn{\{1, \dots, d\}}{{1, ..., d}}.
#' The index set that is not censored.
#' @inheritParams V_HR
#'
#' @return Numeric. The censored exponent measure of the HR distribution.
#'
#' @keywords internal
logdVK_HR <- function(x, K, par) {
  if (any(is_leq(x, 0))) {
    stop("The elements of x must be positive.")
  }

  d <- length(x)
  k <- length(K)
  i <- min(K)
  idxK <- which(K == i)
  G <- par2Gamma(par)

  if (NROW(G) != d) {
    stop("The length of par must be d * (d - 1) / 2.")
  }

  S <- Gamma2Sigma(G, k = i, full = TRUE)
  if (k > 1) {
    SK <- S[K[-idxK], K[-idxK]]
    cholSK <- chol(SK)
    SKm1 <- chol2inv(cholSK)
    logdetSK <- 2 * sum(log(diag(cholSK)))
    idxK <- which(K == i)
    yK <- (log(x[K] / x[i]) + G[K, i] / 2)[-idxK]
    logdvK <- -sum(log(x[K])) - log(x[i]) - ((k - 1) / 2) * log(2 * pi) - 1 / 2 * logdetSK - 1 / 2 * t(yK) %*% SKm1 %*% yK
    SnK <- S[-K, -K]
    SnKK <- S[-K, K[-idxK]]
    SKnK <- t(SnKK)
    muCondK <- -G[-K, i] / 2 + SnKK %*% SKm1 %*% yK
    if (k < d - 1) {
      SCondK <- SnK - SnKK %*% SKm1 %*% SKnK
    }
    if (k == d - 1) {
      SCondK <- SnK - t(SnKK) %*% SKm1 %*% t(SKnK)
    }
    logdvnK <- log(mvtnorm::pmvnorm(upper = c(log(x[-K] / x[i]) - muCondK), sigma = SCondK)[1])
    logdv <- logdvK + logdvnK
  }
  if (k == 1) {
    logdvK <- -2 * log(x[i])
    if (d == 2) {
      logdvnK <- log(stats::pnorm(
        q = c(log(x[-K] / x[i]) + G[-K, i] / 2),
        sd = sqrt(S[-K, -K])
      ))
      names(logdvnK) <- "upper"
    } else {
      logdvnK <- log(mvtnorm::pmvnorm(
        upper = c(log(x[-K] / x[i]) + G[-K, i] / 2),
        sigma = S[-K, -K]
      )[1])
    }
    logdv <- logdvK + logdvnK
  }

  return(logdv)
}


#' Full censored log-likelihood of HR model
#'
#' Computes the full (censored) log-likelihood of HR model.
#'
#' @param data Numeric matrix \eqn{n\times d}{n x d}. It contains
#' observations following a multivariate HR Pareto distribution.
#' @param Gamma Numeric matrix \eqn{n\times d}{n x d}.
#' It represents a variogram matrix \eqn{\Gamma}.
#' @param cens Boolean. If true, then censored log-likelihood is computed.
#' By default, \code{cens = FALSE}.
#'
#' @return Numeric. The full censored log-likelihood of HR model.
#'
#' @keywords internal
logLH_HR <- function(data, Gamma, cens = FALSE) {
  if (is.vector(data)) {
    d <- length(data)
    n <- 1
    data <- t(as.matrix(data))
  }
  if (is.matrix(data)) {
    d <- NCOL(data)
    n <- NROW(data)
  }
  par <- Gamma2par(Gamma)

  # if cens = FALSE (default)
  if (!cens) {
    return(-n * log(V_HR(x = rep(1, times = d), par = par))
           + sum(logdV_HR(x = data, par = par)))
  }

  # if cens = TRUE
  p <- rep(1, d)
  data.p <- censor(data, p)
  r <- nrow(data.p)

  L <- apply(data.p > matrix(p, ncol = d, nrow = r, byrow = TRUE), 1, which)
  I <- which(lapply(L, length) > 0 & lapply(L, length) < d)
  J <- which(lapply(L, length) == d)

  if (length(I) > 0) {
    y1 <- mapply(logdVK_HR,
                 x = as.list(data.frame(t(data.p)))[I], K = L[I],
                 MoreArgs = list(par = par)
    )
  } else {
    y1 <- 0
  }

  if (length(J) > 0) {
    y2 <- logdV_HR(x = data.p[J, ], par = par)
  } else {
    y2 <- 0
  }
  return(sum(y1) + sum(y2) - (length(I) + length(J)) * log(V_HR(p, par = par)))
}


#' Parameter fitting for multivariate Huesler--Reiss Pareto distribution
#'
#' Fits the parameters of a multivariate Huesler--Reiss Pareto distribution
#' using (censored) likelihood estimation.
#'
#'
#' If \code{graph = NULL}, then the parameters of a \eqn{d \times d}{d x d}
#' parameter matrix \eqn{\Gamma} of a Huesler--Reiss Pareto distribution are fitted.
#' If \code{graph} is provided, then the conditional independence
#' structure of this graph is assumed and the parameters on the edges are fitted.
#' In both cases the full likelihood is used and therefore this function should only
#' be used for small dimensions, say, \eqn{d<5}. For models in higher dimensions
#' fitting can be done separately on the cliques; see \code{\link{fmpareto_graph_HR}}.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, \code{cens = FALSE}.
#' @param init Numeric vector. Initial parameter values in the optimization. If
#' \code{graph} is given, then the entries should correspond to the edges of the \code{graph}.
#' @param fixParams Numeric vector. Indices of the parameter vectors that are kept
#' fixed during the optimization. Default is `integer(0)`.
#' @param maxit Positive integer. The maximum number of iterations in the
#' optimization.
#' @param graph Graph object from \code{igraph} package or \code{NULL}.
#' If provided, the \code{graph} must be an undirected block graph, i.e., a decomposable, connected
#' graph with singleton separator sets.
#' @param method String. A valid optimization method used by the function
#' \code{\link[stats]{optim}}. By default, \code{method = "BFGS"}.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{convergence}: Logical. Indicates whether the optimization converged or not.
#' \item \code{par}: Numeric vector. Optimized parameters and fixed parameters.
#' \item \code{par_opt}: Numeric. Optimized parameters.
#' \item \code{Gamma}: Numeric matrix \eqn{d \times d}{d x d}. Fitted variogram
#' matrix.
#' \item \code{nllik}: Numeric. Optimized value of the negative log-likelihood function.
#' \item \code{hessian}: Numeric matrix. Estimated Hessian matrix of the
#' estimated parameters.
#' }
#'
#' @keywords internal
fmpareto_HR <- function(data,
                        p = NULL,
                        cens = FALSE,
                        init,
                        fixParams = integer(0),
                        maxit = 100,
                        graph = NULL,
                        method = "BFGS") {
  if (!is.null(p)) {
    # if p provided -> data not Pareto -> to convert
    data <- data2mpareto(data, p)
  } else {
    # if p not provided -> data already Pareto
    data <- data
  }

  # censoring at 1 since data already normalized
  p <- 1
  d <- ncol(data)
  if (length(p) == 1) {
    p <- rep(p, d)
  }

  # convert vector of fixed parameters to logical if necessary
  if(!is.logical(fixParams)){
    fixParams <- seq_along(init) %in% fixParams
  }

  # negative log likelihood function
  if (cens) {
    # censor below the (multivariate) threshold
    data.p <- censor(data, p)
    r <- nrow(data.p)

    L <- apply(data.p > matrix(p, ncol = d, nrow = r, byrow = TRUE), 1, which)

    if (is.matrix(L)) {
      L <- split(t(L), 1:r)
    }

    I <- which(lapply(L, length) > 0 & lapply(L, length) < d)
    J <- which(lapply(L, length) == d)

    nllik <- function(par) {
      # combine par with fixed parameters
      if(any(fixParams)){
        par_full <- init
        par_full[!fixParams] <- par
        par <- par_full
      }

      if (!is.null(graph)) {
        Gtmp <- complete_Gamma(par, graph, allowed_graph_type = 'decomposable')
        par <- Gtmp[upper.tri(Gtmp)]
      }

      G <- par2Gamma(par)
      S <- Gamma2Sigma(G, k = 1)

      if (any(par <= 0) | !matrixcalc::is.positive.definite(S)) {
        return(10^50)
      } else {
        if (length(I) > 0) {
          y1 <- mapply(logdVK_HR,
                       x = as.list(data.frame(t(data.p)))[I],
                       K = L[I], MoreArgs = list(par = par)
          )
        }
        else {
          y1 <- 0
        }
        if (length(J) > 0) {
          y2 <- logdV_HR(x = data.p[J, ], par = par)
        }
        else {
          y2 <- 0
        }
        y <- sum(y1) + sum(y2) - (length(I) + length(J)) * log(V_HR(p, par = par))
        return(-y)
      }
    }
  }
  else {
    r <- nrow(data)
    L <- apply(data > matrix(p, ncol = d, nrow = r, byrow = TRUE), 1, which)

    if (is.matrix(L)) {
      L <- split(t(L), 1:r)
    }

    I <- which(lapply(L, length) > 0) # 1:r
    nllik <- function(par) {
      # combine par with fixed parameters
      if(any(fixParams)){
        par_full <- init
        par_full[!fixParams] <- par
        par <- par_full
      }

      if (!is.null(graph)) {
        Gtmp <- complete_Gamma(par, graph, allowed_graph_type = 'decomposable')
        par <- Gamma2par(Gtmp)
      }

      G <- par2Gamma(par)
      S <- Gamma2Sigma(G, k = 1)

      if (any(par <= 0) || !matrixcalc::is.positive.definite(S)) {
        return(10^50)
      } else {
        if (length(I) > 0) {
          y1 <- logdV_HR(x = data[I, ], par = par)
        }
        else {
          y1 <- 0
        }
        y <- sum(y1) - length(I) * log(V_HR(p, par = par))
        return(-y)
      }
    }
  }

  init_opt <- init[!fixParams]

  # optimize likelihood
  opt <- stats::optim(
    init_opt,
    nllik,
    hessian = TRUE,
    control = list(maxit = maxit),
    method = method
  )

  par <- init
  par[!fixParams] <- opt$par

  z <- list()
  z$convergence <- opt$convergence
  z$par <- par
  z$par_opt <- opt$par
  if (is.null(graph)) {
    z$Gamma <- par2Gamma(z$par)
  } else {
    z$Gamma <- complete_Gamma(graph = graph, Gamma = z$par)
  }
  z$nllik <- opt$value
  z$hessian <- opt$hessian
  return(z)
}



#' Get Cliques and Separators of a graph
#'
#' Finds all cliques, separators, and (recursively) separators of separators
#' in a graph.
#'
#' @param graph An [igraph::graph] object
#' @return A list of vertex sets that represent the cliques and (recursive)
#' separators of `graph`
#'
#' @keywords internal
get_cliques_and_separators <- function(graph){
  # start with maximal cliques as new cliques:
  newCliques <- igraph::max_cliques(graph)
  cliques <- list()

  while(length(newCliques)>0){
    # compute all separators between the new cliques:
    separators <- list()
    for(i in seq_along(newCliques)){
      for(j in seq_len(i-1)){
        sep <- intersect(newCliques[[i]], newCliques[[j]])
        if(length(sep)>1){
          separators <- c(separators, list(sep))
        }
      }
    }
    # add new cliques to cliques-list:
    # (separators are added after the next iteration)
    cliques <- c(newCliques, cliques)

    # use separators as new cliques:
    newCliques <- setdiff(separators, cliques)
  }

  return(cliques)
}


#' Experimental: Fit graph using empirical Theta matrix
#'
#' Fits a graph to an empirical Gamma matrix by computing the corresponding
#' Theta matrix using [Gamma2Theta()] and greedily chooses `m` edges that
#' correspond to high absolute values in Theta.
#'
#' @param data The (standardized) data from which to compute Gamma
#' @param m The number of edges to add, defaults to the number of edges in a tree
#' @param Gamma_emp The empirical Gamma matrix
#' (can be `NULL` if `data` is given)
#'
#' @return A list containing an [igraph::graph] object and a fitted `Gamma` matrix
#'
#' @keywords internal
fit_graph <- function(data, m=NULL, Gamma_emp=NULL){

  if(is.null(Gamma_emp)){
    Gamma_emp <- emp_vario(data)
  }

  Theta_emp <- Gamma2Theta(Gamma_emp)

  d <- nrow(Gamma_emp)

  if(is.null(m)){
    m <- d-1
  } else if(m < d-1){
    stop('m must be at least d-1!')
  }

  g_c <- igraph::make_full_graph(d)
  weights <- Theta_emp[igraph::as_edgelist(g_c)]
  g <- igraph::mst(g_c, weights)

  missing_edges <- igraph::as_edgelist(igraph::complementer(g))
  missing_weights <- Theta_emp[missing_edges]

  ord <- order(abs(missing_weights), decreasing = TRUE)

  new_edges <- missing_edges[ord[seq_len(m - (d-1))],]

  g <- igraph::add_edges(g, t(new_edges))

  Gamma_g <- complete_Gamma(Gamma_emp, g, N=100000, tol=1e-6, check_tol=100)

  return(list(
    graph = g,
    Gamma = Gamma_g
  ))
}



ml_weight_matrix <- function(data, cens = FALSE, p = NULL){
  ## numeric_matrix boolean numeric -> list
  ## produce a named list made of
  ## - llh_hr: loglikelihood values
  ## - est_gamma: estimated parameters

  # Helpers
  llh_cens <- function(x, data, G_emp) {
    ## numeric_vector numeric_matrix numeric_matrix -> numeric_vector
    ## produce parameter estimates and loglikelihood value for censored HR

    fmpareto_obj <- fmpareto_HR(data = data[, x],
                                init = G_emp[x[1], x[2]],
                                cens = cens)
    par.est <- fmpareto_obj$par
    llh_hr <- -(fmpareto_obj$nllik
                - 2 * (sum(log(data[which(data[, x[1]] > 1), x[1]]))
                       + sum(log(data[which(data[, x[2]] > 1), x[2]]))))
    c(par = par.est, llh_hr = llh_hr)

  }

  llh_uncens <- function(x, data, G_emp) {
    ## numeric_vector numeric_matrix numeric_matrix -> numeric_vector
    ## produce parameter estimates and loglikelihood value for uncensored HR

    par.est <- fmpareto_HR(
      data = data[, x], init = G_emp[x[1], x[2]],
      cens = cens
    )$par

    llh_hr <- logLH_HR(data = data[, x], Gamma = par2Gamma(par.est)) +
      2 * (sum(log(data[, x[1]])) + sum(log(data[, x[2]])))

    c(par = par.est, llh_hr = llh_hr)

  }

  # Standardize data
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  # Set up some variables
  d <- ncol(data.std)
  G_emp <- emp_vario(data = data.std)
  res <- which(upper.tri(matrix(nrow = d, ncol = d)), arr.ind = TRUE)

  # Fig loglikelihood
  if (cens) {
    bivLLH <- apply(res[, 1:2], 1, llh_cens, data.std, G_emp)
  } else {
    bivLLH <- apply(res[, 1:2], 1, llh_uncens, data.std, G_emp)
  }

  bivLLH.mat <- -par2Gamma(bivLLH["llh_hr", ])
  est_gamma <- par2Gamma(bivLLH["par", ])

  return(list(
    llh_hr = bivLLH.mat,
    est_gamma = est_gamma
  ))
}
