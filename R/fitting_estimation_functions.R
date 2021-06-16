#' Empirical estimation of extremal correlation \eqn{\chi}
#'
#' Estimates the \eqn{d}-dimensional extremal correlation coefficient \eqn{\chi} empirically.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#'
#' @return Numeric. The empirical \eqn{d}-dimensional extremal correlation coefficient \eqn{\chi}
#' for the \code{data}.
#' @examples
#' n <- 100
#' d <- 2
#' p <- .8
#' G <- cbind(
#'   c(0, 1.5),
#'   c(1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = G)
#' emp_chi_multdim(my_data, p)
#' @export
emp_chi_multdim <- function(data, p = NULL) {
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  d <- ncol(data)
  data <- stats::na.omit(data)

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }



  rowmin <- apply(data.std, 1, min)
  chi <- mean(sapply(1:ncol(data.std), FUN = function(i) {
    mean(rowmin[which(data.std[, i] > 1)] > 1)
  }), na.rm = TRUE)
  return(chi)
}


#' Empirical estimation of extremal correlation matrix \eqn{\chi}
#'
#' Estimates empirically the matrix of bivariate extremal correlation coefficients \eqn{\chi}.
#'
#' @inheritParams emp_chi_multdim
#'
#' @return Numeric matrix \eqn{d\times d}{d x d}. The matrix contains the
#' bivariate extremal coefficients \eqn{\chi_{ij}}, for \eqn{i, j = 1, ..., d}.
#' @examples
#' n <- 100
#' d <- 4
#' p <- .8
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = Gamma)
#' emp_chi(my_data, p)
#' @export
emp_chi <- function(data, p = NULL) {
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  n <- nrow(data)
  d <- ncol(data)

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }


  ind <- data.std > 1
  ind_mat <- matrix(colSums(ind), byrow = TRUE, ncol = d, nrow = d)
  crossprod(ind, ind) / (1 / 2 * (ind_mat + t(ind_mat)))
}


emp_chi_deprecated <- function(data, p) {
  d <- ncol(data)
  res <- as.matrix(expand.grid(1:d, 1:d))
  res <- res[res[, 1] > res[, 2], , drop = FALSE]
  chi <- apply(res, 1, function(x) {
    emp_chi_multdim(cbind(data[, x[1]], data[, x[2]]), p = p)
  })
  chi.mat <- matrix(NA, ncol = d, nrow = d)
  chi.mat[res] <- chi
  chi.mat[res[, 2:1, drop = FALSE]] <- chi
  diag(chi.mat) <- 1

  return(chi.mat)
}



#' Estimation of the variogram matrix \eqn{\Gamma} of the Huesler--Reiss distribution
#'
#' Estimates the variogram of the Huesler--Reiss distribution empirically.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param k Integer between 1 and \eqn{d}. Component of the multivariate
#' observations that is conditioned to be larger than the threshold \code{p}.
#' If \code{NULL} (default), then an average over all \code{k} is returned.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d}. The estimated
#' variogram of the Huesler--Reiss distribution.
#' @export
emp_vario <- function(data, k = NULL, p = NULL) {

  # helper ####
  G.fun <- function(i, data) {
    idx <- which(data[, i] > 1)
    if (length(idx) > 1) {
      xx <- Sigma2Gamma(stats::cov(log(data[idx, ])), full = TRUE)
    } else {
      xx <- matrix(NA, d, d)
    }
    return(xx)
  }

  # body ####
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  d <- ncol(data)
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  if (!is.null(k)) {
    G <- G.fun(k, data.std)

    if (any(is.na(G))) {
      warning(paste(
        "Produced NA matrix since there are no exceedances in the component k =",
        k
      ))
    }
  } else {

    # take the average
    row_averages <- rowMeans(sapply(1:d, FUN = function(i) {
      G.fun(i, data.std)
    }), na.rm = TRUE)
    G <- matrix(row_averages, nrow = d, ncol = d)
  }

  return(G)
}



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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' \item \code{par}: Numeric vector. Optimized parameters.
#' \item \code{Gamma}: Numeric matrix \eqn{d \times d}{d x d}. Fitted variogram
#' matrix.
#' \item \code{nllik}: Numeric. Optimized value of the negative log-likelihood function.
#' \item \code{hessian}: Numeric matrix. Estimated Hessian matrix of the
#' estimated parameters.
#' }
#'
#' @noRd
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
    fixParams <- (1:d) %in% fixParams
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



#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on block graphs
#'
#' Fits the parameters of a multivariate Huesler--Reiss Pareto distribution using (censored) likelihood estimation.
#' Fitting is done separately on the cliques of the block graph. If  \code{edges_to_add}
#' are provided, then these edges are added in a greedy search to the original \code{graph},
#' such that in each step the likelihood is improved maximally and the new graph stays in the
#' class of block graphs. See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, \code{cens = FALSE}.
#' @param graph Graph object from \code{igraph} package. The \code{graph} must be an undirected block graph, i.e., a decomposable, connected
#' graph with singleton separator sets.
#' @param edges_to_add Numeric matrix \eqn{m\times 2}{m x 2}, where \eqn{m} is
#' the number of edges that are tried to be added in the greedy search.
#' By default, \code{edges_to_add = NULL}.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{graph}: Graph object from \code{igraph} package. If \code{edges_to_add} are provided,
#' then this is a list of the resulting graphs in each step of the greedy search.
#' \item \code{Gamma}: Numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}.
#' If \code{edges_to_add} are provided,
#' then this is a list of the estimated variogram matrices in each step of the greedy search.
#' \item \code{AIC}: (only if \code{edges_to_add} are provided) List of AIC values of the fitted models
#' in each step of the greedy search.
#' \item \code{edges_added}: (only if \code{edges_to_add} are provided) Numeric matrix \eqn{m'\times 2}{m' x 2}, where
#' the \eqn{m'\leq m}{m'<=m} rows contain the edges that were added in the greedy search.
#' }
#'
#' @examples
#' ## Fitting a 4-dimensional HR distribution
#'
#' my_graph <- igraph::graph_from_adjacency_matrix(
#'   rbind(
#'     c(0, 1, 0, 0),
#'     c(1, 0, 1, 1),
#'     c(0, 1, 0, 0),
#'     c(0, 1, 0, 0)
#'   ),
#'   mode = "undirected"
#' )
#' n <- 100
#' Gamma_vec <- c(.5, 1.4, .8)
#' complete_Gamma(Gamma = Gamma_vec, graph = my_graph) ## full Gamma matrix
#' edges_to_add <- rbind(c(1, 3), c(1, 4), c(3, 4))
#'
#' set.seed(123)
#' my_data <- rmpareto_tree(n, "HR", tree = my_graph, par = Gamma_vec)
#' my_fit <- fmpareto_graph_HR(my_data,
#'   graph = my_graph,
#'   p = NULL, cens = FALSE, edges_to_add = edges_to_add
#' )
#' @references
#'  \insertAllCited{}
#' @export
fmpareto_graph_HR <- function(data, graph, p = NULL, cens = FALSE, edges_to_add = NULL) {

  # set up main variables
  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # check if it is directed
  if (igraph::is_directed(graph)) {
    warning("The given graph is directed. Converted to undirected.")
    graph <- igraph::as.undirected(graph)
  }

  # check if it is connected
  is_connected <- igraph::is_connected(graph)

  if (!is_connected) {
    stop("The given graph is not connected.")
  }

  # check if graph is decomposable
  is_decomposable <- igraph::is_chordal(graph)$chordal
  if (!is_decomposable) {
    stop("The given graph is not decomposable (i.e., chordal).")
  }

  # check if it is block graph
  cli <- igraph::max_cliques(graph)
  ncli <- length(cli)
  min_sep <- 0

  for (i in 1:ncli) {
    cli1 <- cli[[i]]
    for (j in 1:ncli) {
      if (j <= i) {
        next
      }
      cli2 <- cli[[j]]

      min_sep <- max(min_sep, length(intersect(cli1, cli2)))

      if (min_sep > 1) {
        break
      }
    }
  }

  if (min_sep > 1) {
    stop("The given graph is not a block graph.")
  }

  # check if the number of nodes in the graph matches the number
  # of variables in the data matrix
  nnodes <- igraph::vcount(graph)
  if (nnodes != NCOL(data)) {
    stop(paste(
      "The number of nodes in the graph doesn't match with the number",
      "of variables (i.e., columns) in the data matrix."
    ))
  }

  # check if you need to rescale data or not
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  l <- 1

  graph.cur <- list()
  graph.cur[[l]] <- graph
  Ghat <- list()
  Ghat[[l]] <- matrix(NA, nrow = nnodes, ncol = nnodes)

  # loop through all cliques
  for (i in 1:ncli) {
    # pick the curren cliques
    cli.idx <- cli[[i]]
    # how many nodes in the current cliques?
    cli.len <- length(cli.idx)
    # compute marginal pareto, on the nodes of the current clique
    data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)

    G.est <- emp_vario(data = data.cli)
    init <- Gamma2par(G.est)
    Ghat[[l]][cli.idx, cli.idx] <- fmpareto_HR(
      data = data.cli,
      init = init, cens = cens
    )$Gamma
  }

  Ghat[[l]] <- complete_Gamma(graph = graph.cur[[l]], Gamma = Ghat[[l]])

  # if you want to add some edges
  if (!is.null(edges_to_add)) {
    # check if edges_to_add is vector
    if (is.vector(edges_to_add)) edges_to_add <- t(as.matrix(edges_to_add))

    # check if any proposed edge is already in the given graph
    adj_mat <- igraph::as_adjacency_matrix(graph, sparse = FALSE) > 0

    m <- nrow(edges_to_add)
    check_new_edges <- 0
    for (k in 1:m) {
      current_edge <- edges_to_add[k, ]

      is_already_edge <- adj_mat[current_edge[1], current_edge[2]] |
        adj_mat[current_edge[2], current_edge[1]]

      if (is_already_edge) {
        break
      }
    }

    if (is_already_edge) {
      stop(paste(
        "The argument edges_to_add cannot contain edges already",
        "present in the given graph."
      ))
    }


    stop.flag <- FALSE
    AIC <- 2 * igraph::ecount(graph.cur[[l]]) - 2 * logLH_HR(
      data = data.std,
      Gamma = Ghat[[l]], cens = cens
    )
    edges_added <- c()

    while (length(edges_to_add) != 0 & stop.flag == FALSE) {
      m <- nrow(edges_to_add)
      AIC.tmp <- rep(NA, times = m)
      Ghat.tmp <- list()

      # go through proposed edges one after the other while retaining a block
      # graph
      # m number of proposed edges
      for (k in 1:m) {
        # current temporary graph
        Ghat.tmp[[k]] <- Ghat[[l]]
        # add the current proposed edge to the graph
        graph.tmp <- igraph::add_edges(
          graph = graph.cur[[l]],
          edges = edges_to_add[k, ]
        )

        # if the obtained graph is decomposable
        if (igraph::is_chordal(graph.tmp)$chordal) {
          # find list of max cliques
          cli <- igraph::max_cliques(graph.tmp)
          # find in which clique the new proposed edge is. It can be in at most
          # one clique, otherwise, the original graph were not decomposable.
          intersections <-
            sapply(cli, FUN = function(x) length(intersect(x, edges_to_add[k, ])) == 2)
          ii <- which(intersections == TRUE)


          # only in the clique itself the separator can be of size > 1
          if (sum(sapply(cli, FUN = function(x) {
            length(intersect(x, cli[[ii]])) > 1
          })) == 1) {
            cat("\nTry edge", edges_to_add[k, ])
            cli.idx <- cli[[ii]]
            cli.len <- length(cli.idx)
            data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)

            G.est <- emp_vario(data = data.cli)
            init <- Gamma2par(G.est)
            Ghat.tmp[[k]][cli.idx, cli.idx] <- fmpareto_HR(
              data = data.cli,
              init = init, cens = cens
            )$Gamma
            Ghat.tmp[[k]] <- complete_Gamma(graph = graph.tmp, Gamma = Ghat.tmp[[k]])
            AIC.tmp[k] <- 2 * igraph::ecount(graph.tmp) -
              2 * logLH_HR(data = data.std, Gamma = Ghat.tmp[[k]], cens = cens)
          }
        }
      }
      if (!all(is.na(AIC.tmp))) {
        add.idx <- which(AIC.tmp == min(AIC.tmp, na.rm = TRUE))
        cat("\nAdded edge ", edges_to_add[add.idx, ])
        l <- l + 1
        graph.cur[[l]] <-
          igraph::add_edges(graph = graph.cur[[l - 1]], edges = edges_to_add[add.idx, ])
        graph.cur[[l]] <- set_graph_parameters(graph.cur[[l]])
        Ghat[[l]] <- Ghat.tmp[[add.idx]]
        AIC <- c(AIC, AIC.tmp[add.idx])
        edges_added <- rbind(edges_added, t(as.matrix(edges_to_add[add.idx, ])))
        edges_to_add <- edges_to_add[-add.idx, ]
      }
      if (all(is.na(AIC.tmp))) stop.flag <- TRUE
    }
    return(list(
      graph = graph.cur,
      Gamma = Ghat, AIC = AIC, edges_added = edges_added
    ))
  }

  return(list(graph = set_graph_parameters(graph), Gamma = Ghat[[1]]))
}


#' Fitting extremal minimum spanning tree
#'
#' Fits an extremal minimum spanning tree, where the edge weights are:
#' \itemize{
#' \item negative maximized log-likelihoods of the bivariate Huesler--Reiss
#' distributions, if `method = "ML"`. See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#' \item empirical extremal variogram, if `method = "vario"`.
#' \item empirical extremal correlation, if `method = "chi"`.
#' }
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param method One of `"ML", "vario", "chi"`.
#' @param cens Logical. This argument is considered only if `method = "ML"`.
#' If `TRUE`, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{tree}: Graph object from \code{igraph} package. The fitted minimum spanning tree.
#' \item \code{Gamma}: Numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}
#' corresponding to the fitted minimum spanning tree.
#' }
#'
#' @examples
#' ## Fitting a 4-dimensional HR MST tree
#'
#' my_graph <- igraph::graph_from_adjacency_matrix(
#'   rbind(
#'     c(0, 1, 0, 0),
#'     c(1, 0, 1, 1),
#'     c(0, 1, 0, 0),
#'     c(0, 1, 0, 0)
#'   ),
#'   mode = "undirected"
#' )
#' n <- 100
#' Gamma_vec <- c(.5, 1.4, .8)
#' complete_Gamma(Gamma = Gamma_vec, graph = my_graph) ## full Gamma matrix
#'
#' set.seed(123)
#' my_data <- rmpareto_tree(n, "HR", tree = my_graph, par = Gamma_vec)
#' my_fit <- mst_HR(my_data, p = NULL, cens = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
#'
emst <- function(data, p = NULL, method = c("ML", "vario", "chi"),
                 cens = FALSE) {

  # Validate arguments
  method <- match.arg(method)

  # Check if you need to rescale data or not
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  # Estimate weight matrix
  if (method == "ML") {
    res <- ml_weight_matrix(data = data.std, cens = cens)
    weight_matrix <- res$llh_hr
    estimated_gamma <- res$est_gamma

  } else if (method == "vario") {
    weight_matrix <- emp_vario(data = data.std)
    estimated_gamma <- weight_matrix

  } else if (method == "chi") {
    weight_matrix <- - emp_chi(data = data.std)
    estimated_gamma <- chi2Gamma(- weight_matrix)

  }

  # Estimate tree
  graph.full <- igraph::make_full_graph(ncol(data.std))
  mst.tree <- igraph::mst(
    graph = graph.full,
    weights = weight_matrix[igraph::ends(graph.full, igraph::E(graph.full))],
    algorithm = "prim"
  )

  # Set graphical parameters
  mst.tree <- set_graph_parameters(mst.tree)

  # Return tree and completed Gamma
  return(list(
    tree = mst.tree,
    Gamma = complete_Gamma(graph = mst.tree, Gamma = estimated_gamma)
  ))
}

mst_HR <- function(data, p = NULL, cens = FALSE) {

  # check if you need to rescale data or not
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  # set some params
  n <- nrow(data.std)
  d <- ncol(data.std)
  graph.full <- igraph::make_full_graph(d)

  # compute weight matrix
  G.emp <- emp_vario(data = data.std)
  res <- which(upper.tri(matrix(nrow = d, ncol = d)), arr.ind = TRUE)
  if (cens) {
    bivLLH <- apply(res[, 1:2], 1, function(x) {
      fmpareto_obj <- fmpareto_HR(data = data.std[, x],
                                  init = G.emp[x[1], x[2]],
                                  cens = cens)
      par.est <- fmpareto_obj$par
      llh_hr <- -(fmpareto_obj$nllik
        - 2 * (sum(log(data.std[which(data.std[, x[1]] > 1), x[1]]))
        + sum(log(data.std[which(data.std[, x[2]] > 1), x[2]]))))
      c(par = par.est, llh_hr = llh_hr)
    })
  }


  if (!cens) {
    bivLLH <- apply(res[, 1:2], 1, function(x) {
      par.est <- fmpareto_HR(
        data = data.std[, x], init = G.emp[x[1], x[2]],
        cens = cens
      )$par

      llh_hr <- logLH_HR(data = data.std[, x], Gamma = par2Gamma(par.est)) +
        2 * (sum(log(data.std[, x[1]])) + sum(log(data.std[, x[2]])))

      c(par = par.est, llh_hr = llh_hr)
    })
  }

  bivLLH.mat <- par2Gamma(bivLLH["llh_hr", ])

  # Estimated tree
  mst.tree <- igraph::mst(
    graph = graph.full, weights =
      -bivLLH.mat[igraph::ends(
        graph.full,
        igraph::E(graph.full)
      )],
    algorithm = "prim"
  )

  # set graphical parameters
  mst.tree <- set_graph_parameters(mst.tree)

  # Estimated Gamma
  est_Gamma <- par2Gamma(bivLLH["par", ])

  # return tree
  return(list(
    tree = mst.tree,
    Gamma = complete_Gamma(graph = mst.tree, Gamma = est_Gamma)
  ))
}



ml_weight_matrix <- function(data, cens){
  ## numeric_matrix boolean -> list
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

  # Set up some variables
  d <- ncol(data)
  G_emp <- emp_vario(data = data)
  res <- which(upper.tri(matrix(nrow = d, ncol = d)), arr.ind = TRUE)

  # Fig loglikelihood
  if (cens) {
    bivLLH <- apply(res[, 1:2], 1, llh_cens, data, G_emp)
  } else {
    bivLLH <- apply(res[, 1:2], 1, llh_uncens, data, G_emp)
  }

  bivLLH.mat <- -par2Gamma(bivLLH["llh_hr", ])
  est_gamma <- par2Gamma(bivLLH["par", ])

  return(list(
    llh_hr = bivLLH.mat,
    est_gamma = est_gamma
  ))
}
