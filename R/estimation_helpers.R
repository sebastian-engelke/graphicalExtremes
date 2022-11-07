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
#' @param x Numeric matrix \eqn{n\times d}{n x d} or vector with \eqn{d} elements.
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
#' By default, `cens = FALSE`.
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
#' If `graph = NULL`, then the parameters of a \eqn{d \times d}{d x d}
#' parameter matrix \eqn{\Gamma} of a Huesler--Reiss Pareto distribution are fitted.
#' If `graph` is provided, then the conditional independence
#' structure of this graph is assumed and the parameters on the edges are fitted.
#' In both cases the full likelihood is used and therefore this function should only
#' be used for small dimensions, say, \eqn{d<5}. For models in higher dimensions
#' fitting can be done separately on the cliques; see [fmpareto_graph_HR].
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#' @param init Numeric vector. Initial parameter values in the optimization. If
#' `graph` is given, then the entries should correspond to the edges of the `graph`.
#' @param fixParams Numeric vector. Indices of the parameter vectors that are kept
#' fixed during the optimization. Default is `integer(0)`.
#' @param maxit Positive integer. The maximum number of iterations in the
#' optimization.
#' @param graph Graph object from `igraph` package or `NULL`.
#' If provided, the `graph` must be an undirected block graph, i.e., a decomposable, connected
#' graph with singleton separator sets.
#' @param method String. A valid optimization method used by the function
#' [stats::optim]. By default, `method = "BFGS"`.
#'
#' @return List consisting of:
#' \item{`convergence`}{Logical. Indicates whether the optimization converged or not.}
#' \item{`par`}{Numeric vector. Optimized parameters and fixed parameters.}
#' \item{`par_opt`}{Numeric. Optimized parameters.}
#' \item{`Gamma`}{Numeric matrix \eqn{d \times d}{d x d}. Fitted variogram #' matrix.}
#' \item{`nllik`}{Numeric. Optimized value of the negative log-likelihood function.}
#' \item{`hessian`}{Numeric matrix. Estimated Hessian matrix of the #' estimated parameters.}
#'
#' @keywords internal
fmpareto_HR_MLE <- function(
  data,
  p = NULL,
  cens = FALSE,
  init,
  fixParams = integer(0),
  maxit = 100,
  graph = NULL,
  method = "BFGS"
){
  # if p provided -> data not Pareto -> to convert
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
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
          y1 <- mapply(
            logdVK_HR,
            x = as.list(data.frame(t(data.p)))[I],
            K = L[I], MoreArgs = list(par = par)
          )
        } else {
          y1 <- 0
        }
        if (length(J) > 0) {
          y2 <- logdV_HR(x = data.p[J, ], par = par)
        } else {
          y2 <- 0
        }
        y <- sum(y1) + sum(y2) - (length(I) + length(J)) * log(V_HR(p, par = par))
        return(-y)
      }
    }
  } else {
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
        } else {
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


#' Censor dataset
#'
#' Censors each row of matrix `x` with vector `p`.
#'
#' @param x Numeric matrix \eqn{n \times d}{n x d}.
#' @param p Numeric vector with \eqn{d} elements.
#'
#' @return Numeric matrix \eqn{n \times d}{n x d}.
#'
#' @keywords internal
censor <- function(x, p) {
  f2 <- function(x, p) {
    x_is_less <- x <= p
    y <- x
    y[x_is_less] <- p[x_is_less]
    return(y)
  }
  return(t(apply(x, 1, f2, p)))
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

    fmpareto_obj <- fmpareto_HR_MLE(
      data = data[, x],
      init = G_emp[x[1], x[2]],
      cens = cens
    )
    par.est <- fmpareto_obj$par
    llh_hr <- -(
      fmpareto_obj$nllik
      - 2 * (sum(log(data[which(data[, x[1]] > 1), x[1]]))
      + sum(log(data[which(data[, x[2]] > 1), x[2]])))
    )
    return(c(par = par.est, llh_hr = llh_hr))
  }

  llh_uncens <- function(x, data, G_emp) {
    ## numeric_vector numeric_matrix numeric_matrix -> numeric_vector
    ## produce parameter estimates and loglikelihood value for uncensored HR

    par.est <- fmpareto_HR_MLE(
      data = data[, x], init = G_emp[x[1], x[2]],
      cens = cens
    )$par

    llh_hr <- (
      logLH_HR(data = data[, x], Gamma = par2Gamma(par.est))
      + 2 * (sum(log(data[, x[1]])) + sum(log(data[, x[2]])))
    )

    return(c(par = par.est, llh_hr = llh_hr))
  }

  # Standardize data
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  # Set up some variables
  d <- ncol(data)
  G_emp <- emp_vario(data)
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


glasso_mb <- function(data, lambda){
  # Initialize variables
  dd <- ncol(data)
  data <- scale(data)
  adj.est <- array(NA, dim = c(dd, dd, length(lambda)))
  adj.ic.est <- array(NA, dim = c(dd, dd, 3))
  lambda_order <- order(lambda, decreasing = TRUE)
  lambda_dec <- sort(lambda, decreasing = TRUE)
  # there was some subtlety with the ordering since glmnet always gives back
  # in a certain order, that's why I have this here

  # Loop through variables
  for(i in (1:dd)){
    X <- data[,-i]
    Y <- data[,i]
    lasso_fit <- glmnet::glmnet(x = X, y = Y, family = "gaussian", lambda = lambda_dec)
    if(i==1){
      # ensures same lambda sequence is used for different lasso regressions
      lambda_dec <- lasso_fit$lambda
      null.vote <- array(0, dim = c(dd, dd, length(lambda)))
      null.vote.ic <- array(0, dim = c(dd, dd, 3))
    }
    # make sure consistent with default value
    null.vote[i, -i, ] <- null.vote[i, -i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= 1e-10)
    null.vote[-i, i, ] <- null.vote[-i, i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= 1e-10)

    aic.idx <- which.min((1 - lasso_fit$dev.ratio) * lasso_fit$nulldev + aic(nrow(X), ncol(X) + 1) * lasso_fit$df)
    bic.idx <- which.min((1 - lasso_fit$dev.ratio) * lasso_fit$nulldev + bic(nrow(X), ncol(X) + 1) * lasso_fit$df)
    mbic.idx <- which.min((1 - lasso_fit$dev.ratio) * lasso_fit$nulldev + mbic(nrow(X), ncol(X) + 1) * lasso_fit$df)

    null.vote.ic[i, -i, ] <- null.vote.ic[i, -i,] +
      (abs(as.matrix(lasso_fit$beta[,c(aic.idx, bic.idx, mbic.idx)])) <= 1e-10)
    null.vote.ic[-i, i, ] <- null.vote.ic[-i, i,] +
      (abs(as.matrix(lasso_fit$beta[,c(aic.idx, bic.idx, mbic.idx)])) <= 1e-10)
  }
  adj.est[,,lambda_order] <- null.vote <= 1
  adj.ic.est <- null.vote.ic <= 1
  return(list(adj.est=adj.est, adj.ic.est = adj.ic.est))
}

# traditional criteria
# is consistent for a fixed design, fixed p
aic <- function(n, p) 2
bic <- function(n, p) log(n)

# modified BIC of Wang & Leng, JRSSB 2009
# it has a weird justification in the paper
mbic <- function(n, p) log(n) * log(log(p))


try_complete_Gamma <- function(graph, Gamma, key, val){
  ## igraph numeric_matrix character double -> numeric_matrix | NA
  ## tries to call `complete_Gamma`, if it fails returns NULL
  if(!igraph::is.connected(graph)){
    message(paste0(
      'The estimated graph for ',
      key, ' = ', val,
      '} is not connected,  so it is not possible to complete `Gamma`.\n'
    ))
    return(NULL)
  }

  Gamma_comp <- complete_Gamma(graph = graph, Gamma = Gamma)
  graph_comp <- Gamma2graph(Gamma_comp)

  # Check if completed Gamma matches with given graph
  if (!graphs_equal(graph_comp, graph)) {
    message(paste0(
      'The completed Gamma for ',
      key, ' = ', val,
      ' does not match the estimated graph.\n'
    ))
  }

  # Return completed Gamma
  return(Gamma_comp)
}

