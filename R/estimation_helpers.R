#' Compute exponent measure
#'
#' Computes the exponent measure of HR distribution.
#'
#' @param x Numeric vector with \eqn{d} positive elements
#' where the exponent measure is to be evaluated.
#' @param par \eqn{d \times d}{d x d} variogram matrix, or numeric vector with
#' \eqn{\frac{d(d - 1)}{2}}{d x (d - 1) / 2} elements,
#' representing the upper triangular portion of a
#' variogram matrix \eqn{\Gamma}.
#'
#' @return Numeric. The exponent measure of the HR distribution.
#'
#' @keywords internal
V_HR <- function(x, Gamma = NULL, Theta = NULL) {
  # Check/convert parameters
  if(is.null(Gamma)){
    Theta <- par2Theta()
    Gamma <- Theta2Gamma(Theta)
  } else{
    Gamma <- par2Gamma(Gamma, TRUE)
  }
  if (NROW(Gamma) != d) {
    stop("`par` must be a vector of length `d * (d - 1) / 2` or a d x d matrix.")
  }

  d <- length(x)

  # helper function
  f1 <- function(i, x) {
    S <- Gamma2Sigma(Gamma, k = i)
    return(1 / x[i] * mvtnorm::pmvnorm(
      upper = (log(x / x[i]) + Gamma[, i] / 2)[-i],
      mean = rep(0, d - 1), sigma = S
    )[1])
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
logdV_HR <- function(x, Gamma = NULL, Theta = NULL) {
  # Make sure x is positive
  if (any(is_leq(x, 0))) {
    stop("The elements of x must be positive.")
  }

  # Make sure x is a matrix
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  # Convert parameter-vector to Gamma matrix
  d <- ncol(x)
  G <- par2Gamma(Gamma, TRUE)
  if (NROW(G) != d) {
    stop("The length of par must be d * (d - 1) / 2.")
  }

  # Compute likelihood
  k <- 1
  Sigma_k <- Gamma2Sigma(G, k = k)
  cholS <- chol(Sigma_k)
  Theta_k <- chol2inv(cholS)
  logdetSigma_k <- 2 * sum(log(diag(cholS)))
  if (is.matrix(x)) {
    yTilde_k <- (t(t(log(x / x[, k])) + G[, k] / 2))[, -k, drop = FALSE]
    logdv <- (
      (-1) * rowSums(log(x))
      - log(x[, k])
      - ((d - 1) / 2) * log(2 * pi)
      - 1 / 2 * logdetSigma_k
      - 1 / 2 * fast_diag(yTilde_k, Theta_k)
    )
  }
  return(logdv)
}

#' Fast computation of diag(y %*% M %*% t(y))
#' 
#' @param y Numeric matrix
#' @param M Numeric matrix
#' @return Numeric vector
fast_diag <- function(y, M) {
  n <- nrow(y)
  if(n == 0){
    return(numeric(0))
  }
  sapply(seq_len(n), function(i) {
    u <- y[i, , drop = FALSE]
    u %*% M %*% t(u)
  })
}



#' Compute censored exponent measure
#'
#' Computes the (censored) exponent measure density of HR distribution.
#'
#' @param x Numeric vector with `d` positive elements
#' where the censored exponent measure is to be evaluated.
#' @param K Integer vector, subset of \eqn{\{1, \dots, d\}}{{1, ..., d}},
#' the index set that is not censored.
#' Or: Logical vector of length `d`, indicating entries that are not censored.
#' @inheritParams V_HR
#'
#' @return Numeric. The censored exponent measure of the HR distribution.
#' If no entries are censored, the result of `logdV_HR(x, par` is returned.
#'
#' @keywords internal
logdVK_HR <- function(x, K, par) {
  if (any(is_leq(x, 0))) {
    stop("The elements of x must be positive.")
  }
  
  # Convert logical K to numeric indices
  if(is.logical(K)){
    K <- which(K)
  }

  # return normal density, if no entries are censored
  if(length(K) == length(x)){
    return(logdV_HR(x, par))
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
  # Todo: this can be expressed as a parallel minimum (`pmin`)?
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

    fmpareto_obj <- fmpareto_HR_MLE_Gamma(
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

    par.est <- fmpareto_HR_MLE_Gamma(
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

