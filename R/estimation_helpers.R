#' Compute exponent measure
#'
#' Computes the exponent measure of HR distribution.
#'
#' @param x Numeric vector with `d` positive elements
#' where the exponent measure is to be evaluated.
#' @param Gamma d x d variogram matrix or numeric vector with d(d-1)/2 elements,
#' containing the upper triangular part of a variogram matrix.
#' @param Theta d x d precision matrix or numeric vector with d(d-1)/2 elements,
#' containing the upper triangular part of a precision matrix.
#' 
#' @details Only `Gamma` is needed for the computation. `Theta` is only used to
#' compute `Gamma` if necessary.
#'
#' @return Numeric. The exponent measure of the HR distribution.
#'
#' @keywords internal
V_HR <- function(x, Gamma = NULL, Theta = NULL) {
  # Convert Theta -> Gamma if necessary (Theta is ignored otherwise)
  if(is.null(Gamma)){
    Theta <- par2Theta(Theta, TRUE)
    Gamma <- Theta2Gamma(Theta)
  } else{
    Gamma <- par2Gamma(Gamma, TRUE)
  }
  d <- length(x)
  if (NROW(Gamma) != d) {
    stop("`par` must be a vector of length `d * (d - 1) / 2` or a d x d matrix.")
  }

  # helper function
  f1 <- function(k) {
    Sk <- Gamma2Sigma(Gamma, k = k)
    return(1 / x[k] * mvtnorm::pmvnorm(
      upper = (log(x / x[k]) + Gamma[, k] / 2)[-k],
      mean = rep(0, d - 1),
      sigma = Sk
    )[1])
  }

  return(sum(sapply(1:d, f1)))
}


#' Compute the exponent measure density of HR distribution
#'
#' Computes the exponent measure density of HR distribution.
#'
#' @param x Numeric \nxd matrix or vector with `d` elements.
#' @inheritParams V_HR
#' 
#' @details Both `Gamma` and `Theta` are needed internally, but if one
#' is missing it is computed from the other one.
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
  d <- ncol(x)

  # Check/convert parameters
  k <- 1
  Gamma <- par2Gamma(Gamma, allowMatrix = TRUE, allowNull = TRUE)
  Theta <- par2Theta(Theta, allowMatrix = TRUE, allowNull = TRUE)
  if(is.null(Gamma) && is.null(Theta)){
    stop('Specify at least one of Gamma, Theta.')
  }
  if(is.null(Theta)){ # -> Gamma must be specified
    Sigma_k <- Gamma2Sigma(Gamma, k = k)
    cholS <- chol(Sigma_k)
    Theta_k <- chol2inv(cholS)
    logdetSigma_k <- 2 * sum(log(diag(cholS)))
  } else{
    Theta_k <- Theta[-k, -k, drop=FALSE]
    tmp <- determinant(Theta_k, logarithm = TRUE)
    logdetSigma_k <- (-1) * c(tmp$modulus) # `c()` removes attributes, tmp$sign is always +1
  }
  if(is.null(Gamma)){ # -> Theta must be specified
    Gamma <- Theta2Gamma(Theta)
  }
  if (NROW(Gamma) != d) {
    stop("`par` must be a vector of length `d * (d - 1) / 2` or a d x d matrix.")
  }

  # Compute likelihood
  yTilde_k <- (t(t(log(x / x[, k])) + Gamma[, k] / 2))[, -k, drop = FALSE]
  logdv <- (
    - rowSums(log(x))
    - log(x[, k])
    - ((d - 1) / 2) * log(2 * pi)
    - 1 / 2 * logdetSigma_k
    - 1 / 2 * fast_diag(yTilde_k, Theta_k)
  )
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
#' @param K Integer vector, subset of 1, ..., `d`, the index set that is not censored.
#' Or logical vector of length `d`, indicating entries that are not censored.
#' @inheritParams V_HR
#'
#' @return Numeric. The censored exponent measure of the HR distribution.
#' If no entries are censored, the result of `logdV_HR(x, par` is returned.
#'
#' @keywords internal
logdVK_HR <- function(x, K, Gamma) {
  ## TODO: this function can probably be optimized by allowing calls with multiple observations

  # Convert logical K to numeric indices
  if(is.logical(K)){
    K <- which(K)
  }

  if (any(is_leq(x, 0))) {
    stop("The elements of x must be positive.")
  }

  # return normal density, if no entries are censored
  if(length(K) == length(x)){
    return(logdV_HR(x, Gamma))
  }


  d <- length(x)
  k <- length(K)
  i <- min(K)
  idxK <- which(K == i)

  Gamma <- par2Gamma(Gamma, allowMatrix = TRUE)
  if (NROW(Gamma) != d) {
    stop("The length of par must be d * (d - 1) / 2.")
  }

  S <- Gamma2Sigma(Gamma, k = i, full = TRUE)
  if (k > 1) {
    SK <- S[K[-idxK], K[-idxK]]
    cholSK <- chol(SK)
    SKm1 <- chol2inv(cholSK)
    logdetSK <- 2 * sum(log(diag(cholSK)))
    idxK <- which(K == i)
    yK <- (log(x[K] / x[i]) + Gamma[K, i] / 2)[-idxK]
    logdvK <- -sum(log(x[K])) - log(x[i]) - ((k - 1) / 2) * log(2 * pi) - 1 / 2 * logdetSK - 1 / 2 * t(yK) %*% SKm1 %*% yK
    SnK <- S[-K, -K]
    SnKK <- S[-K, K[-idxK]]
    SKnK <- t(SnKK)
    muCondK <- -Gamma[-K, i] / 2 + SnKK %*% SKm1 %*% yK
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
        q = c(log(x[-K] / x[i]) + Gamma[-K, i] / 2),
        sd = sqrt(S[-K, -K])
      ))
      names(logdvnK) <- "upper"
    } else {
      logdvnK <- log(mvtnorm::pmvnorm(
        upper = c(log(x[-K] / x[i]) + Gamma[-K, i] / 2),
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
#' @param data Numeric \nxd matrix, containing
#' observations following a multivariate HR Pareto distribution.
#' @param Gamma Numeric \dxd matrix, representing a variogram matrix \eGamma.
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
    return(-n * log(V_HR(x = rep(1, times = d), Gamma = Gamma))
           + sum(logdV_HR(x = data, Gamma = Gamma)))
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
#' @param x Numeric \nxd matrix.
#' @param p Numeric vector with `d` elements.
#'
#' @return Numeric \nxd matrix.
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
      data = data[, x],
      init = G_emp[x[1], x[2]],
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

