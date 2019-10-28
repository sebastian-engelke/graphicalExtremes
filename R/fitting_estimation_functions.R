#' Estimate \eqn{\chi}
#'
#' Estimates the chi coefficient empirically.
#'
#' The \eqn{\chi} coefficient is a scalar coefficient that represents
#' the extremal correlation between two variables.
#'
#' @param data Numeric matrix \eqn{n \times 2}{n x 2}. A data matrix
#' containing two variables.
#' @param u Numeric, between 0 and 1. It is the probability threshold used to
#' compute the \eqn{\chi} coefficient.
#' @param pot Boolean. if TRUE, then pot-type estimation of EC is used.
#' By default, \code{pot = FALSE}.
#'
#' @return Numeric. The empirical \eqn{\chi} coefficient between the 2 variables
#' in \code{data}.
#'
chi.est <- function(data, u, pot=FALSE){
  if (NCOL(data) > 2){
    warning("The data matrix should contain only two columns.")
  }

  data <- na.omit(data)
  n <- nrow(data)
  data <-  apply(data, 2, unif)
  rowmax <- apply(data, 1, max)
  rowmin <- apply(data, 1, min)
  eps <- .Machine$double.eps^0.5
  qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)

  qlim <- qlim2
  nq <- length(u)
  cu <- cbaru <- numeric(nq)
  for (i in 1:nq) cu[i] <- mean(rowmax < u[i])
  for (i in 1:nq) cbaru[i] <- mean(rowmin > u[i])
  if(pot) chiu <- cbaru / (1-u)
  if(!pot) chiu <- 2 - log(cu)/log(u)
  return(chiu)
}



#' Estimate matrix of \eqn{\chi}
#'
#' Estimates empirically the extremal \eqn{\chi} correlation coefficient
#' given the dataset \code{data}.
#'
#' @param data Numeric matrix \eqn{n \times d}{n x d}. A data matrix
#' containing \eqn{d} variables.
#' @inheritParams chi.est
#'
#' @return Numeric matrix \eqn{d\times d}{d x d}. The matrix containing the
#' bivariate extremal coefficientes \eqn{\chi_{ij}}, for \eqn{i, j = 1, ..., d}.
#'
chi_mat <- function(data, u, pot=FALSE){
  d <- ncol(data)
  res <- as.matrix(expand.grid(1:d,1:d))
  res <- res[res[,1]>res[,2],,drop=FALSE]
  chi <- apply(res, 1, function(x){
    chi.est(cbind(data[,x[1]], data[,x[2]]), u=u, pot=pot)
  })[1]
  chi.mat <- matrix(NA, ncol=d, nrow=d)
  chi.mat[res] <- chi
  chi.mat[res[,2:1]] <- chi
  diag(chi.mat) <- 1

  return(chi.mat)
}



#' Compute theoretical \eqn{\chi} in 3D
#'
#' Computes the theoretical \eqn{\chi} coefficient in 3 dimensions.
#'
#' @param Gamma Numeric matrix \eqn{3\times 3}{3 x 3}.
#'
#' @return The 3-dimensional \eqn{\chi} coefficient, i.e.,
#' the extremal correlation coefficient for the HR distribution. Note that
#' \eqn{0 \leq \chi \leq 1}.
#'
Gamma2Chi_HR = function(Gamma){
  d <- NROW(Gamma)
  if (d != 3){
    stop("Gamma must be a 3 x 3 matrix.")
  }
  res = 3 - V_HR(x=rep(1, times=2),par= Gamma2par(Gamma[c(1,2),c(1,2)])) -
    V_HR(x=rep(1, times=2),par= Gamma2par(Gamma[c(1,3),c(1,3)])) -
    V_HR(x=rep(1, times=2),par= Gamma2par(Gamma[c(2,3),c(2,3)])) +
    V_HR(x=rep(1, times=3),par= Gamma2par(Gamma))
  return(res)
}



#' Estimate \eqn{\Gamma}
#'
#' Estimates the variogram of the Huesler-Reiss distribution empirically.
#'
#' @param data Numeric matrix \eqn{n\times d}{n x d}. Data matrix of
#' observations following a Huesler-Reiss distribution.
#' @param k Integer between 1 and \eqn{d}. Component of the multivariate
#' observations that is conditioned to be larger than the threshold \code{p}.
#' @param p Numeric between 0 and 1. Probability threshold for the
#' the components \code{k}.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d}. The estimated
#' variogram of the Huesler-Reiss distribution.
#'
vario.est <- function(data, k=NULL, p=NULL){
  # helper ####
  G.fun = function(i, data){
    idx = which(data[,i]>1)
    if(length(idx) > 1)
      xx = Sigma2Gamma(cov(log(data[idx,])), full=TRUE)
    else{
      xx = matrix(0,d,d)
    }
    return(xx)
  }

  # body ####
  d <- ncol(data)
  if(!is.null(p)){
    data.std = data2mpareto(data, p)
  } else {
    data.std <- data
  }

  if(!is.null(k)){

    G <- G.fun(k, data.std)

  } else {

    # take the average
    row_averages <- rowMeans(sapply(1:d, FUN = function(i){
      G.fun(i, data.std)
    }))
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
V_HR <- function(x, par){
  # helper function ####
  f1 <- function(i,x){
    S <- Gamma2Sigma(G, k=i)
    return(1/x[i]*mvtnorm::pmvnorm(upper=(log(x/x[i])+G[,i]/2)[-i],
                                   mean=rep(0,d-1),sigma= S)[1])
  }

  # function body ####
  if (any(x <= 0)) {
  stop("The elements of x must be positive.")
  }

  d <- length(x)
  G = par2Gamma(par)

  if (NROW(G) != d){
    stop("The length of par must be d * (d - 1) / 2.")
  }

  return(sum(apply(cbind(1:d),1,f1,x=x)))
}



#' Compute censored exponent measure
#'
#' Computes the censored exponent measure density of HR distribution.
#'
#' @param x Numeric vector with \eqn{d} positive elements
#' where the censored exponent measure is to be evaluated.
#' @param K Integer vector, subset of \eqn{\{1, \dots, d\}}{{1, ..., d}}.
#' The index set that is \strong{not} censored.
#' @inheritParams V_HR
#'
#' @return Numeric. The censored exponent measure of the HR distribution.
#'
logdVK_HR <- function(x, K, par){
  d <- length(x)
  k <- length(K)
  i <- min(K)
  idxK <- which(K == i)
  G = par2Gamma(par)
  S <- Gamma2Sigma(G, k=i, full=TRUE)
  if(k>1){
    SK <- S[K[-idxK],K[-idxK]]
    cholSK <- chol(SK)
    SKm1 <- chol2inv(cholSK)
    logdetSK <- 2*sum(log(diag(cholSK)))
    idxK <- which(K == i)
    yK <- (log(x[K]/x[i])+ G[K,i]/2)[-idxK]
    logdvK <- - sum(log(x[K])) - log(x[i]) -((k-1)/2)*log(2 *pi) - 1/2*logdetSK - 1/2 * t(yK)%*%SKm1%*%yK
    SnK <- S[-K,-K]
    SnKK <- S[-K,K[-idxK]]
    SKnK <- t(SnKK)
    muCondK <- -G[-K,i]/2 + SnKK %*% SKm1 %*% yK
    if(k < d-1)
      SCondK <- SnK - SnKK %*% SKm1 %*% SKnK
    if(k == d-1)
      SCondK <- SnK - t(SnKK) %*% SKm1 %*% t(SKnK)
    logdvnK <- log(mvtnorm::pmvnorm(upper=c(log(x[-K]/x[i])-muCondK),sigma=SCondK)[1])
    logdv <- logdvK + logdvnK
  }
  if(k==1){
    logdvK <- - 2*log(x[i])
    logdvnK <- log(mvtnorm::pmvnorm(upper=c(log(x[-K]/x[i]) + G[-K,i]/2),sigma=S[-K,-K])[1])
    logdv <- logdvK + logdvnK
  }

  return(logdv)
}
