
fmpareto_HR_MLE_Theta <- function(
  data,
  p = NULL,
  cens = FALSE,
  init = NULL,
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
  d <- ncol(data)
  n <- nrow(data)
  oneVec <- rep(1, d)

  # no graph => complete graph (does not incur much computational cost here)
  if(is.null(graph)){
    graph <- igraph::make_full_graph(d)
  }
  
  # Helper function to convert parameter vector to Theta matrix
  edgeIndices <- getEdgeIndices(graph, 'upper')
  parToTheta <- function(par){
    Theta <- matrix(0, d, d)
    Theta[edgeIndices] <- par
    Theta <- Theta + t(Theta)
    diag(Theta) <- -rowSums(Theta)
    return(Theta)
  }

  # use emp_vario if no init provided
  if(is.null(init)){
    G0 <- emp_vario(data)
    Theta0 <- ensure_symmetry(Gamma2Theta(G0))
    init <- getEdgeEntries(Theta0, graph, type = 'upper')
  }

  # convert vector of fixed parameters to logical if necessary
  if(!is.logical(fixParams)){
    fixParams <- seq_along(init) %in% fixParams
  }
  
  if(cens){
    # TODO!
    stop('Censoring not implemented yet!')
  } else{
    nllik <- function(par){
      # Combine par with fixed parameters
      if(any(fixParams)){
        par_full <- init
        par_full[!fixParams] <- par
        par <- par_full
      }

      # Complete parameters according to graph structure
      Theta <- parToTheta(par)
      
      # Check if parameters are valid
      if(!is_valid_Theta(Theta)){
        return(10^50)
      }
      
      # Compute likelihood
      G <- Theta2Gamma(Theta)
      y1 <- logdV_HR(x = data, par = G)
      y <- sum(y1) - n * log(V_HR(oneVec, par = G))
      return(-y)
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

  # Prepare return values
  par <- init
  par[!fixParams] <- opt$par
  
  Theta <- parToTheta(par)
  Gamma <- Theta2Gamma(Theta)

  ret <- list(
    convergence = opt$convergence,
    par = par,
    par_opt = opt$par,
    Theta = Theta,
    Gamma = Gamma,
    nllik = opt$value,
    hessian = opt$hessian
  )
  return(ret)
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
fmpareto_HR_MLE_Gamma <- function(
  data,
  p = NULL,
  cens = FALSE,
  init = NULL,
  fixParams = integer(0),
  useTheta = FALSE,
  maxit = 100,
  graph = NULL,
  method = "BFGS"
){
  # if p provided -> data not Pareto -> to convert
  if(!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  # censoring at 1 since data already normalized
  d <- ncol(data)
  n <- nrow(data)
  oneVec <- rep(1, d)

  # use emp_vario if no init provided
  if(is.null(init)){
    Gamma0 <- emp_vario(data)
    if(useTheta){
      Theta0 <- Gamma2Theta(Gamma0)
      init <- getEdgeEntries(Theta0, graph, type = 'upper')
    } else{
      init <- getEdgeEntries(Gamma0, graph, type = 'upper')
    }
  }

  # Prepare helper function to convert (partial) params to Gamma/Theta:
  parToMatrices <- parToMatricesFactory(
    d = d,
    init = init,
    fixParams = fixParams,
    parIsTheta = useTheta,
    graph = graph,
    checkValidity = TRUE
  )

  # Create negative log likelihood function, maybe using censoring
  if(cens) {
    # Since the data is standardized, censor below (1, 1, ...)
    data_cens <- censor(data, oneVec)

    # Check for each entry (of each observation) if it is censored
    censored_entries <- (data_cens <= 1)
    
    # Make sure no observation is completely censored
    obs_completely_censored <- apply(censored_entries, 1, all)
    if(any(obs_completely_censored)){
      stop('Make sure the data is properly standardized (i.e. Inf-norm > 1)!')
    }

    # Get indices of (non-)censored observations
    obs_censored <- apply(censored_entries, 1, any)
    obs_not_censored <- apply(!censored_entries, 1, all)

    nllik <- function(par) {
      # Convert to Gamma/Theta.
      # - Is NULL if par implies an invalid matrix.
      # - Guaranteed to contain Gamma, might contain Theta (if useTheta==TRUE).
      matrices <- parToMatrices(par, forceGamma = TRUE)
      if(is.null(matrices)){
        return(10^50)
      }
      
      ## Compute likelihood
      # Compute censored densities
      y_censored <- vapply(which(obs_censored), FUN.VALUE = 0, function(i){
        logdVK_HR(
          x = data_cens[i,],
          K = which(!censored_entries[i,]),
          Gamma = matrices$Gamma
        )
      })

      # Compute uncensored densities (faster than using `logdVK_HR`)
      y_not_censored <- logdV_HR(
        x = data_cens[obs_not_censored, , drop=FALSE],
        Gamma = matrices$Gamma,
        Theta = matrices$Theta
      )
      
      # Compute combined likelihood
      y <- sum(y_censored) + sum(y_not_censored) - n * log(V_HR(oneVec, par))
      return(-y)
    }
  } else {
    nllik <- function(par) {
      # Convert to Gamma/Theta.
      # - Is NULL if par implies an invalid matrix.
      # - Otherwise contains (at least) one of Gamma, Theta (depends on useTheta)
      matrices <- parToMatrices(par)
      if(is.null(matrices)){
        return(10^50)
      }

      # Compute likelihood
      y1 <- logdV_HR(
        x = data,
        Gamma = matrices$Gamma,
        Theta = matrices$Theta
      )
      y <- sum(y1) - n * log(V_HR(oneVec, par = Gamma))
      return(-y)
    }
  }

  # Actual optimization
  init_opt <- init[!fixParams]
  opt <- stats::optim(
    init_opt,
    nllik,
    hessian = TRUE,
    control = list(maxit = maxit),
    method = method
  )

  # Interpret results
  par <- fillFixedParams(opt$par, init, fixParams)
  matrices <- parToMatrices(opt$par, forceGamma = TRUE, forceTheta = TRUE)
  convergence <- opt$convergence && !is.null(matrices)

  ret <- list(
    convergence = convergence,
    Gamma = matrices$Gamma,
    Theta = matrices$Theta,
    par = par,
    par_opt = opt$par,
    nllik = opt$value,
    hessian = opt$hessian
  )
  return(ret)
}


#' Helper function to combine par with fixed params (in init)
fillFixedParams <- function(par, init, fixParams){
  if(!is.logical(fixParams)){
    fixParams <- seq_along(init) %in% fixParams
  }
  init[!fixParams] <- par
  return(init)
}


# Creates a helper function to convert parameter vector to Gamma/Theta matrix
parToMatricesFactory <- function(
  d,
  init = NULL,
  fixParams = integer(0),
  parIsTheta = FALSE,
  defaultForceTheta = FALSE,
  defaultForceGamma = FALSE,
  graph = NULL,
  checkValidity = TRUE
){
  # Ignore graph if it's the complete graph
  if(igraph::ecount(graph) == d*(d-1)/2){
    graph <- NULL
  }

  # Make sure fixParams is boolean
  if(!is.logical(fixParams)){
    fixParams <- seq_along(init) %in% fixParams
  }

  # Get indices of par in the matrix (according to edges in the graph)
  if(is.null(graph)){
    edgeIndices <- which(upper.tri(matrix(NA, d, d)))
  } else{
    edgeIndices <- getEdgeIndices(graph, 'upper')
  }
  transposedEdgeIndices <- getTransposedIndices(edgeIndices)
  
  # Create parToMatrices(), depending on whether par represents Theta or Gamma
  if(parIsTheta){
    parToMatrices <- function(
      par,
      forceGamma = defaultForceGamma,
      forceTheta = defaultForceTheta
    ){
      # Fill fixed params
      par <- fillFixedParams(par, init, fixParams)
      
      ## Make Theta
      Theta <- matrix(0, d, d)
      Theta[edgeIndices] <- par
      Theta[transposedEdgeIndices] <- par
      diag(Theta) <- (-1)*rowSums(Theta)
      
      # Return NULL if par implies an invalid Theta
      if(checkValidity && !is_valid_Theta(Theta)){
        return(NULL)
      }

      # Compute Gamma if specified
      if(forceGamma){
        Gamma <- Theta2Gamma(Theta)
      } else{
        Gamma <- NULL
      }

      return(list(Gamma = Gamma, Theta = Theta))
    }
  } else{
    parToMatrices <- function(
      par,
      forceGamma = defaultForceGamma,
      forceTheta = defaultForceTheta
    ){
      # Fill fixed params
      par <- fillFixedParams(par, init, fixParams)
      
      # Return NULL if par is non-positive (cheap check before making Gamma)
      if(checkValidity && any(par <= 0)){
        return(NULL)
      }

      # Compute (partial) Gamma:
      Gamma <- matrix(NA, d, d)    
      Gamma[edgeIndices] <- par
      Gamma[transposedEdgeIndices] <- par
      diag(Gamma) <- 0

      # Complete according go graph:
      if(!is.null(graph)){
        Gamma <- complete_Gamma(Gamma, graph)
      }
      
      # Return NULL if par implies an invalid Gamma
      if(checkValidity && !is_sym_cnd(Gamma)){
        return(NULL)
      }
      
      # Compute Theta if specified
      if(forceTheta){
        Theta <- Gamma2Theta(Gamma)
      } else{
        Theta <- NULL
      }
      
      return(list(Gamma = Gamma, Theta = Theta))
    }
  }

  return(parToMatrices)
}

