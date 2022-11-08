
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
