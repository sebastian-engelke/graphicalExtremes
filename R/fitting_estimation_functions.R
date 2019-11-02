#' Estimate \eqn{\chi}
#'
#' Estimates the chi coefficient empirically.
#'
#' The \eqn{\chi} coefficient is a scalar coefficient that represents
#' the extremal correlation between two variables.
#'
#' @param data Numeric matrix \eqn{n \times d}{n x d}. A data matrix
#' containing \eqn{d} variables.
#' @param p Numeric, between 0 and 1. It is the probability threshold used to
#' compute the \eqn{\chi} coefficient.
#' @param pot Boolean. if TRUE, then pot-type estimation of EC is used.
#' By default, \code{pot = FALSE}.
#'
#' @return Numeric. The empirical \eqn{\chi} coefficient between the \eqn{d}
#' variables in \code{data}.
#'
emp_chi <- function(data, p, pot=FALSE){
  if (!is.matrix(data)){
    stop("The data should be a matrix")
  }

  if (ncol(data) <= 1){
    stop("The data should be a matrix with at least two columns.")
  }

  data <- na.omit(data)
  n <- nrow(data)
  data <-  apply(data, 2, unif)
  rowmax <- apply(data, 1, max)
  rowmin <- apply(data, 1, min)
  eps <- .Machine$double.eps^0.5
  qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)

  qlim <- qlim2
  nq <- length(p)
  cu <- cbaru <- numeric(nq)
  for (i in 1:nq) cu[i] <- mean(rowmax < p[i])
  for (i in 1:nq) cbaru[i] <- mean(rowmin > p[i])
  if(pot) chiu <- cbaru / (1-p)
  if(!pot) chiu <- 2 - log(cu)/log(p)
  return(chiu)
}



#' Estimate matrix of \eqn{\chi}
#'
#' Estimates empirically the extremal \eqn{\chi} correlation coefficient
#' given the dataset \code{data}.
#'
#' @param data Numeric matrix \eqn{n \times d}{n x d}. A data matrix
#' containing \eqn{d} variables.
#' @inheritParams emp_chi
#'
#' @return Numeric matrix \eqn{d\times d}{d x d}. The matrix containing the
#' bivariate extremal coefficientes \eqn{\chi_{ij}}, for \eqn{i, j = 1, ..., d}.
#'
emp_chi_mat <- function(data, p, pot=FALSE){
  d <- ncol(data)
  res <- as.matrix(expand.grid(1:d,1:d))
  res <- res[res[,1]>res[,2],,drop=FALSE]
  chi <- apply(res, 1, function(x){
    emp_chi(cbind(data[,x[1]], data[,x[2]]), p=p, pot=pot)
  })
  chi.mat <- matrix(NA, ncol=d, nrow=d)
  chi.mat[res] <- chi
  chi.mat[res[,2:1]] <- chi
  diag(chi.mat) <- 1

  return(chi.mat)
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
#' the components \code{k}. If \code{NULL},
#' it is assumed that \code{data} has already multivariate Pareto margins.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d}. The estimated
#' variogram of the Huesler-Reiss distribution.
#'
emp_vario <- function(data, k=NULL, p=NULL){
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
logdV_HR <- function(x,par){
  if (any(x <= 0)) {
    stop("The elements of x must be positive.")
  }

  if (is.vector(x)){d <- length(x)}
  if (is.matrix(x)){d <- ncol(x)}

  i <- 1
  G = par2Gamma(par)

  if (NROW(G) != d){
    stop("The length of par must be d * (d - 1) / 2.")
  }

  S <- Gamma2Sigma(G, k=i)
  cholS <- chol(S)
  Sm1 <- chol2inv(cholS)
  logdetS <- 2*sum(log(diag(cholS)))
  if (is.vector(x)){
    y <- (log(x/x[i])+ G[,i]/2)[-i]
    logdv <- - sum(log(x)) - log(x[i]) -((d-1)/2)*log(2*pi) -1/2*logdetS  - 1/2 * t(y)%*%Sm1%*%y
  }
  if (is.matrix(x)){
    y <- (t(t(log(x/x[,i])) + G[,i]/2))[,-i]
    logdv <- - apply(log(x),1,sum) - log(x[,i]) -((d-1)/2)*log(2*pi) -1/2*logdetS  - 1/2 * diag(y%*%Sm1%*%t(y))
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
#' The index set that is \strong{not} censored.
#' @inheritParams V_HR
#'
#' @return Numeric. The censored exponent measure of the HR distribution.
#'
logdVK_HR <- function(x, K, par){

  if (any(x <= 0)) {
    stop("The elements of x must be positive.")
  }

  d <- length(x)
  k <- length(K)
  i <- min(K)
  idxK <- which(K == i)
  G = par2Gamma(par)

  if (NROW(G) != d){
    stop("The length of par must be d * (d - 1) / 2.")
  }

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



#' Full censored log-likelihood of HR model
#'
#' Computes the full (censored) log-likelihood of HR model.
#'
#' @param data Numeric matrix \eqn{n\times d}{n x d}. It contains
#' observations following a multivariate HR Pareto distribution.
#' @param Gamma Numeric matrix \eqn{n\times d}{n x d}.
#' It represents a variogram matrix \eqn{\Gamma}.
#' @param cens Boolean. If TRUE, then censored log-likelihood is computed.
#' By default, \code{cens = FALSE}.
#'
#' @return Numeric. The full censored log-likelihood of HR model.
#'
logLH_HR <- function(data, Gamma, cens = FALSE){
  # helper function
  censor <- function(x,p){
    f2 <- function(x,p){
      y <- c()
      y[x>p] <- x[x>p]
      y[x<=p] <- p[x<=p]
      return(y)
    }
    return(t(apply(x, 1, f2, p)))
  }

  # function body
  if (is.vector(data)){
    d <- length(data)
    n = 1
    data <- t(as.matrix(data))
    # !!! put mistake because n = 1??
  }
  if (is.matrix(data)){
    d <- NCOL(data)
    n = NROW(data)
  }
  par = Gamma2par(Gamma)

  # if cens = FALSE (default)
  if(!cens){
    return(-n*log(V_HR(x=rep(1, times=d),par=par))
           + sum(logdV_HR(x=data,par=par)))
  }

  # if cens = TRUE
  p <- rep(1,d)
  data.p <- censor(data,p)
  r <- nrow(data.p)
  L <- apply(data.p>matrix(p,ncol=d,nrow=r,byrow=TRUE),1,which)
  I <- which(lapply(L,length)>0 & lapply(L,length)<d)
  J <- which(lapply(L,length)==d)

  if (length(I)>0){
    y1 <- mapply(logdVK_HR,x=as.list(data.frame(t(data.p)))[I],K=L[I],
                 MoreArgs=list(par=par))
  } else {
    y1 <- 0
  }

  if (length(J)>0){
    y2 <- logdV_HR(x=data.p[J,],par=par)
  } else {
    y2 <- 0
  }
  return(sum(y1)+sum(y2) - (length(I)+length(J))*log(V_HR(p,par=par)))
}



#' Fit parameters of multivariate HR Pareto distribution
#'
#' This function fits the parameters of a multivariate HR Pareto distribution
#'  using (censored) likelihood estimation.
#'
#'  If \code{graph} is given, it assumes the conditional independence
#'  structure of this graph and fits only the parameters on the edges,
#'  but with the full likelihood. This should only be used for small dimensions.
#'  If you give \code{data}, and not \code{graph}, then the function fits HR
#'  with the whole matrix. This is computationally heavy and works only for
#'  small dimensions, e.g., 3 or 4.
#'
#' @param data Numeric matrix \eqn{n\times d}{n x d}. A dataset containing
#' \eqn{n} observations following a \eqn{d}-variate HR Pareto distribution.
#' @param p Numeric between 0 and 1. Threshold probability. If \code{NULL},
#' it is assumed that \code{data} is already distributed as multivariate Pareto.
#' @param cens Boolean. If TRUE, then censored log-likelihood is computed.
#' By default, \code{cens = FALSE}.
#' @param init Numeric vector. Initial parameter values.
#' @param maxit Positive integer. The maximum number of iterations in the
#' optimization.
#' @param graph Graph object from \code{igraph} package or
#' \code{NULL}.
#' If \code{graph} is not \code{NULL}, then only edge parameters are fit.
#' @param method String. A valid optimization method used by the function
#' \code{\link[stats]{optim}}. By default, \code{method = "BFGS"}.
#'
#' @return List. The list made of:
#' \itemize{
#' \item \code{convergence} Boolean. Whether the optimization converged or not.
#' \item \code{par} Numeric vector. Optimized parameters.
#' \item \code{Gamma} Numeric matrix \eqn{d \time d}{d x d}. Fitted variogram
#' matrix.
#' \item \code{nllik} Numeric. Optimum value of the likelihood function.
#' \item \code{hessian} Numeric matrix. Estimated Hessian matrix of the
#' estimated parameters.
#' }
# !!! check if graph has specific form?
fmpareto_HR <- function(data,
                       p = NULL,
                       cens = FALSE,
                       init,
                       maxit = 100,
                       graph = NULL,
                       method = "BFGS"){

  if(!is.null(p)){
    # if p provided -> data not Pareto -> to convert
    data.std = data2mpareto(data, p)
  } else {
    # if p not provided -> data already Pareto
    data.std <- data
  }

  # censoring at 1 since data already normalized
  p = 1
  d <- ncol(data)
  if (length(p)==1){p <- rep(p,d)}

  # negative log likelihood function
  if(cens){
    # censor below the (multivariate) threshold
    censor <- function(x,p){
      f2 <- function(x,p){
        y <- c()
        y[x>p] <- x[x>p]
        y[x<=p] <- p[x<=p]
        return(y)
      }
      return(t(apply(x,1,f2,p)))
    }
    data.p <- censor(data,p)
    r <- nrow(data.p)

    L <- apply(data.p>matrix(p,ncol=d,nrow=r,byrow=TRUE),1,which)
    I <- which(lapply(L,length)>0 & lapply(L,length)<d)
    J <- which(lapply(L,length)==d)

    nllik <- function(par){
      if(!is.null(graph)){
        Gtmp = complete_Gamma(graph = graph, Gamma = par)
        par = Gtmp[upper.tri(Gtmp)]
      }

      G = par2Gamma(par)
      S <- Gamma2Sigma(G, k=1)

      if (any(par <= 0) | !matrixcalc::is.positive.definite(S)){return(10^50)}

      else {
        if (length(I)>0){y1 <- mapply(logdVK_HR,
                                      x=as.list(data.frame(t(data.p)))[I],
                                      K=L[I],MoreArgs=list(par=par))}
        else {y1 <- 0}
        if (length(J)>0){y2 <- logdV_HR(x=data.p[J,],par=par)}
        else {y2 <- 0}
        y <- sum(y1)+sum(y2) - (length(I)+length(J))*log(V_HR(p,par=par))
        return(-y)
      }
    }
  }
  else{
    r <- nrow(data)
    L <- apply(data>matrix(p,ncol=d,nrow=r,byrow=TRUE),1,which)
    I <- which(lapply(L,length)>0) #1:r
    nllik <- function(par){
      if(!is.null(graph)){
        Gtmp = complete_Gamma(graph = graph, Gamma = par)
        par = Gamma2par(Gtmp)
      }

      G = par2Gamma(par)
      S <- Gamma2Sigma(G, k=1)

      if (any(par <= 0) | !matrixcalc::is.positive.definite(S)){return(10^50)}
      else {
        if (length(I)>0){y1 <- logdV_HR(x=data[I,],par=par)}
        else {y1 <- 0}
        y <- sum(y1) - length(I)*log(V_HR(p,par=par))
        return(-y)
      }
    }
  }

  # optimize likelihood
  opt <- optim(init, nllik, hessian = TRUE,
               control = list(maxit = maxit), method = method)

  z <- list()
  z$convergence <- opt$convergence
  z$par <- opt$par
  if(is.null(graph)) z$Gamma <- par2Gamma(z$par)
  else z$Gamma <- complete_Gamma(graph = graph, Gamma = z$par)
  z$nllik <- opt$value
  z$hessian <- opt$hessian
  return(z)
}



#' Maximum likelihood parameters of a HR block graph.
#'
#' Estimates the parameters of a HR block graph by maximizing the
#' (censored) log-likelihood.
#'
#' @inheritParams fmpareto_HR
#' @param graph Graph object from \code{igraph} package.
#' An undirected block graph, i.e., a decomposable, connected
#' graph where the minimal separators of the cliques have size at most one.
#' @param edges_to_add Numeric matrix \eqn{m\times 2}{m x 2}, where \eqn{m} is
#' the number of edges that are tried to be added in the forward selection.
#' By default, \code{edges_to_add = NULL}.
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{graph} --- Graph object from \code{igraph} package. It represents
#' the block graph passed to the function.
#' \item \code{Gamma} --- Numeric matrix \eqn{d\times d}{d x d}. It represents
#' the estimated variogram matrix \eqn{\Gamma}.
#' }
#'
fmpareto_graph_HR = function(graph, data, p = NULL, cens = FALSE, edges_to_add = NULL){

  # set up main variables
  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # check if it is directed
  if (igraph::is_directed(graph)){
    warning("The given graph is directed. Converted to undirected.")
    graph <- igraph::as.undirected(graph)
  }

  # check if it is connected
  is_connected <- igraph::is_connected(graph)

  if (!is_connected){
    stop("The given graph is not connected.")
  }

  # check if graph is decomposable
  is_decomposable <- igraph::is_chordal(graph)$chordal
  if (!is_decomposable){
    stop("The given graph is not decomposable (i.e., chordal).")
  }

  # check if it is block graph
  cli = igraph::max_cliques(graph)
  ncli = length(cli)
  min_sep <- 0

  for (i in 1:ncli){
    cli1 <- cli[[i]]
    for (j in 1:ncli){
      if (j <= i) {next}
      cli2 <- cli[[j]]

      min_sep <- max(min_sep, length(intersect(cli1, cli2)))

      if (min_sep > 1) {break}
    }
  }

  if (min_sep > 1){
    stop("The given graph is not a block graph.")
  }

  # check if the number of nodes in the graph matches the number
  # of variables in the data matrix
  nnodes = igraph::vcount(graph)
  if (nnodes != NCOL(data)){
    stop(paste("The number of nodes in the graph doesn't match with the number",
               "of variables (i.e., columns) in the data matrix"))
  }

  # check if you need to rescale data or not
  if(!is.null(p)){
    data.std = data2mpareto(data, p)
  } else {
    data.std <- data
  }

  l = 1

  graph.cur = list()
  graph.cur[[l]] = graph
  Ghat = list()
  Ghat[[l]] = matrix(NA, nrow=nnodes, ncol=nnodes)

  # loop through all cliques
  for(i in 1:ncli){
    # pick the curren cliques
    cli.idx = cli[[i]]
    # how many nodes in the current cliques?
    cli.len = length(cli.idx)
    # compute marginal pareto, on the nodes of the current clique
    data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)

    G.est <- emp_vario(data = data.cli)
    init = Gamma2par(G.est)
    Ghat[[l]][cli.idx, cli.idx] = fmpareto_HR(data=data.cli,
                                             init=init, cens=cens)$Gamma

  }

  Ghat[[l]] = complete_Gamma(graph=graph.cur[[l]], Gamma=Ghat[[l]])

  # if you want to add some edges
  if(!is.null(edges_to_add)){
    stop.flag = FALSE
    AIC = 2*ecount(graph.cur[[l]]) - 2 * logLH_HR(data=data.std,
                                                  Gamma = Ghat[[l]], cens=cens)
    added.edges = c()

    while(length(edges_to_add)!=0 & stop.flag==FALSE){
      if(is.vector(edges_to_add)) edges_to_add = t(as.matrix(edges_to_add))
      m = nrow(edges_to_add)
      AIC.tmp = rep(NA,times=m)
      Ghat.tmp = list()

      # go through proposed edges one after the other while retaining a block
      # graph
      # m number of proposed edges
      for(k in 1:m){
        # current temporary graph
        Ghat.tmp[[k]] = Ghat[[l]]
        # add the current proposed edge to the graph
        graph.tmp = igraph::add_edges(graph = graph.cur[[l]],
                                      edges = edges_to_add[k,])

        # if the obtained graph is decomposable
        if(is_chordal(graph.tmp)$chordal){
          # find list of max cliques
          cli = max_cliques(graph.tmp)
          # find in which clique the new proposed edge is. It can be in at most
          # one clique, otherwise, the original graph were not decomposable.
          intersections <-
            sapply(cli, FUN=function(x) length(intersect(x, edges_to_add[k,]))==2)
          ii <-  which(intersections == TRUE)


          # only in the clique itself the separator can be of size > 1
          if(sum(sapply(cli, FUN=function(x)
            length(intersect(x, cli[[ii]])) > 1))==1){
            cat("\nTry edge", edges_to_add[k,])
            cli.idx = cli[[ii]]
            cli.len = length(cli.idx)
            data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)

            G.est <- emp_vario(data = data.cli)
            init = Gamma2par(G.est)
            Ghat.tmp[[k]][cli.idx, cli.idx] = fmpareto_HR(data=data.cli,
                                                         init=init, cens=cens)$Gamma
            Ghat.tmp[[k]] = complete_Gamma(graph=graph.tmp, Gamma=Ghat.tmp[[k]])
            AIC.tmp[k] = 2*igraph::ecount(graph.tmp) -
              2 * logLH_HR(data = data.std, Gamma = Ghat.tmp[[k]], cens=cens)
          }
        }
      }
      if(!all(is.na(AIC.tmp))){
        add.idx = which(AIC.tmp == min(AIC.tmp, na.rm = TRUE))
        cat("\nAdded edge ", edges_to_add[add.idx,])
        l = l+1
        graph.cur[[l]] =
          add_edges(graph = graph.cur[[l-1]], edges = edges_to_add[add.idx,])
        graph.cur[[l]] <- set_graph_parameters(graph.cur[[l]])
        Ghat[[l]] = Ghat.tmp[[add.idx]]
        AIC = c(AIC, AIC.tmp[add.idx])
        added.edges = rbind(added.edges, t(as.matrix(edges_to_add[add.idx,])))
        edges_to_add = edges_to_add[-add.idx,]
      }
      if(all(is.na(AIC.tmp))) stop.flag=TRUE
    }
    return(list(graph=graph.cur,
                Gamma=Ghat, AIC=AIC, added.edges=added.edges))
  }

  return(list(graph=set_graph_parameters(graph), Gamma=Ghat[[1]]))
}



#' Estimate the HR minimum spanning tree
#'
#' Estimates the minimum spanning tree as the HR tree that maximizes
#' the (censored) log-likelihood.
#'
#' @inheritParams logLH_HR
#' @inheritParams fmpareto_HR
#'
#' @return Graph object from \code{igraph} package. A tree.
#'
mst_HR = function(data, p = NULL, cens = FALSE){

  # check if you need to rescale data or not
  if(!is.null(p)){
    data.std = data2mpareto(data, p)
  } else {
    data.std <- data
  }

  n <- nrow(data)
  d = ncol(data)
  graph.full <- make_full_graph(d)
  G.emp = emp_vario(data=data)
  res <- as.matrix(expand.grid(1:d,1:d))
  res <- res[res[,1]>res[,2],,drop=FALSE]
  if(cens)
    bivLLH <- apply(res[,1:2], 1, function(x) -(fmpareto_HR(data=data[,x], init=G.emp[x[1],x[2]], cens=cens)$nllik -
                                                  2*(sum(log(data[which(data[,x[1]] > 1),x[1]])) + sum(log(data[which(data[,x[2]] > 1),x[2]])))))
  if(!cens)
    bivLLH <- apply(res[,1:2], 1, function(x) {par.est = fmpareto_HR(data=data[,x], init=G.emp[x[1],x[2]], cens=cens)$par;
    logLH_HR(data=data[,x], Gamma=par2Gamma(par.est)) + 2*(sum(log(data[,x[1]])) + sum(log(data[,x[2]]))) })

  bivLLH.mat <- matrix(NA, ncol=d, nrow=d)
  bivLLH.mat[res] <- bivLLH
  bivLLH.mat[res[,2:1]] <- bivLLH
  diag(bivLLH.mat) <- 0
  mst.tree = igraph::mst(graph=graph.full, weights = -bivLLH.mat[ends(graph.full,E(graph.full))], algorithm = "prim")

  # set graphical parameters
  mst.tree <- set_graph_parameters(mst.tree)

  # return tree
  return(mst.tree)
}
