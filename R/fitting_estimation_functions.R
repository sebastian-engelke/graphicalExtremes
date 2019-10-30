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
#' the components \code{k}. If \code{NULL},
#' it is assumed that \code{data} has already multivariate Pareto margins.
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
#' @param data Numeric matrix eqn{n\times d}{n x d}. It contains
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
  censor <- function(x,u){
    f2 <- function(x,u){
      y <- c()
      y[x>u] <- x[x>u]
      y[x<=u] <- u[x<=u]
      return(y)
    }
    return(t(apply(x,1,f2,u)))
  }

  # function body
  if (is.vector(data)){
    d <- length(data)
    n = 1
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
  u <- rep(1,d)
  data.u <- censor(data,u)
  r <- nrow(data.u)
  L <- apply(data.u>matrix(u,ncol=d,nrow=r,byrow=TRUE),1,which)
  I <- which(lapply(L,length)>0 & lapply(L,length)<d)
  J <- which(lapply(L,length)==d)

  if (length(I)>0){
    y1 <- mapply(logdVK_HR,x=as.list(data.frame(t(data.u)))[I],K=L[I],
                 MoreArgs=list(par=par))
  } else {
    y1 <- 0
  }

  if (length(J)>0){
    y2 <- logdV_HR(x=data.u[J,],par=par)
  } else {
    y2 <- 0
  }
  return(sum(y1)+sum(y2) - (length(I)+length(J))*log(V_HR(u,par=par)))
}



#' Maximum likelihood parameters of a HR block graph.
#'
#' Estimates the parameters of a HR block graph by maximizing the
#' (censored) log-likelihood.
#'
#' @param graph Graph object from \code{igraph} package.
#' An undirected block graph, i.e., a decomposable, connected
#' graph where the minimal separators of the cliques have size at most one.
#' @param data Numeric matrix \eqn{n\times d}{n x d}. It contains \eqn{n}
#' observations following a \eqn{d}-variate HR Pareto distribution.
#' @param p Numeric between 0 and 1. Threshold probability. If \code{NULL},
#' it is assumed that \code{data} has already multivariate Pareto margins.
#' @param cens Boolean. If TRUE, then censored log-likelihood is computed.
#' By default, \code{cens = FALSE}.
#' @param sel.edges Numeric matrix \eqn{m\times 2}{m x 2}, where \eqn{m} is
#' the number of edges that are tried to be added in the forward selection.
#' By default, \code{sel.edges = NULL}.
#'
#' @return List. The list is made of:
#' \itemize{
#' \item \code{graph} --- Graph object from \code{igraph} package. It represents
#' the block graph passed to the function.
#' \item \code{Gamma} --- Numeric matrix \eqn{d\times d}{d x d}. It represents
#' the estimated variogram matrix \eqn{\Gamma}.
#' }
#'
estGraph_HR = function(graph, data, p = NULL, cens = FALSE, sel.edges = NULL){

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
  browser()
  cli = igraph::max_cliques(graph)
  ncli = length(cli)
  min_sep <- 0

  for (i in 1:ncli){
    cli1 <- cli[[i]]
    for (j in 1:ncli){
      if (j <= i) {next}
      cat(i, j)
      cli2 <- cli[[j]]

      min_sep <- max(min_sep, length(intersect(cli1, cli2)))

      if (min_sep > 1) {break}
    }
  }

  if (min_sep > 1){
    stop("The given graph is not a block graph.")
  }

  # check if the number of nodes in the graph matches with the number
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

  for(i in 1:ncli){
    #cat("Fitting clique ",i," with nodes ", cli[[i]]," \n")
    cli.idx = cli[[i]]
    cli.len = length(cli.idx)
    data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)
    # data.cli = data[,cli.idx]

    # remove from here
    if(!is.null(q))  quant = quantile(data.cli,q)
    if(!is.null(thr)) quant = thr
    data.thr = data.cli[which(apply(data.cli, 1, max) > quant),]/quant
    # to here
    data.thr <- data2mpareto(data.cli, p = q) # not needed as fpareto_HR accepts
    # raw data + p

    G.est <- vario.est(data = data.cli, p = q)
    row_averages <- rowMeans(sapply(1:cli.len,
                                    FUN=function(i) vario.est(data=data.thr, k=i)))
    G.est <- matrix(row_means,cli.len,cli.len)
    init = Gamma2par(G.est)
    Ghat[[l]][cli.idx, cli.idx] = fpareto_HR(data=data.cli, init=init, cens=cens)$Gamma
  }
  Ghat[[l]] = fullGamma(graph=graph.cur[[l]], Gamma=Ghat[[l]])

  if(!is.null(sel.edges)){
    # maybe data2mpareto?
    d <- ncol(data)
    if(!is.null(p)){
      data.full.thr = data2mpareto(data, p)
    } else {
      data.full.thr <- data
    }

    # replace from here
    if(!is.null(q))  quant = quantile(data,q)
    if(!is.null(thr)) quant = thr
    data.full.thr = data[which(apply(data, 1, max) > quant),]/quant
    # to here


    stop.flag = FALSE
    AIC = 2*ecount(graph.cur[[l]]) - 2 * logLH_HR(data=data.full.thr,
                                                  Gamma = Ghat[[l]], cens=cens)
    added.edges = c()

    while(length(sel.edges)!=0 & stop.flag==FALSE){
      if(is.vector(sel.edges)) sel.edges = t(as.matrix(sel.edges))
      m = nrow(sel.edges)
      AIC.tmp = rep(NA,times=m)
      Ghat.tmp = list()

      for(k in 1:m){
        Ghat.tmp[[k]] = Ghat[[l]]
        graph.tmp = add_edges(graph = graph.cur[[l]], edges = sel.edges[k,])
        if(is_chordal(graph.tmp)$chordal){
          cli = max_cliques(graph.tmp)
          ii = which(sapply(cli, FUN=function(x) length(intersect(x,sel.edges[k,]))==2)==TRUE)
          if(sum(sapply(cli, FUN=function(x) length(intersect(x,cli[[ii]])) > 1))==1){
            cat("Try edge", sel.edges[k,]," \n")
            cli.idx = cli[[ii]]
            cli.len = length(cli.idx)
            data.cli = data[,cli.idx]
            # change as before
            if(!is.null(q))  quant = quantile(data.cli,q)
            if(!is.null(thr)) quant = thr
            data.thr = data.cli[which(apply(data.cli, 1, max) > quant),]/quant
            G.est <- matrix(rowMeans(sapply(1:cli.len, FUN=function(i) vario.est(data=data.thr, k=i))),cli.len,cli.len)
            ###
            G.est <- vario.est(data = data.cli, p = p)
            init = G.est[upper.tri(G.est)]
            Ghat.tmp[[k]][cli.idx, cli.idx] = fpareto_HR(data=data.thr, init=init, cens=cens)$Gamma
            Ghat.tmp[[k]] = fullGamma(graph=graph.tmp, Gamma=Ghat.tmp[[k]])
            AIC.tmp[k] = 2*ecount(graph.tmp) - 2 * logLH_HR(data=data.full.thr, Gamma = Ghat.tmp[[k]], cens=cens)
          }
        }
      }
      if(!all(is.na(AIC.tmp))){
        add.idx = which(AIC.tmp == min(AIC.tmp, na.rm = TRUE))
        cat("Added edge ", sel.edges[add.idx,]," \n")
        l = l+1
        graph.cur[[l]] = add_edges(graph = graph.cur[[l-1]], edges = sel.edges[add.idx,])
        Ghat[[l]] = Ghat.tmp[[add.idx]]
        AIC = c(AIC, AIC.tmp[add.idx])
        added.edges = rbind(added.edges, t(as.matrix(sel.edges[add.idx,])))
        sel.edges = sel.edges[-add.idx,]
      }
      if(all(is.na(AIC.tmp))) stop.flag=TRUE
    }
    return(list(graph=graph.cur, Gamma=Ghat, AIC=AIC, added.edges=added.edges))
  }

  return(list(graph=graph, Gamma=Ghat[[1]]))
}



#' Estimate the HR minimum spanning tree
#'
#' Estimates the minimum spanning tree as the HR tree that maximizes
#' the (censored) log-likelihood.
#'
#' @inheritParams logLH_HR
#'
#' @return Graph object from \code{igraph} package. A tree.
#'
mst_HR = function(data, cens = FALSE){
  n <- nrow(data)
  d = ncol(data)
  graph.full <- make_full_graph(d)
  G.emp = vario.est(data=data)
  res <- as.matrix(expand.grid(1:d,1:d))
  res <- res[res[,1]>res[,2],,drop=FALSE]
  if(cens)
    bivLLH <- apply(res[,1:2], 1, function(x) -(fpareto_HR(data=data[,x], init=G.emp[x[1],x[2]], cens=cens)$nllik -
                                                  2*(sum(log(data[which(data[,x[1]] > 1),x[1]])) + sum(log(data[which(data[,x[2]] > 1),x[2]])))))
  if(!cens)
    bivLLH <- apply(res[,1:2], 1, function(x) {par.est = fpareto_HR(data=data[,x], init=G.emp[x[1],x[2]], cens=cens)$par;
    logLH_HR(data=data[,x], Gamma=par2Gamma(par.est)) + 2*(sum(log(data[,x[1]])) + sum(log(data[,x[2]]))) })

  bivLLH.mat <- matrix(NA, ncol=d, nrow=d)
  bivLLH.mat[res] <- bivLLH
  bivLLH.mat[res[,2:1]] <- bivLLH
  diag(bivLLH.mat) <- 0
  mst.tree = igraph::mst(graph=graph.full, weights = -bivLLH.mat[ends(graph.full,E(graph.full))], algorithm = "prim")
  return(mst.tree)
}
