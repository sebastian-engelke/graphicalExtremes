## Graphical Models functions ####

tr <- function(m) sum(diag(m))

unif <- function(x){rank(x)/(length(x)+1)}


### This function simulates exact samples of multivariate Pareto distributions
#model = the parametric model type; one of "HR", "logistic", "neglogistic", "dirichlet"
#d: dimension of the multivariate Pareto distribution
#no.simu: number of simulations
#par: the respective parameter for the model; a dxd variogram matrix for HR

rmpareto <- function(model, d, no.simu=1, par) {

  stopifnot((d==round(d)) & (d>=1))
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))
  stopifnot(model %in% c("HR", "logistic", "neglogistic", "dirichlet"))

  if (model=="HR") {
    stopifnot(is.matrix(par))
    Gamma = par
    stopifnot(nrow(Gamma) == d & ncol(Gamma) == d)
    cov.mat <- sapply(1:d, function(i) sapply(1:d, function(j)
      (Gamma[i,1] + Gamma[j,1] - Gamma[i,j])/2))
    cov.mat <- cov.mat + 1e-3 ##add constant random effect to avoid numerical problems
    chol.mat <- chol(cov.mat)
  } else if (model=="logistic") {
    stopifnot(length(par) == 1 & 1e-12 < par & par < 1 - 1e-12)
    theta = par
  } else if (model=="neglogistic") {
    stopifnot(par > 1e-12)
    theta = par
  } else if (model=="dirichlet") {
    alpha = par
    stopifnot(length(alpha) == d)
    stopifnot(all(alpha>1e-12))
  }

  counter <- 0
  res <- numeric(0)
  n.total <- 0
  while (n.total < no.simu){
    counter <- counter + 1
    shift <- sample(1:d, no.simu, replace=TRUE)
    for(k in 1:d){
      if (model == "HR") {
        trend <- sapply(1:d, function(j) Gamma[j,k]/2)
      }
      n.k <- sum(shift==k)

      if(n.k>0){
        proc <- switch(model,
                       "HR"           = simu_px_HR(no.simu=n.k, idx=k, trend=trend, chol.mat=chol.mat),
                       "logistic"     = simu_px_logistic(no.simu=n.k, idx=k, N=d, theta=theta),
                       "neglogistic"  = simu_px_neglogistic(no.simu=n.k, idx=k, N=d, theta=theta),
                       "dirichlet"    = simu_px_dirichlet(no.simu=n.k, idx=k, N=d, alpha=alpha)
        )
        stopifnot(dim(proc)==c(n.k, d))
        proc <- proc/rowSums(proc) / (1-runif(nrow(proc)))
        idx.sim <- which(apply(proc,1,max) > 1)
        res <- rbind(res, proc[idx.sim,])
        n.total <- nrow(res)
      }
    }
  }
  return(list(res=res[sample(1:nrow(res), no.simu, replace=FALSE),], counter=counter))
}



### Internal: simulates HR extremal functions
simu_px_HR <- function(no.simu=1, idx, trend, chol.mat) {
  stopifnot(length(idx)==1)
  d <- nrow(chol.mat)
  proc <- t(chol.mat)%*%matrix(rnorm(d*no.simu), ncol=no.simu) - trend
  proc <- exp(t(proc) - proc[idx,])
  return(proc)
}


### Internal: simulates logistic extremal functions
simu_px_logistic <- function(no.simu=1, idx, N, theta) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- matrix(1/gamma(1-theta)*(-log(runif(no.simu*N)))^(-theta), nrow=no.simu, ncol=N)
  res[cbind(1:no.simu,idx)] <- 1/gamma(1-theta)*rgamma(no.simu,shape=1-theta)^(-theta)
  return(res/res[cbind(1:no.simu,idx)])
}

### Internal: simulates negative logistic extremal functions
simu_px_neglogistic <- function(no.simu=1, idx, N, theta) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- matrix(rweibull(no.simu*N, shape=theta, scale=1/gamma(1+1/theta)), nrow=no.simu, ncol=N)
  res[cbind(1:no.simu,idx)] <- 1/gamma(1+1/theta)*rgamma(no.simu,shape=1+1/theta)^(1/theta)
  return(res/res[cbind(1:no.simu,idx)])
}


### Internal: simulates Dirichlet extremal functions
simu_px_dirichlet <- function(no.simu, idx, N, alpha) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  shape <- alpha
  shape[idx] <- alpha[idx] + 1
  shape.mat <- matrix(shape, nrow=N, ncol=no.simu)
  rate.mat <- matrix(alpha, nrow=N, ncol=no.simu)
  res <- t(matrix(rgamma(N*no.simu, shape=shape.mat, rate=rate.mat), nrow=N, ncol=no.simu))
  return(res/res[cbind(1:no.simu,idx)])
}


### Internal: simulates Dirichlet mixture extremal functions
simu_px_dirichlet_mix <- function(no.simu, idx, N, weights, alpha, norm.alpha) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  if (length(idx)==1) {
    k <- sample(1:length(weights), no.simu, replace=TRUE, prob=N*weights*norm.alpha[idx,])
  } else {
    k <- sapply(1:no.simu, function(i) sample(1:length(weights), 1, prob=N*weights*norm.alpha[idx[i],]))
  }
  shape.mat <- alpha[,k,drop=FALSE]
  shape.mat[cbind(idx,1:no.simu)] <- shape.mat[cbind(idx,1:no.simu)]+1
  res <- t(matrix(rgamma(N*no.simu, shape=shape.mat), nrow=N, ncol=no.simu))
  return(res/res[cbind(1:no.simu,idx)])
}


# !!! this is used for fitting
### This function selects all possible edges that can be added to graph
### while still remaining in the class of graphs described in the paper
#graph: the initial graph object
selectEdges = function(graph){
  d = vcount(graph)
  sel.edges = matrix(0, nrow=0, ncol=2)
  for(i in 1:(d-1)) for(j in (i+1):d) if(is_chordal(add_edges(graph = graph, edges = c(i,j)))$chordal & length(as.vector(shortest_paths(graph, from=i, to=j)$vpath[[1]])) !=2) sel.edges = rbind(sel.edges,c(i,j))
  return(sel.edges)
}


### Transforms data empirically to multivariate Pareto scale
#data: nxd data matrix
#p: probability that is used for quantile to threshold the data
data2mpareto <- function(data, p){
  xx <- 1/(1-apply(data, 2, unif))
  q <- 1 / (1 - p) # !!! use theoretical quantile from Pareto
  idx <- which(apply(xx, 1, max) > q)
  return(xx[idx,] / q)
}

### Transforms Gamma matrix to graph and plots it
#Gamma: the parameter matrix
Gamma2Graph <- function(Gamma, to_plot = T){
  null.mat <- matrix(0, nrow=nrow(Gamma), ncol=ncol(Gamma))
  for(i in 1:nrow(Gamma)){
    null.mat[-i,-i] <- null.mat[-i,-i] + (abs(solve(Gamma2Sigma(Gamma, i))) < 1e-6)
  }
    graph = igraph::graph_from_adjacency_matrix(null.mat==0, diag =FALSE, mode="undirected")
  igraph::V(graph)$color <- "cyan2"
  igraph::V(graph)$size <- 15
  igraph::E(graph)$width <- 2
  igraph::E(graph)$color <- "darkgrey"
  if (to_plot){
    igraph::plot.igraph(graph)
  }
  return(graph)
}

### Transforms Sigma^1 matrix to Gamma matrix
#S: Sigma^1 matrix
#full: if TRUE then the dxd Sigma^1 matrix must be supplied
Sigma2Gamma <- function(S, full=FALSE){
  # !!! add argument k = index that is missing, with default = 1
  # !!! full = F -> give S d-1 x d-1, we specify which element is missing.
  # !!! First "complete" Sigma -> create GG
  # !!! G is always d x d
  One <- rep(1, times=ncol(S))
  D <- diag(S)
  if(!full)
    Gamma <- cbind(c(0, D), rbind(t(diag(S)),  One%*%t(D) + D%*%t(One) - 2*S))
  if(full)
    Gamma <- One%*%t(D) + D%*%t(One) - 2*S # !!! use this once Sigma is "completed"
  return(Gamma)
}

### Transforms Gamma matrix Sigma^k matrix
#Gamma: Gamma matrix
#k: index of which Sigma^k matrix should be computed
#full: if TRUE then the dxd Sigma^k matrix is returned
Gamma2Sigma <- function(Gamma,k=1,full=FALSE){
  d <- ncol(Gamma)
  if(full)
    1/2 * (matrix(rep(Gamma[,k],d), ncol=d,nrow=d) + t(matrix(rep(Gamma[,k],d), ncol=d,nrow=d)) - Gamma)
  else
    1/2 * (matrix(rep(Gamma[-k,k],d-1), ncol=d-1,nrow=d-1) + t(matrix(rep(Gamma[-k,k],d-1), ncol=d-1,nrow=d-1)) - Gamma[-k,-k])
}

### Transform extremal coefficient into HR parameter
Theta2Gamma <- function(theta) (2*qnorm(theta/2))^2

### Transform HR parameter into extremal coefficient
Gamma2Theta <- function(gamma) 2*pnorm(sqrt(gamma)/2)



### Estimates the chi and chibar coefficients empirically (modified version from evd::chiplot function)
#data: nx2 data matrix
#u: probability threshold
#pot: if TRUE, then pot-type estimation of EC is used
chi.est <- function(data, u, pot=FALSE)
{
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
  chibaru <- (2 * log(1 - u))/log(cbaru) - 1
  return(c(chiu, chibaru))
}


### Estimates empirically the extremal coefficient
#data: nxd data matrix
#u: probability threshold for chi.est
#Gtrue: if supplied then the estimated EC are plotted against the once implied from this HR matrix
#pot: if TRUE, then pot-type estimation of EC is used
est.theta <- function(data, u, Gtrue=NULL, pot=FALSE){
  d <- ncol(data)
  res <- as.matrix(expand.grid(1:d,1:d))
  res <- res[res[,1]>res[,2],,drop=FALSE]
  theta <- 2 - apply(res[,1:2], 1, function(x) chi.est(cbind(data[,x[1]], data[,x[2]]), u=u, pot=pot))[1,]
  theta.mat <- matrix(NA, ncol=d, nrow=d)
  theta.mat[res] <- theta
  theta.mat[res[,2:1]] <- theta
  diag(theta.mat) <- 1
  if(!is.null(Gtrue)){
    plot(Gamma2Theta(Gtrue[res]), theta, xlim=c(1,2), ylim=c(1,2))
    abline(0,1, xlim=c(1,2))
  }
  return(theta.mat)
}

### Computes the theoretical chi coefficient in 3 dimensions
#Gamma: dxd parameter matrix
chi3D = function(Gamma){
  res = 3 - V(x=rep(1, times=2),par= Gamma2par(Gamma[c(1,2),c(1,2)])) -
    V(x=rep(1, times=2),par= Gamma2par(Gamma[c(1,3),c(1,3)])) -
    V(x=rep(1, times=2),par= Gamma2par(Gamma[c(2,3),c(2,3)])) +
    V(x=rep(1, times=3),par= Gamma2par(Gamma))
  return(res)

}

### Estimates empirically the chi coefficient in 3 dimensions
#data: nxd data matrix
#triplets: mx3 matrix with locations out of 1:d to be evaluated
#u: probability threshold for chi.est
#Gtrue: if supplied then the estimated chi are plotted against the once implied from this HR matrix
#pot: if TRUE, then pot-type estimation of chi is used
est.chi3D <- function(data, triplets, u, Gtrue=NULL, pot=FALSE, main=""){
  d <- ncol(data)
  chi <- apply(triplets, 1, function(x) chi.est(data[,x], u=u, pot=pot))[1,]
  if(!is.null(Gtrue)){
    chi.theo = apply(triplets, 1, function(x) chi3D(Gamma=Gtrue[x,x]))

    par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
        mar = c(5,5,4,2) +.1)
    plot(chi, chi.theo, xlim = c(0.1,.9), ylim = c(0.1,.9), main=main, xlab="Fitted Model", ylab="Empirical")
    abline(0,1, xlim=c(1,2))
  }
  return(chi)
}




### This function estimates the variogram of the Huesler-Reiss distribution empirically
#graph: data matrix of observations of a Huesler-Reiss distribution
#k: component that is conditioned to be larger than u
#p: probability threshold for the kth components
vario.est <- function(data, k=NULL, p=NULL){
  d <- ncol(data)
  if(!is.null(p)){
    xx <- 1/(1-apply(data, 2, unif))
    q <- quantile(xx[,1], p)
    data.std = xx[which(apply(xx, 1, max) > q),]/q
  }
  else
    data.std = data

  if(!is.null(k)){
    idx <- which(data.std[,k]>1)
    if(length(idx) > 1)
      G = Sigma2Gamma(cov(log(data.std[idx,])), full=TRUE)
    else{
      G = matrix(0,d,d)
      cat("### no exceedance ###")
    }
  }
  else{
    G.fun = function(i){
      idx = which(data.std[,i]>1)
      if(length(idx) > 1)
        xx = Sigma2Gamma(cov(log(data.std[idx,])), full=TRUE)
      else{
        xx = matrix(0,d,d)
        cat("### no exceedance ###")
      }
      return(xx)
    }
    G = matrix(rowMeans(sapply(1:d, FUN=function(i)  G.fun(i))) ,d,d)
  }
  return(G)
}

### Estimates chi for multivariate Pareto distributions empirically
chi.mpd.est <- function(data, pot=FALSE)
{
  d <- ncol(data)
  chi.mpd <- sapply(1:d, function(i) sapply(1:d, function(j)
    2*mean(data[,i]>1 & data[,j]>1) / (mean(data[,i]>1) + mean(data[,j]>1))))
  return(chi.mpd)
}

### Plots two chi-estimates against each other
plotChi <- function(Chi.emp,
                    Chi.theo,
                    is.con,
                    main="",
                    PDF = FALSE,
                    filename = "")
{

  if(PDF) pdf(filename, width = 7)
  par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
      mar = c(5,5,4,2) +.1)
  upper.idx <- upper.tri(Chi.emp)
  con.idx <- as.logical(upper.idx * is.con)
  plot(Chi.theo[upper.idx], Chi.emp[upper.idx], xlim = c(0.23,1), ylim = c(0.23,1),main=main,xlab="Fitted Model", ylab="Empirical")
  abline(b = 1, a = 0)
  points(Chi.theo[con.idx],Chi.emp[con.idx],col="blue")
  #legend("topleft", col=c("blue","black"), pch=19, legend=c("flow-connected pairs","flow-unconnected pairs"), bty = "n", cex = 1.5, pt.cex = 1.5)
  #print(mean((Chi.theo - Chi.emp)^2))
  if(PDF) dev.off()
}



### This function simulates tree graphical models, either multivariate Pareto or max-stable distributions
#tree: a graph object that must be a tree
#model: either "HR" or "logistic" for HR or logistics tree model, respectively
#method: either "mpareto" or "maxstable"
#loc, scale, shape: if method="maxstable", output is transformed to general GEV margins
#no.simu: number of simulations
#Gamma: parameter matrix if model="HR"
#theta: parameter if model="logsitic"
simu_tree <- function(tree, model, method, no.simu=1, Gamma=NULL, theta=NULL, alpha.mat=NULL, loc=1, scale=1, shape=1) {
  require("igraph")
  adj =  as.matrix(as_adj(tree))
  d <- nrow(adj)
  e <- ecount(tree)
  ends.mat = ends(tree, E(tree))

  stopifnot(model %in% c("logistic", "HR", "dirichlet"))
  stopifnot((d==round(d)) & (d>=1))
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))

  if (length(loc)  ==1) loc   <- rep(loc  , times=d)
  if (length(scale)==1) scale <- rep(scale, times=d)
  if (length(shape)==1) shape <- rep(shape, times=d)
  stopifnot(all(scale>1e-12))

  if (model=="logistic") {
    stopifnot(1e-12 < theta & theta < 1 - 1e-12)
  } else if (model=="HR") {
    par.vec = Gamma[ends.mat]
  } else if (model=="dirichlet") {
    stopifnot(nrow(alpha.mat) == d-1 & ncol(alpha.mat) == 2)
  }

  ## Define a matrix A[[k]] choosing the paths from k to other vertices
  idx.e <- matrix(0, nrow=d, ncol=d)
  idx.e[ends.mat] = 1:e
  idx.e = idx.e + t(idx.e)

  A <- e.start <- e.end <- list() #e.start[[k]][h] gives the index (1 or 2) of the starting node in the h edge in the tree rooted at k

  for (k in 1:d) {
    A[[k]] <- matrix(0, nrow=d, ncol=e)
    e.start[[k]] = e.end[[k]] = numeric(e)
    short.paths <- shortest_paths(tree, from = k, to=1:d)
    for(h in 1:d){
      path = short.paths$vpath[[h]]
      idx.tmp = idx.e[cbind(path[-length(path)], path[-1])]
      A[[k]][h,idx.tmp] <- 1
      e.start[[k]][idx.tmp] = apply(ends.mat[idx.tmp,] == matrix(path[-length(path)], nrow = length(idx.tmp), ncol=2), MARGIN=1, FUN = function(x) which(x==TRUE)) #path[-length(path)]
      e.end[[k]][idx.tmp] = apply(ends.mat[idx.tmp,] == matrix(path[-1], nrow = length(idx.tmp), ncol=2), MARGIN=1, FUN = function(x) which(x==TRUE))  #path[-1]
    }
  }

  if(method=="mpareto")
  {
    counter <- 0
    res <- numeric(0)
    n.total <- 0
    while (n.total < no.simu) {
      counter <- counter + 1
      shift <- sample(1:d, no.simu, replace=TRUE)
      for(k in 1:d){
        n.k <- sum(shift==k)
        if(n.k>0){
          proc <- switch(model,
                         "HR" = simu_px_tree_HR(no.simu=n.k, G.vec=par.vec, A = A[[k]]),
                         "logistic"     = simu_px_tree_logistic(no.simu=n.k, idx=k, nb.edges=e, theta=theta, A=A),
                         "dirichlet"     = simu_px_tree_dirichlet(no.simu=n.k, alpha.start = alpha.mat[cbind(1:e, e.start[[k]])],
                                                                  alpha.end = alpha.mat[cbind(1:e, e.end[[k]])], A=A[[k]])
          )
          stopifnot(dim(proc)==c(n.k, d))
          proc <- proc/rowSums(proc) / (1-runif(nrow(proc)))
          idx.sim <- which(apply(proc,1,max) > 1)
          res <- rbind(res, proc[idx.sim,])
          n.total <- nrow(res)
        }
      }
    }
  }else if(method=="maxstable"){
    counter <- rep(0, times=no.simu)
    res <- matrix(0, nrow=no.simu, ncol=d)
    for (k in 1:d) {
      poisson <- rexp(no.simu)

      while (any(1/poisson > res[,k])) {
        ind <- (1/poisson > res[,k])
        n.ind <- sum(ind)
        idx <- (1:no.simu)[ind]
        counter[ind] <- counter[ind] + 1
        proc <- switch(model,
                       "HR" = simu_px_tree_HR(no.simu=n.ind, G.vec=par.vec, A = A[[k]]),
                       "logistic"     = simu_px_tree_logistic(no.simu=n.ind, idx=k, nb.edges=e, theta=theta, A=A),
                       "dirichlet"     = simu_px_tree_dirichlet(no.simu=n.ind, alpha.start = alpha.mat[cbind(1:e, e.start[[k]])],
                                                                alpha.end = alpha.mat[cbind(1:e, e.end[[k]])], A=A[[k]])
        )
        stopifnot(dim(proc)==c(n.ind, d))
        if (k==1) {
          ind.upd <- rep(TRUE, times=n.ind)
        } else {
          ind.upd <- sapply(1:n.ind, function(i)
            all(1/poisson[idx[i]]*proc[i,1:(k-1)] <= res[idx[i],1:(k-1)]))
        }
        if (any(ind.upd)) {
          idx.upd <- idx[ind.upd]
          res[idx.upd,] <- pmax(res[idx.upd,], 1/poisson[idx.upd]*proc[ind.upd,])
        }
        poisson[ind] <- poisson[ind] + rexp(n.ind)
      }
    }
    res <- sapply(1:d, function(i) {
      if (abs(shape[i]<1e-12)) {
        return(log(res[,i])*scale[i] + loc[i])
      } else {
        return(1/shape[i]*(res[,i]^shape[i]-1)*scale[i] + loc[i])
      }
    })

  }
  return(list(res=res[sample(1:nrow(res), no.simu, replace=FALSE),], counter=counter))
}


### This function takes the parameters (upper triangular Gamma matrix) and returns full Gamma
#par: upper triangular Gamma matrix
par2Gamma = function(par){
  d = 1/2 + sqrt(1/4 + 2*length(par))
  stopifnot(round(d)==d)
  G = matrix(0, nrow=d, ncol=d)
  G[upper.tri(G)] = par
  return(G + t(G))
}

### Inverse of par2Gamma
Gamma2par = function(Gamma){
  if(is.matrix(Gamma))
    return(Gamma[upper.tri(Gamma)])
  else
    return(Gamma)
}



### Exponent measure of HR distribution
#x: vector of dimension d where the exponent measure is to be evaluated
#par: upper triangular Gamma matrix
V <- function(x, par){
  d <- length(x)
  G = par2Gamma(par)
  stopifnot(nrow(G)==d)
  f1 <- function(i,x){
    S <- Gamma2Sigma(G, k=i)
    return(1/x[i]*pmvnorm(upper=(log(x/x[i])+G[,i]/2)[-i],mean=rep(0,d-1),sigma= S)[1])
  }
  return(sum(apply(cbind(1:d),1,f1,x=x)))
}

### Exponent measure density of HR distribution
#x: vector of dimension d or matrix of dimension nxd where the exponent measure density is to be evaluated
#par: upper triangular Gamma matrix
logdV <- function(x,par){
  if (is.vector(x)){d <- length(x)}
  if (is.matrix(x)){d <- ncol(x)}
  i <- 1
  G = par2Gamma(par)
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

### Censored exponent measure density of HR distribution
#x: vector of dimension d where the censored exponent measure density is to be evaluated
#K: the index set that is NOT censored
#par: upper triangular Gamma matrix
logdVK <- function(x,K,par){
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
    logdvnK <- log(pmvnorm(upper=c(log(x[-K]/x[i])-muCondK),sigma=SCondK)[1])
    logdv <- logdvK + logdvnK
  }
  if(k==1){
    logdvK <- - 2*log(x[i])
    logdvnK <- log(pmvnorm(upper=c(log(x[-K]/x[i]) + G[-K,i]/2),sigma=S[-K,-K])[1])
    logdv <- logdvK + logdvnK
  }

  return(logdv)
}

### Full (censored) log-likelihood of HR model
#data: data nxd matrix of observation of multivariate HR Pareto distribution
#Gamma: parameter matrix
#cens: if TRUE then censored log-likelihood is computed
logLH_HR <- function(data,Gamma,cens=FALSE){
  if (is.vector(data)){
    d <- length(data)
    n = 1
  }
  if (is.matrix(data)){
    d <- ncol(data)
    n = nrow(data)
  }
  par = Gamma2par(Gamma)
  if(!cens) return(-n*log(V(x=rep(1, times=d),par=par)) + sum(logdV(x=data,par=par)))
  if(cens){
    u <- rep(1,d)
    censor <- function(x,u){
      f2 <- function(x,u){
        y <- c()
        y[x>u] <- x[x>u]
        y[x<=u] <- u[x<=u]
        return(y)
      }
      return(t(apply(x,1,f2,u)))
    }

    data.u <- censor(data,u)
    r <- nrow(data.u)
    L <- apply(data.u>matrix(u,ncol=d,nrow=r,byrow=TRUE),1,which)
    I <- which(lapply(L,length)>0 & lapply(L,length)<d)
    J <- which(lapply(L,length)==d)

    if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
    else {y1 <- 0}
    if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
    else {y2 <- 0}
    return(sum(y1)+sum(y2) - (length(I)+length(J))*log(V(u,par=par)))
  }
}




### This function fits the parameters of a multivariate HR Pareto distribution using (censored) likelihood estimation.
### If a graph is given, it assumes the cond. independence structure of this graph fits only the parameters on the edges,
##  but with the full likelihood. This should only be used for small dimensions; for graphs use fpareto_graph
#data: data nxd matrix of observation of multivariate HR Pareto distribution
#cens: if TRUE then censored log-likelihood is computed
#init: initial parameter value
#maxit: maximum number of iteration in the optimization
#graph: if not null, then only edge parameters are fit
#method: optimization method
fpareto_HR <- function(data,
                       cens=TRUE,
                       init,
                       maxit = 100,
                       graph=NULL,
                       method = "BFGS"){

  require("mvtnorm")
  ti <- proc.time()
  u = 1  #censoring at 1 since data already normalized
  d <- ncol(data)
  if (length(u)==1){u <- rep(u,d)}

  # negative log likelihood function
  if(cens){
    # censor below the (multivariate) threshold
    censor <- function(x,u){
      f2 <- function(x,u){
        y <- c()
        y[x>u] <- x[x>u]
        y[x<=u] <- u[x<=u]
        return(y)
      }
      return(t(apply(x,1,f2,u)))
    }
    data.u <- censor(data,u)
    r <- nrow(data.u)

    L <- apply(data.u>matrix(u,ncol=d,nrow=r,byrow=TRUE),1,which)
    I <- which(lapply(L,length)>0 & lapply(L,length)<d)
    J <- which(lapply(L,length)==d)

    nllik <- function(par){
      if(!is.null(graph)){
        Gtmp = fullGamma(graph = graph, Gamma = par)
        par = Gtmp[upper.tri(Gtmp)]
      }

      G = par2Gamma(par)
      S <- Gamma2Sigma(G, k=1)

      if (any(par <= 0) | !is.positive.definite(S)){return(10^50)}

      else {
        if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
        else {y1 <- 0}
        if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
        else {y2 <- 0}
        y <- sum(y1)+sum(y2) - (length(I)+length(J))*log(V(u,par=par))
        return(-y)
      }
    }
  }
  else{
    r <- nrow(data)
    L <- apply(data>matrix(u,ncol=d,nrow=r,byrow=TRUE),1,which)
    I <- which(lapply(L,length)>0) #1:r
    nllik <- function(par){
      if(!is.null(graph)){
        Gtmp = fullGamma(graph = graph, Gamma = par)
        par = Gtmp[upper.tri(Gtmp)]
      }

      G = par2Gamma(par)
      S <- Gamma2Sigma(G, k=1)

      if (any(par <= 0) | !is.positive.definite(S)){return(10^50)}
      else {
        if (length(I)>0){y1 <- logdV(x=data[I,],par=par)}
        else {y1 <- 0}
        y <- sum(y1) - length(I)*log(V(u,par=par))
        return(-y)
      }
    }
  }

  # optimize likelihood
  opt <- optim(init,nllik,hessian=TRUE,control = list(maxit = maxit), method = method)

  z <- list()
  z$convergence <- opt$convergence
  z$par <- opt$par
  if(is.null(graph)) z$Gamma <- par2Gamma(z$par)
  else z$Gamma <- fullGamma(graph = graph, Gamma = z$par)
  z$nllik <- opt$value
  z$hessian <- opt$hessian
  z$time <- (proc.time()-ti)[3]
  return(z)
}



# !!! the graph should be decomposable. It is only made for block graphs
# Block graph
### This function takes a graph and Gamma matrix specified only on the edges/cliques
### of this graph and returns the full Gamma matrix implied by the conditional indepdencies
#graph: graph object from igraph package
#Gamma: the Gamma with entries only inside the cliques; or vector with weights for each
#       edge in the same order as in graph object
fullGamma = function(graph, Gamma){ # !!! block_gamma_completion

  if(is.vector(Gamma)){
    G = matrix(0,d,d)
    G[ends(graph,E(graph))] = Gamma
    G = G + t(G)
  }
  else G = Gamma
  cli = max_cliques(graph)
  ncli = length(cli)
  cli.selected = 1
  idx1 = cli[[1]]
  V = 1:ncli

  for(i in 1:(ncli-1)){
    cli.idx = min(V[which(sapply(V, function(j) length(intersect(idx1, cli[[j]])) > 0) == 1 & !is.element(V, cli.selected))])
    idx2 = cli[[cli.idx]]
    l1 = length(idx1)
    l2 = length(idx2)
    k0 = intersect(idx1, idx2)
    G[setdiff(idx1, k0), setdiff(idx2, k0)] = matrix(rep(G[setdiff(idx1, k0),k0], times=l2-1), l1-1, l2-1) +
      t(matrix(rep(G[setdiff(idx2, k0),k0], times=l1-1), l2-1, l1-1))
    G[setdiff(idx2, k0), setdiff(idx1, k0)] = t(G[setdiff(idx1, k0), setdiff(idx2, k0)])
    cli.selected = c(cli.selected, cli.idx)
    idx1 = union(idx1, idx2)
  }
  return(G)
}


### This function estimates the minimum spanning tree as the HR tree that maximizes the (censored) LLH
#data: nxd data matrix
#cens: censored LLH if TRUE
mst_HR = function(data, cens=TRUE){
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

### This function estimates the parameters of a HR block graph maximizes the (censored) LLH
#graph: graph object from igraph package
#data: nxd data matrix, where d is the number of nodes of the graph and n the sample size
#q: the theshold probability, e.g. 0.9
#thr: alternative to q, the absolute threshold
#cens: logical, whether censored estimation is performed
#sel.edges: if provided, then it must be mx2 matrix, where m is the number of edges that are tried to add in forward selection
estGraph_HR = function(graph, data, q=NULL, thr=NULL, cens=TRUE, sel.edges=NULL){
  stopifnot((is.null(q) + is.null(thr)) == 1)
  cli = max_cliques(graph)
  ncli = length(cli)
  nnodes = vcount(graph)
  stopifnot(nnodes==ncol(data))
  l = 1

  graph.cur = list()
  graph.cur[[l]] = graph
  Ghat = list()
  Ghat[[l]] = matrix(NA, nrow=nnodes, ncol=nnodes)

  for(i in 1:ncli){
    #cat("Fitting clique ",i," with nodes ", cli[[i]]," \n")
    cli.idx = cli[[i]]
    cli.len = length(cli.idx)
    data.cli = data[,cli.idx]
    if(!is.null(q))  quant = quantile(data.cli,q)
    if(!is.null(thr)) quant = thr
    data.thr = data.cli[which(apply(data.cli, 1, max) > quant),]/quant
    G.est <- matrix(rowMeans(sapply(1:cli.len, FUN=function(i) vario.est(data=data.thr, k=i))),cli.len,cli.len)
    init = G.est[upper.tri(G.est)]
    Ghat[[l]][cli.idx, cli.idx] = fpareto_HR(data=data.thr, init=init, cens=cens)$Gamma
  }
  Ghat[[l]] = fullGamma(graph=graph.cur[[l]], Gamma=Ghat[[l]])

  if(!is.null(sel.edges)){
    if(!is.null(q))  quant = quantile(data,q)
    if(!is.null(thr)) quant = thr
    data.full.thr = data[which(apply(data, 1, max) > quant),]/quant
    stop.flag = FALSE
    AIC = 2*ecount(graph.cur[[l]]) - 2 * logLH_HR(data=data.full.thr, Gamma = Ghat[[l]], cens=cens)
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
            if(!is.null(q))  quant = quantile(data.cli,q)
            if(!is.null(thr)) quant = thr
            data.thr = data.cli[which(apply(data.cli, 1, max) > quant),]/quant
            G.est <- matrix(rowMeans(sapply(1:cli.len, FUN=function(i) vario.est(data=data.thr, k=i))),cli.len,cli.len)
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




