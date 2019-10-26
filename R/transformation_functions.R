#' fullGamma
#'
#' Given a \code{graph} and \code{Gamma} matrix with entries only on the
#' edges/cliques of the \code{graph}, it returns the full Gamma matrix
#' implied by the conditional independencies.
#'
#' @param graph Graph object from \code{igraph} package.
#' An undirected block graph, i.e., a decomposable
#' graph where the minimal separators of the cliques have size at most one.
#' @param Gamma Numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram, with entries only inside the cliques. Alternatively, can be a
#' vector containing the entries for each edge in the same order as in
#' graph object.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d} representing the completed
#' Gamma matrix.
fullGamma = function(graph, Gamma){

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

  # transform Gamma if needed
  if(is.vector(Gamma)){

    if (length(Gamma) != e){
      stop(paste("The argument Gamma must be a symmetric d x d matrix,",
                 "or a vector with as many entries as the number of edges",
                 "in the graph."))
    }
    G = matrix(0,d,d)
    G[igraph::ends(graph,igraph::E(graph))] = Gamma
    G = G + t(G)
  }  else{
    G = Gamma

    # check that Gamma is d x d

    if (NROW(G) != d | NCOL(G) != d){
      stop(paste("The argument Gamma must be a symmetric d x d matrix,",
                 "or a vector with as many entries as the number of edges",
                 "in the graph."))
    }

    # check that Gamma is symmetric
    if (any(G != t(G))) {
      stop(paste("The argument Gamma must be a symmetric d x d matrix,",
                 "or a vector with as many entries as the number of edges",
                 "in the graph."))
    }

  }

  # computes cliques
  cli = igraph::max_cliques(graph)
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

    if (length(k0) > 1){
      stop("The given graph is not a block graph.")
    }
    G[setdiff(idx1, k0), setdiff(idx2, k0)] = matrix(rep(G[setdiff(idx1, k0),k0], times=l2-1), l1-1, l2-1) +
      t(matrix(rep(G[setdiff(idx2, k0),k0], times=l1-1), l2-1, l1-1))
    G[setdiff(idx2, k0), setdiff(idx1, k0)] = t(G[setdiff(idx1, k0), setdiff(idx2, k0)])
    cli.selected = c(cli.selected, cli.idx)
    idx1 = union(idx1, idx2)
  }
  return(G)
}



#' Gamma2Graph
#'
#' Transforms Gamma matrix to an \code{igraph} object plots it (optionally).
#'
#' @param Gamma Numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram.
#' @param to_plot Boolean. If \code{TRUE} (default), it plots the produced
#' graph.
#'
#' @return Graph object from \code{igraph} package. An undirected graph.
#'
Gamma2Graph <- function(Gamma, to_plot = T){
  null.mat <- matrix(0, nrow=nrow(Gamma), ncol=ncol(Gamma))
  for(i in 1:nrow(Gamma)){
    null.mat[-i,-i] <- null.mat[-i,-i] +
      (abs(solve(Gamma2Sigma(Gamma, i))) < 1e-6)
  }
  graph = igraph::graph_from_adjacency_matrix(null.mat==0, diag =FALSE,
                                              mode="undirected")
  igraph::V(graph)$color <- "cyan2"
  igraph::V(graph)$size <- 15
  igraph::E(graph)$width <- 2
  igraph::E(graph)$color <- "darkgrey"
  if (to_plot){
    igraph::plot.igraph(graph)
  }
  return(graph)
}



#' data2mpareto
#'
#' Transfroms the data empirically to multivariate Pareto scale.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}.
#' @param p Numeric between 0 and 1. Probability used for the quantile to
#' threshold the data.
#'
#' @return Numeric matrix \eqn{m \times d}{m x d}, where \eqn{m} is the number
#' of rows in the original \code{data} matrix that are above the threshold.
data2mpareto <- function(data, p){
  xx <- 1/(1-apply(data, 2, unif))
  q <- 1 / (1 - p)
  idx <- which(apply(xx, 1, max) > q)
  return(xx[idx,] / q)
}



#' Sigma2Gamma
#'
#' Transforms \eqn{\Sigma^(k)} to the respective Gamma matrix.
#'
#'  @param S Numeric matrix \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)}.
#'  It represents the \eqn{\Sigma^(k)} matrix as defined in equation (10) in
#'  the paper of Engelke, S., and Hitz, A.,
#'  \url{https://arxiv.org/abs/1812.01734}.
#'  @param k Integer between \code{1} (the default value) and \code{d}.
#'  It represents the index that is missing in \eqn{\Sigma^(k)}.
#'  @param full Boolean. If true, then \code{S} must be a
#'  \eqn{d \times d}{d x d} matrix. By default, \code{full = FALSE}.
#'
#'  @return Numeric matrix \eqn{d\times d}{d x d}. It represents a variogram
#'  matrix.
Sigma2Gamma <- function(S, k = 1, full = FALSE){
  # complete S
  if (!full){
    d <- NROW(S)
    S_full <- rbind(rep(0, d + 1), cbind(rep(0, d), S))

    shuffle <- 1:(d + 1)
    shuffle[shuffle <= k] <- shuffle[shuffle <= k] - 1
    shuffle[1] <- k
    shuffle <- order(shuffle)

    S_full <- S_full[shuffle, shuffle]
  } else {
    S_full <- S
  }

  # compute Gamma
  One <- rep(1, times=ncol(S_full))
  D <- diag(S_full)
  Gamma <- One%*%t(D) + D%*%t(One) - 2*S_full

  return(Gamma)
}



#' Gamma2Sigma
#'
#' Transforms Gamma matrix to \eqn{\Sigma^(k)} matrix.
#'
#' @param Gamma Numeric matrix \eqn{d\times d}{d x d}. It represents a variogram
#' matrix.
#' @param k Integer between \code{1} (the default value) and \code{d}.
#' It represents the index that is missing in \eqn{\Sigma^(k)}.
#' @param full Boolean. If true, then \code{S} must be a
#' \eqn{d \times d}{d x d} matrix. By default, \code{full = FALSE}.
#'
#' @return S Numeric matrix \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)}.
#' It represents the \eqn{\Sigma^(k)} matrix as defined in equation (10) in
#' the paper of Engelke, S., and Hitz, A.,
#' \url{https://arxiv.org/abs/1812.01734}.
#'
Gamma2Sigma <- function(Gamma,k=1,full=FALSE){
  d <- ncol(Gamma)
  if(full)
    1/2 * (matrix(rep(Gamma[,k],d), ncol=d,nrow=d) +
             t(matrix(rep(Gamma[,k],d), ncol=d,nrow=d)) - Gamma)
  else
    1/2 * (matrix(rep(Gamma[-k,k],d-1), ncol=d-1,nrow=d-1) +
             t(matrix(rep(Gamma[-k,k],d-1), ncol=d-1,nrow=d-1)) - Gamma[-k,-k])
}








