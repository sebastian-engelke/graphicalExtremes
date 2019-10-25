# Gamma2Graph
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


# data2mpareto
# Sigma2Gamma

#' fullGamma
#'
#' Given a \code{graph} and \code{Gamma} matrix with entries only on the
#' edges/cliques of the \code{graph}, it returns the full Gamma matrix
#' implied by the conditional independencies.
#'
#' @param graph Graph object from \code{igraph} package.
#' An undirected block graph, i.e., a decomposable
#' graph where the minimal separators of the cliques have size at most one.
#' @param Gamma numeric matrix representing a \eqn{d \times d}{d x d}
#' variogram, with entries only inside the cliques. Alternatively, can be a
#' vector containing the entries for each edge in the same order as in
#' graph object.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d} representing the completed
#' Gamma matrix.
fullGamma = function(graph, Gamma){ # !!! change name -> block_gamma_completion

  # check if it is directed
  if (igraph::is_directed(graph)){
    warning("The given graph is directed. Converted to undirected.")
    graph <- igraph::as.undirected(graph)
  }

  # check if it is connected
  # !!!

  # check if graph is decomposable
  is_decomposable <- igraph::is_chordal(graph)$chordal
  if (!is_decomposable){
    stop("The given graph is not decomposable (i.e., chordal).")
  }

  # transform Gamma if needed
  if(is.vector(Gamma)){
    G = matrix(0,d,d)
    G[ends(graph,igraph::E(graph))] = Gamma
    G = G + t(G)
  }
  else{
    G = Gamma
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

# par2Gamma
# Gamma2par


# unif
# selectEdges
# chi.est
# est.theta
# chi3D
# est.chi3D
# vario.est
# chi.mpd.est
# V
# logdV
# logdVK
# logLH_HR
# fpareto_HR
# mst_HR
# estGraph_HR

