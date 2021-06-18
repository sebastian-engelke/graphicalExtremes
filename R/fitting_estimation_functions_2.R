

#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on non-decomposable graphs
#'
#' Works by computing an empirical Gamma matrix and then fitting this to the given `graph`.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param graph Connected graph object from \code{igraph} package.
#' @param ... Further arguments passed to [complete_gamma_general()] if `graph`
#' is not decomposable
#'
#'
#' @export
fmpareto_graph_HR_general <- function(data, graph, p = NULL, ...) {
  # number of data columns (and expected vertices in graph)
  d <- NCOL(data)

  # make sure graph is valid
  graph <- check_graph(graph, graph_type = 'general', nVertices = d)

  # compute the empirical variogram
  Gamma_emp <- emp_vario(data, p = p)

  # fit the empirical variogram to the graph
  Gamma_graph <- complete_Gamma(Gamma_emp, graph, allowed_graph_type = 'general', ...)

  return(list(
    Gamma = Gamma_graph,
    graph = graph
  ))
}

#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on decomposable graphs
#'
#' Similar to [fmpareto_graph_HR()]. Differences are:
#' * Works with any decomposable graph, not just block graphs
#' * Does not support `edges_to_add`
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, \code{cens = FALSE}.
#' @param graph Decomposable graph object from \code{igraph} package.
#'
#' @export
fmpareto_graph_HR_decomposable <- function(data, graph, p = NULL, cens = FALSE) {
  # number of data columns (and expected vertices in graph)
  d <- NCOL(data)

  # make sure graph is valid
  graph <- check_graph(graph, graph_type = 'decomposable', nVertices = d)

  # rescale data if necessary:
  if(!is.null(p)){
    data <- data2mpareto(data, p)
  }

  # compute cliques:
  cliques <- get_cliques_and_separators(graph)

  # initialize variables:
  Ghat <- matrix(NA, d, d)
  fixParams <- logical(d)

  # loop through cliques and estimate Ghat:
  for(clique in cliques){
    # compute marginal pareto, on the nodes of the current clique
    data.cli <- mparetomargins(data = data, set_indices = clique)

    # find Ghat-entries that are already fixed by a previous clique/separator:
    G.cli <- Ghat[clique, clique]
    par.cli <- Gamma2par(G.cli)
    fixParams <- !is.na(par.cli)

    # get initial parameters, keeping the fixed ones:
    G.est0 <- emp_vario(data = data.cli)
    G.est <- replaceGammaSubMatrix(G.est0, G.cli)
    init.cli <- Gamma2par(G.est)

    # estimate parameters
    par.cli <- fmpareto_HR(
      data = data.cli,
      init = init.cli,
      fixParams = fixParams,
      cens = cens
    )$par

    # update Ghat
    G.cli <- par2Gamma(par.cli)
    Ghat[clique, clique] <- G.cli
  }

  # fill entries that don't correspond to edges:
  Ghat <- complete_Gamma(Ghat, allowed_graph_type = 'decomposable')

  return(list(graph = graph, Gamma = Ghat))
}


#' Get Cliques and Separators of a graph
#'
#' Finds all cliques, separators, and (recursively) separators of separators
#' in a graph.
#'
#' @param graph An [igraph::graph] object
#' @return A list of vertex sets that represent the cliques and (recursive)
#' separators of `graph`
get_cliques_and_separators <- function(graph){
  # start with maximal cliques as new cliques:
  newCliques <- igraph::max_cliques(graph)
  cliques <- list()

  while(length(newCliques)>0){
    # compute all separators between the new cliques:
    separators <- list()
    for(i in seq_along(newCliques)){
      for(j in seq_len(i-1)){
        sep <- intersect(newCliques[[i]], newCliques[[j]])
        if(length(sep)>1){
          separators <- c(separators, list(sep))
        }
      }
    }
    # add new cliques to cliques-list:
    # (separators are added after the next iteration)
    cliques <- c(newCliques, cliques)

    # use separators as new cliques:
    newCliques <- setdiff(separators, cliques)
  }

  return(cliques)
}


#' Experimental: Fit graph using empirical Theta matrix
#'
#' Fits a graph to an empirical Gamma matrix by computing the corresponding
#' Theta matrix using [Gamma2Theta()] and greedily chooses `m` edges that
#' correspond to high absolute values in Theta.
#'
#' @param data The (standardized) data from which to compute Gamma
#' @param m The number of edges to add, defaults to the number of edges in a tree
#' @param Gamma_emp The empirical Gamma matrix
#' (can be `NULL` if `data` is given)
#'
#' @return A list containing an [igraph::graph] object and a fitted `Gamma` matrix
fit_graph <- function(data, m=NULL, Gamma_emp=NULL){

  if(is.null(Gamma_emp)){
    Gamma_emp <- emp_vario(data)
  }

  Theta_emp <- Gamma2Theta(Gamma_emp)

  d <- nrow(Gamma_emp)

  if(is.null(m)){
    m <- d-1
  } else if(m < d-1){
    stop('m must be at least d-1!')
  }

  g_c <- igraph::make_full_graph(d)
  weights <- Theta_emp[igraph::as_edgelist(g_c)]
  g <- igraph::mst(g_c, weights)

  missing_edges <- igraph::as_edgelist(igraph::complementer(g))
  missing_weights <- Theta_emp[missing_edges]

  ord <- order(abs(missing_weights), decreasing = TRUE)

  new_edges <- missing_edges[ord[seq_len(m - (d-1))],]

  g <- igraph::add_edges(g, t(new_edges))

  Gamma_g <- complete_Gamma(Gamma_emp, g, N=100000, tol=1e-6, check_tol=100)

  return(list(
    graph = g,
    Gamma = Gamma_g
  ))
}


