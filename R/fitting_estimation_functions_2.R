

#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on decomposable graphs
#' 
#' Similar to [fmpareto_graph_HR()]. Differences are:
#' * Works with any decomposable graph, not just block graphs
#' * Does not support `edges_to_add`
fmpareto_graph_HR_2 <- function(data, graph, p = NULL, cens = FALSE) {
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
  
  return(list(graph = set_graph_parameters(graph), Gamma = Ghat))
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


