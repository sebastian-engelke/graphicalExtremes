


#' @export 
ar_gauss <- function(Sigma, mu=NULL, graph=NULL, tol=1e-6){
  # Check input
  if(is.null(graph)){
    P <- solve(Sigma)
    A <- abs(P) > tol
    graph <- igraph::graph_from_adjacency_matrix(A, mode='undirected')
  }
  if(!igraph::is.chordal(graph)$chordal || !igraph::is.connected(graph)){
    stop('graph must be connected, decomposable.')
  }
  d <- nrow(Sigma)
  if(is.null(mu)){
    mu <- numeric(d)
  }
  
  # Make separators etc.
  cliques <- order_cliques(igraph::maximal.cliques(graph))
  tmp <- c()
  JList <- list()
  for(i in seq_along(cliques)){
    tmp <- union(tmp, cliques[[i]])
    JList[[i]] <- tmp
  }
  separators <- lapply(seq_len(length(cliques)-1)+1, function(i){
    intersect(cliques[[i]], JList[[i-1]])
  })
  EList <- lapply(seq_len(length(cliques)-1)+1, function(i){
    setdiff(cliques[[i]], JList[[i-1]])
  })
  
  # Compute vectors/matrices
  cli1 <- cliques[[1]]
  mu1 <- mu[cli1]
  Sigma1 <- Sigma[cli1, cli1]
  
  SigmaList <- list()
  MList <- list()
  muList <- list()
  for(i in seq_len(length(cliques) - 1)){
    E <- EList[[i]]
    D <- separators[[i]]
    SDDinv <- solve(Sigma[D, D])
    SigmaList[[i]] <- Sigma[E,E] - Sigma[E,D] %*% SDDinv %*% Sigma[D,E]
    muList[[i]] <- mu[E] - Sigma[E,D] %*% SDDinv %*% mu[D]
    MList[[i]] <- Sigma[E,D] %*% SDDinv
  }
  
  return(list(
    cliques = cliques,
    mu1 = mu1,
    Sigma1 = Sigma1,
    SigmaList = SigmaList,
    MList = MList,
    muList = muList
  ))
}


