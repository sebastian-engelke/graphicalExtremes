


#' @export 
ar_gauss <- function(Sigma, mu=NULL, graph=NULL, tol=1e-6){
  # Check input
  if(is.null(graph)){
    P <- solve(Sigma)
    A <- abs(P) > tol
    diag(A) <- 0
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
  tmp <- computeVertexSets(cliques)
  separators <- tmp$separators
  JList <- tmp$JList
  EList <- tmp$EList

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
    M <- Sigma[E,D] %*% solve(Sigma[D, D])
    # or: M <- t(solve(Sigma[D,D], t(Sigma[E,D])))
    SigmaList[[i]] <- Sigma[E,E] - M %*% Sigma[D,E]
    muList[[i]] <- mu[E] - M %*% mu[D]
    MList[[i]] <- M
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

ar2Sigma_gauss <- function(cliques, mu1, Sigma1, SigmaList, MList, muList){
  d <- max(do.call(`c`, cliques))
  Sigma <- matrix(0, d, d)
  mu <- numeric(d)
  
  # Compute separators etc.
  tmp <- computeVertexSets(cliques)
  separators <- tmp$separators
  JList <- tmp$JList
  EList <- tmp$EList

  # Reconstruct Sigma, mu
  cli1 <- cliques[[1]]
  mu[cli1] <- mu1
  Sigma[cli1, cli1] <- Sigma1
  
  for(i in seq_len(length(cliques) - 1)){
    D <- separators[[i]]
    E <- EList[[i]]
    J <- JList[[i]]
    M <- MList[[i]]
    mu[E] <- muList[[i]] + M %*% mu[D]
    Sigma[E, E] <- M %*% Sigma[D, D] %*% t(M) + SigmaList[[i]]
    Sigma[J, E] <- Sigma[J, D] %*% t(M)
    Sigma[E, J] <- t(Sigma[J, E])
  }
  
  return(list(
    mu = mu,
    Sigma = Sigma
  ))
}



computeVertexSets <- function(cliques){
  # Make separators etc.
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
  return(list(
    cliques = cliques,
    separators = separators,
    JList = JList,
    EList = EList
  ))
}


