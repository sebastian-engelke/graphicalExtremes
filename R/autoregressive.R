


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


ar_HR <- function(Gamma, graph=NULL, tol=1e-6, GammaSub=TRUE){
  # Check input
  if(is.null(graph)){
    P <- Gamma2Theta(Gamma)
    A <- abs(P) > tol
    diag(A) <- 0
    graph <- igraph::graph_from_adjacency_matrix(A, mode='undirected')
  }
  if(!igraph::is.chordal(graph)$chordal || !igraph::is.connected(graph)){
    stop('graph must be connected, decomposable.')
  }
  d <- nrow(Gamma)

  # Make separators etc.
  cliques <- order_cliques(igraph::maximal.cliques(graph))
  tmp <- computeVertexSets(cliques, makeSep1=TRUE)
  separators <- tmp$separators
  JList <- tmp$JList
  EList <- tmp$EList
  
  kList <- lapply(separators, function(sep) sep[1])

  # Compute vectors, matrices
  SigmaList <- list()
  MList <- list()
  for(i in seq_along(cliques)){
    # k <- kList[[i]]
    MList[[i]] <- list()
    SigmaList[[i]] <- list()
    for(j in seq_along(separators[[i]])){
      k <- separators[[i]][[j]]
      D <- separators[[i]]
      D_k <- setdiff(separators[[i]], k) # "D without k"
      E_k <- setdiff(EList[[i]], k) # "E without k"

      J <- JList[[i]]

      if(GammaSub){
        k2 <- which(J == k)
        Sk <- matrix(0, d, d)
        Sk[J, J] <- Gamma2Sigma(Gamma[J, J], k=k2, full=TRUE)
      } else{
        Sk <- Gamma2Sigma(Gamma, k=k, full=TRUE)
      }
      if(length(D_k) > 0){
        B <- Sk[E_k,D_k] %*% solve(Sk[D_k,D_k])
        SigmaList[[i]][[j]] <- Sk[E_k,E_k] - B %*% Sk[D_k,E_k]
        M <- matrix(0, nrow(B), ncol(B) + 1)
        M[, D != k] <- B
        M[, D == k] <- 1 - rowSums(B)
        MList[[i]][[j]] <- M
      } else{
        SigmaList[[i]][[j]] <- Sk[E_k, E_k]
        MList[[i]][[j]] <- matrix(1, length(E_k), 1)
      }
    }
  }
  
  return(list(
    cliques = cliques,
    SigmaList = SigmaList,
    MList = MList
  ))
}


computeVertexSets <- function(cliques, makeSep1=FALSE){
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
  if(makeSep1){
    k1 <- cliques[[1]][1]
    separators <- c(list(k1), separators)
    EList <- c(list(setdiff(cliques[[i]], k1)), EList)
  }
  return(list(
    cliques = cliques,
    separators = separators,
    JList = JList,
    EList = EList
  ))
}


