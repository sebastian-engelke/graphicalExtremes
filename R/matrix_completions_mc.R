

complete_Gamma_general_mc <- function(Gamma, graph, N = 1000, tol=0, check_tol=100){
  # wip
  graph <- setPids(graph)
  invSubGraphs <- split_graph(graph)
  subMatrices <- lapply(invSubGraphs, getSubMatrixForSubgraph, fullMatrix = Gamma)
  
  completedSubMatrices <- mapply(
    complete_Gamma_general,
    subMatrices,
    invSubGraphs,
    MoreArgs = list(
      N = N,
      tol = tol,
      check_tol = check_tol
    ),
    SIMPLIFY = FALSE
  )
  
  GammaDecomp <- NA * Gamma
  for(i in seq_along(invSubGraphs)){
    inds <- getIdsForSubgraph(invSubGraphs[[i]], graph)
    GammaDecomp[inds, inds] <- completedSubMatrices[[i]]
  }
  
  A <- 1*(!is.na(GammaDecomp))
  diag(A) <- 0
  gDecomp <- igraph::graph_from_adjacency_matrix(A, mode='undirected')
  Gamma <- complete_Gamma_decomposable(GammaDecomp, gDecomp)
  
  return(Gamma)
}


iterateIndList <- function(Gamma, indList, nonEdgeIndices, N, tol, check_tol){
  for (n in seq_len(N)) {
    # Read indices of the two cliques (A, B) and separator (C) to be used
    t <- (n - 1) %% m + 1
    inds <- indList[[t]]
    vA <- inds$vA
    vB <- inds$vB
    vC <- inds$vC
    vC_Sigma <- inds$vC_Sigma
    k0 <- inds$k0
    
    # Separators of size 1 don't happen if graph is decomposed first!
    
    ## Compute 1 iteration:
    # Compute Sigma
    Sigma <- Gamma2Sigma(Gamma, k = k0, full = TRUE)
    # Invert Sigma_CC
    R <- chol(Sigma[vC_Sigma, vC_Sigma, drop=FALSE])
    SigmaCCinv <- chol2inv(R)
    # Compute completion:
    SigmaAB <- Sigma[vA, vC_Sigma, drop=FALSE] %*% SigmaCCinv %*% Sigma[vC_Sigma, vB, drop=FALSE]
    # Store completion and convert back:
    Sigma[vA, vB] <- SigmaAB
    Sigma[vB, vA] <- t(SigmaAB)
    Gamma <- Sigma2Gamma(Sigma)
    
    # Continue iteration if no tol-check due
    if(check_tol <= 0 || n %% check_tol != 0){
      next
    }
    
    # Check if tolerance has been reached
    P <- Gamma2Theta(Gamma)
    err <- max(abs(P[nonEdgeIndices]))
    if(err <= tol){
      break
    }
  }
  return(Gamma)
}

