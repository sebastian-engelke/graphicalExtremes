

complete_Gamma_general_mc <- function(
  Gamma,
  graph,
  N = 10000,
  tol = 1e-12,
  check_tol = 100,
  mc.cores = 1,
  final_tol = 1e-9
){
  # wip
  graph <- setPids(graph)
  invSubGraphs <- split_graph(graph)
  subMatrices <- lapply(invSubGraphs, getSubMatrixForSubgraph, fullMatrix = Gamma)

  needsCompletion <- sapply(invSubGraphs, function(g){
    d <- igraph::vcount(g)
    maxEcount <- d * (d-1) / 2
    return(igraph::ecount(g) < maxEcount)
  })

  completedSubMatrices <- parallel::mcmapply(
    mc.cores = mc.cores,
    complete_Gamma_general_sc,
    subMatrices[needsCompletion],
    invSubGraphs[needsCompletion],
    MoreArgs = list(
      N = N,
      tol = tol,
      check_tol = check_tol
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE
  )

  subMatrices[needsCompletion] <- completedSubMatrices

  GammaDecomp <- NA * Gamma
  for(i in seq_along(invSubGraphs)){
    inds <- getIdsForSubgraph(invSubGraphs[[i]], graph)
    GammaDecomp[inds, inds] <- subMatrices[[i]]
  }

  gDecomp <- partialMatrixToGraph(GammaDecomp)
  Gamma <- complete_Gamma_decomposable(GammaDecomp, gDecomp)

  if(final_tol >= 0){
    Theta <- Gamma2Theta(Gamma)
    err <- max(abs(getNonEdgeEntries(Theta, graph)))
    if(err > tol){
      warning('Matrix completion did not converge (err = ', err, ')')
    }
  }

  return(Gamma)
}

complete_Gamma_general_sc <- function(Gamma, graph, N=10000, tol=1e-12, check_tol=100){
  
  if(is_complete_graph(graph)){
    return(Gamma)
  }

  detailedSepList <- make_sep_list(graph)

  nonEdgeIndices <- getNonEdgeIndices(graph, 'upper', doWhich = TRUE)

  GammaComp <- iterateIndList(Gamma, detailedSepList, nonEdgeIndices, N, tol, check_tol)

  return(GammaComp)
}

iterateIndList <- function(Gamma, sepList, nonEdgeIndices, N, tol, check_tol){
  m <- length(sepList)
  n <- 0
  while(n < N){
    # Read of the cliques (`parts`), separator (`sep`), etc. to be used
    n <- n+1
    t <- (n - 1) %% m + 1
    inds <- sepList[[t]]
    parts <- inds$parts
    sep <- inds$sep
    sep_Sigma <- inds$sepWithoutK
    k <- inds$k
    partPairs <- inds$partPairs

    ## Compute 1 iteration. Check case #sep=1 vs #sep>1:
    if(length(sep) == 1){
      # Separator is of size 1 -> we can simply add variogram entries
      for(partPair in partPairs){
        vA <- partPair[[1]]
        vB <- partPair[[2]]
        GammaAB <- outer(Gamma[vA, k], Gamma[k, vB], `+`)
        Gamma[vA, vB] <- GammaAB
        Gamma[vB, vA] <- t(GammaAB)
      }
    } else{
      # Separator is of size >1 -> we need to work with Sigma^{(k)}
      # Compute Sigma
      Sigma <- Gamma2Sigma(Gamma, k = k, full = TRUE)
      # Invert Sigma_CC
      R <- chol(Sigma[sep_Sigma, sep_Sigma, drop=FALSE])
      SigmaCCinv <- chol2inv(R)
      # Compute completion for each pair of separated parts
      for(partPair in partPairs){
        vA <- partPair[[1]]
        vB <- partPair[[2]]
        SigmaAB <- Sigma[vA, sep_Sigma, drop=FALSE] %*% SigmaCCinv %*% Sigma[sep_Sigma, vB, drop=FALSE]
        # Store completion
        Sigma[vA, vB] <- SigmaAB
        Sigma[vB, vA] <- t(SigmaAB)
      }
      # Convert back to Gamma
      Gamma <- Sigma2Gamma(Sigma)
    }

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
  Gamma <- ensure_symmetry(Gamma) # todo: set tolerance = Inf?
  return(Gamma)
}






#' Completion of non-decomposable Gamma matrices (demo-version)
#'
#' Given a `graph` and variogram matrix `Gamma`, returns the full `Gamma`
#' matrix implied by the conditional independencies.
#' This function uses a convergent iterative algorithm.
#' DEMO VERSION: Returns a lot of details and allows specifying the graph list
#' that is used. Is way slower than `complete_Gamma_general`.
#'
#' @param Gamma A complete variogram matrix (without any graphical structure)
#' @param graph An [igraph::graph] object
#' @param N The maximal number of iterations of the algorithm
#' @param tol The tolerance to use when checking for zero entries in `Theta`
#'
#' @return A nested list, containing the following details: TODO!!!
#' 
#' The corresponding \eTheta matrix produced by [Gamma2Theta()] has values
#' close to zero in the remaining entries (how close depends on the input
#' and the number of iterations).
#'
#' @family Matrix completions
complete_Gamma_general_demo <- function(Gamma, graph = NULL, N = 1000, tol=0, gList=NULL) {
  # Compute gList if not provided:
  if(is.null(gList)){
    sepDetails <- make_sep_list(graph, details=TRUE)
    gList <- lapply(sepDetails, function(dets) dets$graph)
  }
  m <- length(gList)

  # Initialize ret-list with initial Gamma, Theta, etc.:
  Theta <- Gamma2Theta(Gamma)
  iterations <- list()
  ret <- list(
    graph = graph,
    Gamma0 = Gamma,
    Theta0 = Theta,
    err0 = max(abs(getNonEdgeEntries(Theta, graph))),
    gList = gList,
    tol = tol,
    N = N
  )

  # Iterate over gList:
  n <- 0
  while(n < N) {
    n <- n + 1
    t <- (n - 1) %% m + 1
    g <- gList[[t]]
    Gamma <- complete_Gamma_decomposable(Gamma, g)

    GammaList <- c(GammaList, list(Gamma))

    # Compute Theta
    Theta <- Gamma2Theta(Gamma)
    err <- max(abs(getNonEdgeEntries(Theta, graph)))

    # Store results
    iterations <- c(iterations, list(list(
      n = n,
      t = t,
      g = g,
      Gamma = Gamma,
      Theta = Theta,
      err = err
    )))

    # Check if tolerance has been reached
    if(err <= tol){
      break
    }
  }

  ret$iterations <- iterations

  return(ret)
}