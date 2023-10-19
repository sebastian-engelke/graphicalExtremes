


#' Non-decomposable completion of variogram matrices
#'
#' Given a non-decomposable `graph`, and (non-graphical) variogram matrix `Gamma`,
#' modifies `Gamma` in non-edge entries, such that the resulting matrix is a
#' variogram matrix with graphical structure described by `graph`.
#' Does so by splitting `graph` at complete separators into smaller subgraphs,
#' and calling `complete_Gamma_general` for each subgraph/submatrix,
#' using multiple cores if available.
#'
#' @param Gamma Numeric \dxd variogram matrix.
#' @param graph `igraph::graph()` object.
#' @param N Maximum number of iterations.
#' @param sub_tol Numeric scalar. Tolerance to be used when completing submatrices.
#' Should be smaller than `final_tol`.
#' @param check_tol Numeric/integer scalar. How often to check the tolerance when completing submatrices.
#' @param mc_cores_overwrite `NULL` or numeric/integer scalar. Maximal number of cores to use.
#' @param final_tol Numeric scalar. Check convergence of the final result with this tolerance.
#' Skipped if this value is < 0.
#'
#' @return A completed \dxd variogram matrix.
#' 
#' @family matrixCompletions
#' @export
complete_Gamma_general_split <- function(
  Gamma,
  graph,
  N = 10000,
  sub_tol = get_large_tol() * 1e-3,
  check_tol = 100,
  mc_cores_overwrite = NULL,
  final_tol = get_large_tol()
){
  # Check/find mc_cores
  mc_cores <- get_mc_cores(mc_cores_overwrite)
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
    mc.cores = mc_cores,
    complete_Gamma_general,
    subMatrices[needsCompletion],
    invSubGraphs[needsCompletion],
    MoreArgs = list(
      N = N,
      tol = sub_tol,
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
    Theta <- Gamma2Theta(Gamma, check = FALSE)
    err <- max(abs(getNonEdgeEntries(Theta, graph)))
    if(err > final_tol){
      warning('Matrix completion did not converge (err = ', err, ')')
    }
  }

  return(Gamma)
}


#' Non-decomposable completion of variogram matrices
#'
#' Given a non-decomposable `graph`, and (non-graphical) variogram matrix `Gamma`,
#' modifies `Gamma` in non-edge entries, such that the resulting matrix is a
#' variogram matrix with graphical structure described by `graph`.
#'
#' @param Gamma Numeric \dxd variogram matrix.
#' @param graph `igraph::graph()` object.
#' @param N Maximum number of iterations.
#' @param tol Numeric scalar. Tolerance to be used when completing submatrices.
#' @param check_tol Numeric/integer scalar. How often to check the tolerance when completing submatrices.
#'
#' @return A completed \dxd variogram matrix.
#' 
#' @family matrixCompletions
#' @export
complete_Gamma_general <- function(Gamma, graph, N=10000, tol=get_large_tol(), check_tol=100){
  
  if(is_complete_graph(graph)){
    return(Gamma)
  }

  detailedSepList <- make_sep_list(graph)

  nonEdgeIndices <- getNonEdgeIndices(graph, 'upper', doWhich = TRUE)

  GammaComp <- iterateIndList(Gamma, detailedSepList, nonEdgeIndices, N, tol, check_tol)

  return(GammaComp)
}

# Workhorse of `complete_Gamma_general()`
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
      Sigma <- Gamma2Sigma(Gamma, k = k, full = TRUE, check = FALSE)
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
      Gamma <- Sigma2Gamma(Sigma, check = FALSE)
    }

    # Continue iteration if no tol-check due
    if(check_tol <= 0 || n %% check_tol != 0){
      next
    }

    # Check if tolerance has been reached
    P <- Gamma2Theta(Gamma, check = FALSE)
    err <- max(abs(P[nonEdgeIndices]))
    if(err <= tol){
      break
    }
  }
  Gamma <- ensure_matrix_symmetry(Gamma)
  return(Gamma)
}






#' DEMO-VERSION: Completion of non-decomposable Gamma matrices
#'
#' Given a `graph` and variogram matrix `Gamma`, returns the full `Gamma`
#' matrix implied by the conditional independencies.
#' DEMO VERSION: Returns a lot of details and allows specifying the graph list
#' that is used. Is way slower than other functions.
#'
#' @param Gamma A complete variogram matrix (without any graphical structure).
#' @param graph An [`igraph::graph`] object.
#' @param N The maximal number of iterations of the algorithm.
#' @param tol The tolerance to use when checking for zero entries in `Theta`.
#' @param gList A list of graphs to be used instead of the output from [make_sep_list()].
#'
#' @return A nested list, containing the following details.
#' The "error term" is the maximal absolute value of `Theta` in a non-edge entry.
#' \item{graph, N, tol}{As in the input}
#' \item{gList}{As in the input or computed by [make_sep_list()].}
#' \item{Gamma0, Theta0, err0}{Initial `Gamma`, `Theta`, and error term.}
#' \item{iterations}{
#'  A nested list, containing the following infos for each performed iteration:
#'  \describe{
#'    \item{`n`}{Number of the iteration}
#'    \item{`t`}{Corresponding index in `gList`}
#'    \item{`g`}{The graph used}
#'    \item{`Gamma`, `Theta`, `err`}{The value of `Gamma`, `Theta`, and error term after the iteration}
#'  }
#' }
#' 
#' @family matrixCompletions
#' @export
complete_Gamma_general_demo <- function(Gamma, graph, N = 1000, tol=0, gList=NULL) {
  # Compute gList if not provided:
  if(is.null(gList)){
    sepDetails <- make_sep_list(graph, details=TRUE)
    gList <- lapply(sepDetails, function(dets) dets$graph)
  }
  m <- length(gList)

  # Initialize ret-list with initial Gamma, Theta, etc.:
  Theta <- Gamma2Theta(Gamma, check = FALSE)
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
    Theta <- Gamma2Theta(Gamma, check = FALSE)
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
