
#' Completion of Gamma matrices
#' 
#' Given a `graph` and a (partial) variogram matrix `Gamma`, returns a full
#' variogram matrix that agrees with `Gamma` in entries corresponding to edges
#' of `graph` and whose corresponding precision matrix, obtained by
#' [Gamma2Theta()], has zeros in entries corresponding to non-edges of `graph`.
#' For results on the existence and uniqueness of this completion, see
#' \insertCite{hen2022;textual}{graphicalExtremes}.
#' 
#' @param Gamma Numeric \dxd variogram matrix.
#' @param graph `NULL` or [`igraph::graph`] object. If `NULL`, the graph
#' is implied by non-edge entries in `Gamma` being `NA`. Must be connected, undirected.
#' @param ... Further arguments passed to [complete_Gamma_general_split()] if `graph`
#' is not decomposable
#'
#' @return Completed \dxd `Gamma` matrix.
#'
#' @details
#' If `graph` is decomposable, `Gamma` only needs to be specified on
#' the edges of the graph, other entries are ignored.
#' If `graph` is not decomposable, the graphical completion algorithm requires
#' a fully specified (but non-graphical) variogram matrix `Gamma` to begin with.
#' If necessary, this initial completion is computed using [edmcr::npf()].
#'
#' @examples
#' ## Block graph:
#' Gamma <- rbind(
#'   c(0, .5, NA, NA),
#'   c(.5, 0, 1, 1.5),
#'   c(NA, 1, 0, .8),
#'   c(NA, 1.5, .8, 0)
#' )
#'
#' complete_Gamma(Gamma)
#'
#' ## Alternative representation of the same completion problem:
#' my_graph <- igraph::graph_from_adjacency_matrix(rbind(
#'   c(0, 1, 0, 0),
#'   c(1, 0, 1, 1),
#'   c(0, 1, 0, 1),
#'   c(0, 1, 1, 0)
#' ), mode = "undirected")
#' Gamma_vec <- c(.5, 1, 1.5, .8)
#' complete_Gamma(Gamma_vec, my_graph)
#'
#' ## Decomposable graph:
#' G <- rbind(
#' c(0, 5, 7, 6, NA),
#' c(5, 0, 14, 15, NA),
#' c(7, 14, 0, 5, 5),
#' c(6, 15, 5, 0, 6),
#' c(NA, NA, 5, 6, 0)
#' )
#'
#' complete_Gamma(G)
#'
#' ## Non-decomposable graph:
#' G <- rbind(
#' c(0, 5, 7, 6, 6),
#' c(5, 0, 14, 15, 13),
#' c(7, 14, 0, 5, 5),
#' c(6, 15, 5, 0, 6),
#' c(6, 13, 5, 6, 0)
#' )
#' g <- igraph::make_ring(5)
#'
#' complete_Gamma(G, g)
#'
#'
#' @references \insertAllCited{}
#'
#' @seealso [Gamma2Theta()]
#' @family Matrix completions
#'
#' @export
complete_Gamma <- function(
  Gamma,
  graph = NULL,
  ...
){
  tmp <- check_Gamma_and_graph(Gamma, graph)
  Gamma <- tmp$Gamma
  graph <- tmp$graph

  # Return completion if graph is decomposable (=chordal)
  if(igraph::is_chordal(graph)$chordal){
    return(complete_Gamma_decomposable(Gamma, graph))
  }

  # Compute initial non-graphical completion if necessary:
  if(any(is.na(Gamma))){
    A <- 1*!is.na(Gamma)
    tmp <- edmcr::npf(Gamma, A, d = NROW(Gamma)-1)
    Gamma <- tmp$D
    if(!is_sym_cnd(Gamma)){
      stop('Did not find an initial non-graphical completion (using edmcr::npf)!')
    }
  }

  # Compute non-decomposable completion
  return(complete_Gamma_general_split(Gamma, graph, ...))
}



#' Completion of decomposable Gamma matrices
#'
#' Given a decomposable `graph` and incomplete variogram matrix `Gamma`,
#' returns the full `Gamma` matrix implied by the conditional independencies.
#'
#' @param Gamma A variogram matrix that is specified on the edges of `graph`
#' and the diagonals. All other entries are ignored (if `graph` is specified),
#' or should be `NA` to indicate non-edges in `graph`.
#' @param graph `NULL` or a decomposable \[`igraph::graph`\] object. If `NULL`, the
#' structure of `NA` entries in `Gamma` is used instead.
#'
#' @return A complete variogram matrix that agrees with `Gamma` on the entries
#' corresponding to edges in `graph` and the diagonals.
#' The corresponding \eTheta matrix produced by [Gamma2Theta()] has zeros
#' in the remaining entries.
#'
#' @family Matrix completions
#'
#' @export
complete_Gamma_decomposable <- function(Gamma, graph = NULL) {
  # Compute graph if not specified
  if(is.null(graph)){
    graph <- partialMatrixToGraph(Gamma)
  }

  # compute cliques and order by the running intersection property
  cliques <- igraph::max_cliques(graph)
  cliques <- order_cliques(cliques)

  # Start with first clique
  oldVertices <- cliques[[1]]

  # loop through remaining cliques. Skipped if only one clique in graph.
  for (p in seq_along(cliques)[-1]) {
    newVertices <- cliques[[p]]
    vC <- intersect(oldVertices, newVertices)
    vA <- setdiff(oldVertices, vC)
    vB <- setdiff(newVertices, vC)

    vACB <- c(vA, vC, vB)

    GammaACB <- Gamma[vACB, vACB]
    GammaACB <- complete_Gamma_one_step(GammaACB, length(vA), length(vC), length(vB))
    Gamma[vACB, vACB] <- GammaACB

    oldVertices <- union(oldVertices, newVertices)
  }
  Gamma <- ensure_symmetry(Gamma)
  return(Gamma)
}


#' Completion of two-clique decomposable Gamma matrices
#'
#' Given a decomposable `graph` consisting of two cliques and incomplete
#' variogram matrix `Gamma`, returns the full `Gamma` matrix implied by the
#' conditional independencies. The rows/columns of `Gamma` must be ordered
#' such that the clique of size `nA` (excluding separator) comes first, then
#' the separator of size `nC`, and then the remaining `nB` vertices.
#'
#' @keywords internal
complete_Gamma_one_step <- function(Gamma, nA, nC, nB) {
  # Check arguments
  n <- nA + nB + nC
  if (nrow(Gamma) != n || ncol(Gamma) != n) {
    stop("Make sure that nrow(Gamma) == ncol(Gamma) == nA+nB+nC")
  }

  # Prepare index sets
  vA <- seq_len(nA) # first nA entries
  vC <- seq_len(nC) + nA # next nC entries
  vB <- seq_len(nB) + nA + nC # remaining nB entries
  k0 <- vC[1] # condition on first entry in vC.
  vC_Sigma <- vC[-1]

  if(nC == 1){
    ## Separator of size 1 -> use additive property of block matrix completion
    GammaAB <- outer(Gamma[vA, k0], Gamma[k0, vB], `+`)
    Gamma[vA, vB] <- GammaAB
    Gamma[vB, vA] <- t(GammaAB)
    return(Gamma)
  }

  ## Larger separator -> convert to Sigma, complete, convert back
  Sigma <- Gamma2Sigma(Gamma, k = k0, full = TRUE)

  # Invert separator submatrix
  R <- chol(Sigma[vC_Sigma, vC_Sigma, drop=FALSE])
  SigmaCCinv <- chol2inv(R)
  # Compute completion
  SigmaAB <- Sigma[vA, vC_Sigma] %*% SigmaCCinv %*% Sigma[vC_Sigma, vB]

  # Fill in completed values
  Sigma[vA, vB] <- SigmaAB
  Sigma[vB, vA] <- t(SigmaAB)

  # Convert back
  Gamma <- Sigma2Gamma(Sigma)
  return(Gamma)
}



