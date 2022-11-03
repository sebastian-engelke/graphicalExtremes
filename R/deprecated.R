
#' DEPRECATED: Completion of non-decomposable Gamma matrices
#'
#' Given a `graph` and variogram matrix `Gamma`, returns the full `Gamma`
#' matrix implied by the conditional independencies.
#' This function uses a convergent iterative algorithm.
#'
#' @param Gamma A complete variogram matrix (without any graphical structure)
#' @param graph An [igraph::graph] object
#' @param N The maximal number of iterations of the algorithm
#' @param tol The tolerance to use when checking for zero entries in `Theta`
#' @param check_tol After how many iterations to check the tolerance in `Theta`
#'
#' @return A matrix that agrees with `Gamma` on the entries corresponding to
#' edges in `graph` and the diagonals.
#' The corresponding \eqn{\Theta} matrix produced by [Gamma2Theta] has values
#' close to zero in the remaining entries (how close depends on the input
#' and the number of iterations).
#'
#' @family Matrix completions
#' @export
DEPRECATED_complete_Gamma_general <- function(Gamma, graph, N = 1000, tol=0, check_tol=100) {

  tmp <- make_graph_list(graph)
  partitionList <- tmp$partitions
  gList <- tmp$graphs

  indList <- lapply(partitionList, function(AB) list(
    vC = intersect(AB$A, AB$B),
    vA = setdiff(AB$A, AB$B),
    vB = setdiff(AB$B, AB$A)
  ))
  m <- length(indList)

  if(length(gList) == 0){
    N <- 0
  }

  for (n in seq_len(N)) {
    t <- (n - 1) %% m + 1
    vABC <- indList[[t]]
    vA <- vABC$vA
    vB <- vABC$vB
    vC <- vABC$vC

    if(length(vC) == 1){
      k0 <- vC[1]
      GammaAB <- outer(Gamma[vA, k0], Gamma[k0, vB], '+')
      Gamma[vA, vB] <- GammaAB
      Gamma[vB, vA] <- t(GammaAB)
    } else{
      k0 <- vC[1]
      vC_Sigma <- vC[-1]
      Sigma <- Gamma2Sigma(Gamma, k = k0, full = TRUE)
      R <- chol(Sigma[vC_Sigma, vC_Sigma, drop=FALSE])
      SigmaCCinv <- chol2inv(R)
      SigmaAB <- Sigma[vA, vC_Sigma, drop=FALSE] %*% SigmaCCinv %*% Sigma[vC_Sigma, vB, drop=FALSE]
      Sigma[vA, vB] <- SigmaAB
      Sigma[vB, vA] <- t(SigmaAB)
      Gamma <- Sigma2Gamma(Sigma)
    }

    # Check if tolerance has been reached
    if(check_tol > 0 && n %% check_tol == 0){
      P <- Gamma2Theta(Gamma)
      A <- igraph::as_adjacency_matrix(graph, sparse=FALSE)
      diag(A) <- 1
      err <- max(abs(P[A == 0]))
      if(err <= tol){
        break
      }
    }
  }

  return(Gamma)
}


