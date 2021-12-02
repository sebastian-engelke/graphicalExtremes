
##Functions+settings for emtp2
settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-10, eps_rel = 1e-10)

#' @export
start.positive.dual <- function(S){
  Z <- golazo::Zmatrix(S)
  return(Sigma2Gamma(Z))
}


#' Performs Gaussian likelihood optimization under Laplacian matrix constraints.
#'
#' This function implements a block descent algorithm to find the maximum of the 
#' Gaussian log-likelihood under the constraint that the concentration matrix is a Laplacian matrix.
#' @param G conditionally negative semidefinite matrix. This will be typically the empirical variogram matrix.
#' @param tol The convergence tolerance (default tol=1e-7). The algorithm terminates when the dual gap (guaranteed to be nonnegative) is less than tol.
#' @param initial_point if TRUE (default), the algorithm will construct an initial point before the iteration steps.
#' @param verbose if TRUE (default) the output will be printed.
#' @return `G_emtp2` the optimal value of the variogram matrix
#' @return `it` the number of iterations
#' @keywords block descent, concentration matrix, Laplacian matrix.
#' @export
emtp2 <- function(G,tol=1e-6, verbose = TRUE, initial_point = TRUE){
  d <- nrow(G)
  if (verbose==TRUE){
    cat("** The function maximizes the log-likelihood function under Laplacian matrix constraints.\n")
  }
  Gam <- G
  if (initial_point==TRUE){
    P <- diag(d)-matrix(1,d,d)/(d)
    Gam <- start.positive.dual(P%*%(-G/2)%*%P)
  }
  it <- 0
  if (verbose==TRUE){
    cat("\n The algorithm will stop when the duality gap is below: ",tol,"\b.\n\n")
    cat("Iteration | Gap\n")
  }
  gap <- Inf
  while (gap>tol){
    Gam0 <- Gam
    for (i in 1:d) {
      A <- solve((-Gam/2)[-i,-i])
      Dmat <- 2*(A%*%matrix(1,d-1,d-1)%*%A-sum(A)*A)
      dvec <- -2*A%*%rep(1,d-1)
      bvec <- (-G/2)[-i,i]
      y <- osqp::solve_osqp(P = Dmat, q = dvec, A = diag(d-1), l = bvec, u = rep(0,d-1), settings)
      Gam[-i,i] <- Gam[i,-i] <- -2*y$x
    }
    it <- it+1
    gap <- sum(abs(Gam-Gam0))
    if (verbose==TRUE){
      cat(it,"\t  | ",gap,"\n")
    }
  }
  return(list(G_emtp2=Gam,it=it))
}


