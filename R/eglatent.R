

#' Learning extremal graph structure with latent variables
#' 
#' Following the methodology from \insertCite{engelkeTaeb2024;textual}{graphicalExtremes},
#' fits an extremal graph structure with latent variables.
#' 
#' 
#' @param Gamma conditionally negative semidefinite matrix. This will be typically the empirical variogram matrix.
#' 
#' @param lam1_list Numeric vector of non-negative regularization parameters for eglatent.
#' Default is `lam1_list = c(0.1, 0.15, 0.19, 0.205)`.
#' 
#' @param lam2_list Numeric vector of non-negative regularization parameters for eglatent.
#' Default is `lam2_list = c(2)`.
#' 
#' @param refit Logical scalar, if TRUE then the model is refit on the estimated graph to obtain an estimate of the Gamma matrix on that graph.
#' Default is `refit = TRUE`.
#' 
#' @param verbose Logical scalar, indicating whether to print progress updates.
#' 
#' @return The function fits one model for each combination 
#' of values in `lam1_list` and `lam2_list`. All returned objects 
#' have one entry per model. List consisting of:
#' #' \item{`graph`}{
#'   A list of [`igraph::graph`] objects representing the
#'   fitted graphs.
#' }
#' \item{`rk`}{
#'   Numeric vector containing the estimated ranks of the latent
#'   variables.
#' }
#' \item{`G_est`}{
#'   A list of numeric estimated \dxd
#'   variogram matrices \eGamma corresponding to the fitted graphs. 
#' }
#' \item{`G_refit`}{
#'   A list of numeric estimated \dxd
#'   variogram matrices \eGamma refitted with fixed graphs 
#'   corresponding to the fitted graphs.
#' }
#' \item{`lambda_list`}{
#'   A list containing the values of `lam1_list` and `lam2_list` 
#'   used for the model fit.
#' }
#' 
#' @references \insertAllCited{}
#'
#' @family structureEstimation
#' @export
eglatent <- function(
  Gamma,
  lam1_list = c(0.1, 0.15, 0.19, 0.205),
  lam2_list = c(2),
  refit = TRUE,
  verbose = FALSE
){
  # Log function to print updates if specified
  logcat <- function(...){
    if(verbose){
      cat(...)
    }
  }

  d <- nrow(Gamma)
  U <- svd(diag(d) - 1 / d * matrix(1, d, 1) %*% t(matrix(1, d, 1)))$u[, 1:as.numeric(d - 1)]
  r <- 1
  graph_obs <- list()
  rk_vec <- numeric()
  G_est <- list()
  G_refit <- list()
  lambda_list <- list()

  logcat('Iterating over lambda values...\n')
  for (lambda1_iter in 1:length(lam1_list)) {
    logcat("\n\nLambda 1:", lam1_list[lambda1_iter], "(", lambda1_iter, "/", length(lam1_list), ")\n\n\n")
    for (lambda2_iter in 1:length(lam2_list)) {
      logcat("Lambda 2:", lam2_list[lambda2_iter], "(", lambda2_iter, "/", length(lam2_list), ")\n")

      lambda_1 <- lam1_list[lambda1_iter]
      lambda_2 <- lam2_list[lambda2_iter]

      lambda_list[[r]] <- c(lambda_1, lambda_2)

      # run the sparse+low-rank estimator using CVXR
      P <- CVXR::Variable(d, d, PSD = TRUE)
      L <- CVXR::Variable(d, d, PSD = TRUE)
      S <- CVXR::Variable(d, d, PSD = TRUE)
      R <- -CVXR::log_det(t(U) %*% P %*% U) - 1 / 2 * sum(CVXR::diag(P %*% Gamma)) + lambda_1 * (sum(sum(abs(S))) + lambda_2 * sum(CVXR::diag(L)))
      objective <- CVXR::Minimize(R)
      constraints <- list(P == S - L, U %*% t(U) %*% (P) %*% U %*% t(U) == P)
      prob <- CVXR::Problem(objective, constraints)
      result <- CVXR::solve(prob)
      G_est[[r]] <- graphicalExtremes:::ensure_matrix_symmetry(Theta2Gamma(result$getValue(P), check = FALSE))

      rk_vec[r] <- length(which(eigen(result$getValue(L))$values >= 10^(-3)))
      if (rk_vec[r] == 0) {
        subspace_est <- matrix(0, d, 1)
      } else {
        subspace_est <- base::svd(result$getValue(L))$u[, 1:rk_vec[r], drop = FALSE]
      }

      output_1 <- result$getValue(S)
      output_1[which(abs(output_1) <= 10^(-3))] <- 0
      output_1[which(abs(output_1) > 10^(-3))] <- 1
      output_1 <- output_1 - diag(diag(output_1))
      off_support_est <- which(output_1 + diag(d) == 0, arr.ind = TRUE)
      graph_obs[[r]] <- igraph::graph_from_adjacency_matrix(output_1, mode = "undirected")

      if (refit) {
        A <- matrix(0, 1, d^2)
        if (nrow(off_support_est) > 0) {
          A <- matrix(0, nrow(off_support_est), d^2)
          for (i in 1:nrow(off_support_est)) {
            S <- matrix(0, d, d)
            S[off_support_est[i, 1], off_support_est[i, 2]] <- 1
            A[i, ] <- c(S)
          }
        }
        rk <- rk_vec[r]
        if (rk != 0) {
          P <- CVXR::Variable(d, d, PSD = TRUE)
          M <- CVXR::Variable(rk, rk, PSD = TRUE)
          S <- CVXR::Variable(d, d, PSD = TRUE)
          R <- -CVXR::log_det(t(U) %*% P %*% U) - 1 / 2 * sum(CVXR::diag(P %*% Gamma)) 
          objective <- CVXR::Minimize(R)
          constraints <- list(P == S - subspace_est %*% M %*% t(subspace_est), U %*% t(U) %*% (P) %*% U %*% t(U) == P, A %*% reshape_expr(S, c(d^2, 1)) == 0)

          prob <- CVXR::Problem(objective, constraints)
          result <- CVXR::solve(prob)
        } else {
          P <- CVXR::Variable(d, d, PSD = TRUE)
          S <- CVXR::Variable(d, d, PSD = TRUE)
          R <- -CVXR::log_det(t(U) %*% P %*% U) - 1 / 2 * sum(CVXR::diag(P %*% Gamma)) 
          objective <- CVXR::Minimize(R)
          constraints <- list(P == S, U %*% t(U) %*% (P) %*% U %*% t(U) == P, A %*% CVXR::reshape_expr(S, c(d^2, 1)) == 0)

          prob <- CVXR::Problem(objective, constraints)
          result <- CVXR::solve(prob)
        }

        G_refit[[r]] <- ensure_matrix_symmetry(Theta2Gamma(result$getValue(P), check = FALSE))
      }
      r <- r + 1
    }
  }

  return(list(
    graph = graph,
    rk = rk_vec,
    G_est = G_est,
    G_refit = G_refit,
    lambda_list = lambda_list
  ))
}
