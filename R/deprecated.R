
#' Fitting extremal graphical lasso (MH: DEPRECATED?)
#'
#' Fits an extremal minimum spanning tree
#'
#' @param Gamma Numeric matrix \eqn{n\times d}{n x d}.
#' It represents a variogram matrix \eqn{\Gamma}.
#'
#' @param rholist Numeric vector of non-negative regularization parameters
#' for the lasso. For details see [glasso::glassopath].
#'
#' @param reg_method One of `"mb"` and `"glasso"`.
#' Default is `reg_method = "mb"`.
#'
#' @param eps Regularization parameter for the covariance matrix.
#' Default is `eps = 0.5`.
#'
#' @param complete_Gamma Whether you want to complete Gamma matrix.
#' Default is `complete_Gamma = FALSE`.
#'
#'
#' @return List made of:
#' \describe{
#'   \item{`graph`}{A list of [igraph::graph] objects representing the
#'   fitted graphs for each `rho` in `rholist`.}
#'   \item{`Gamma`}{A list of numeric \eqn{d\times d}{d x d} estimated
#'   variogram matrices \eqn{\Gamma} corresponding to the fitted graphs,
#'   for each `rho` in `rholist`.}
#'   \item{`rholist`}{The list of penalty coefficients.}
#' }
#'
#'
#' @export
eglasso <- function(Gamma, rholist= c(0.1, 0.15, 0.19, 0.205),
                    reg_method =  c("mb", "glasso"),
                    eps=0.5, complete_Gamma = FALSE){

  # Check args
  reg_method <- match.arg(reg_method)
  if (any(rholist < 0)) {
    stop("The regularization parameters in `rholist` must be non-negative.",
         call. = FALSE)
  }

  # Set main variables
  r <- length(rholist)
  d <- ncol(Gamma)
  null.vote <- array(
    0,
    dim=c(d, d, length(rholist))) # votes for EXCLUDING the edge

  for(k in 1:d){
    Sk <- Gamma2Sigma(Gamma=Gamma, k=k) #stats::cov2cor(Gamma2Sigma(Gamma=Gamma, k=k)) #Not normalizing seems slighlty more stable
    ###### Same regularization, but does not require Sk to be invertible
    tmp <- solve(diag(ncol(Sk)) + eps*Sk)
    Ck <- tmp %*% Sk

    ###### Using "glasso" package
    approx <- (reg_method == "mb")
    if(reg_method != "glasso" && reg_method != "mb") {
      warning(paste(
        "Method",
        reg_method,
        "not implemended in glasso. Regular glasso was used instead."),
        call. = FALSE)
    }

    invisible(utils::capture.output(
      gl.tmp <- glasso::glassopath(Ck, rholist = rholist, approx=approx)))
    null.vote[-k,-k, ] <-  null.vote[-k,-k, , drop = FALSE] +
      (abs(gl.tmp$wi)<=1e-10)  ## change back to == 0 .. <=1e-4

  }
  adj.est <- (null.vote/(ncol(null.vote)-2)) < .49


  graphs <- list()
  Gammas <- list()
  rhos <- list()

  for(j in 1:r) {
    rho <- rholist[j]
    est_graph <- igraph::graph_from_adjacency_matrix(adj.est[,,j],
                                                     mode="undirected",
                                                     diag=FALSE)

    if (complete_Gamma == FALSE) {
      Gamma_curr <- NA
    } else {
      Gamma_curr <- tryCatch({
        complete_Gamma(graph = est_graph, Gamma = Gamma)

      },
      error = function(e){
        if (e$message == "The given graph is not connected."){
          message(paste0("The estimated graph for rho = ", round(rho, 3),
                         " is not connected, ",
                         "so it is not possible to complete Gamma.\n"))

          NA

        } else {
          stop(e)
        }
      })

      if (all(!is.na(Gamma_curr))) {
        completed_graph <- Gamma2graph(Gamma_curr, to_plot = FALSE)

        if (!(graphs_equal(completed_graph, est_graph))) {
          message(paste0("The completed Gamma for rho = ", round(rho, 3),
                         " does not match the estimated graph.\n"))
        }
      }
    }

    graphs[[j]] <- est_graph
    Gammas[[j]] <- Gamma_curr
    rhos[[j]] <- rho

  }

  return(list(graph = graphs, Gamma = Gammas, rholist = rhos))
}

