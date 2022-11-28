#' Learning extremal graph structure
#'
#' Fits an extremal graph structure using the neighborhood selection approach
#' (see \insertCite{meins2006;textual}{graphicalExtremes}) or graphical lasso
#' (see \insertCite{friedman2008;textual}{graphicalExtremes}).
#'
#' @param data Numeric \nxd matrix, where `n` is the
#' number of observations and `d` is the dimension.
#'
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#'
#' @param rholist Numeric vector of non-negative regularization parameters
#' for the lasso.
#' Default is `rholist = c(0.1, 0.15, 0.19, 0.205)`.
#' For details see [glasso::glassopath].
#'
#' @param reg_method One of `"ns", "glasso"`, for neighborhood selection and
#' graphical lasso, respectively.
#' Default is `reg_method = "ns"`.
#' For details see \insertCite{meins2006;textual}{graphicalExtremes},
#' \insertCite{friedman2008;textual}{graphicalExtremes}.
#'
#' @param complete_Gamma Whether you want to try fto complete Gamma matrix.
#' Default is `complete_Gamma = FALSE`.
#'
#' @return List made of:
#' \item{`graph`}{
#'   A list of [igraph::graph] objects representing the
#'   fitted graphs for each `rho` in `rholist`.
#' }
#' \item{`Gamma`}{
#'   A list of numeric estimated \dxd
#'   variogram matrices \eGamma corresponding to the fitted graphs,
#'   for each `rho` in `rholist`. If `complete_Gamma = FALSE` or the
#'   underlying graph is not connected, it returns `NULL`.
#' }
#' \item{`rholist`}{
#'   The list of penalty coefficients.
#' }
#' \item{`graph_ic`}{
#'   A list of [igraph::graph] objects
#'   representing the optimal graph
#'   according to the `aic`, `bic`, and `mbic` information criteria.
#'   If `reg_method = "glasso"`, it returns a list of `NA`.
#' }
#' \item{`Gamma_ic`}{
#'   A list of numeric \dxd estimated
#'   variogram matrices \eGamma corresponding
#'   to the `aic`, `bic`, and `mbic` information criteria.
#'   If `reg_method = "glasso"`, `complete_Gamma = FALSE`, or the underlying
#'   graph is not connected, it returns a list of `NA`.
#' }
#'
#' @export
#'
eglearn <- function(
  data,
  p = NULL,
  rholist = c(0.1, 0.15, 0.19, 0.205),
  reg_method = c("ns", "glasso"),
  complete_Gamma = FALSE
) {
  # Check arguments
  reg_method <- match.arg(reg_method)
  if (any(rholist < 0)) {
    stop("The regularization parameters in `rholist` must be non-negative.")
  }

  # Normalize data
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  # Initialize variables
  Gamma <- emp_vario(data)
  sel_methods <- c("aic", "bic", "mbic")

  graphs <- list()
  Gammas <- list()
  rhos <- list()
  nRho <- length(rholist)
  d <- ncol(Gamma)

  null.vote <- array(0, dim = c(d, d, length(rholist)))
  null.vote.ic <- array(0, dim = c(d, d, length(sel_methods)))

  # Loop through variables
  for (k in 1:d) {
    if (reg_method == "glasso") {
      Sk <- Gamma2Sigma(Gamma = Gamma, k = k)
      gl.fit <- lapply(1:length(rholist), FUN = function(i) {
        glassoFast::glassoFast(S = Sk, rho = rholist[i], thr = 1e-8, maxIt = 100000)$wi
      })
      gl.tmp <- array(unlist(gl.fit), dim = c(d - 1, d - 1, nRho))
      null.vote[-k, -k, ] <- null.vote[-k, -k, , drop = FALSE] +
        (abs(gl.tmp) <= 1e-5) ## make sure is consistent with default threshold value
    } else if (reg_method == "ns") {
      idx_k <- which(data[, k] > 1)
      X <- log(data[idx_k, -k] / data[idx_k, k])
      gl.tmp <- glasso_mb(data = X, lambda = rholist)
      null.vote[-k, -k, ] <- null.vote[-k, -k, , drop = FALSE] + (!gl.tmp$adj.est)
      null.vote.ic[-k, -k, ] <- null.vote.ic[-k, -k, , drop = FALSE] + (!gl.tmp$adj.ic.est)
    }
  }

  adj.est <- (null.vote / (ncol(null.vote) - 2)) < 0.5
  # only makes sense for MB at the moment
  adj.ic.est <- (null.vote.ic / (ncol(null.vote.ic) - 2)) < 0.5


  # complete Gamma for rholist
  for (j in 1:nRho) {
    rho <- rholist[j]
    est_graph <- igraph::graph_from_adjacency_matrix(
      adj.est[, ,j], mode = "undirected", diag = FALSE)

    if (complete_Gamma == FALSE) {

      Gamma_curr <- NA
    } else {
      Gamma_curr <- try_complete_Gamma(est_graph, Gamma, "rho", round(rho, 3))
    }

    graphs[[j]] <- est_graph
    Gammas[[j]] <- Gamma_curr
    rhos[[j]] <- rho
  }

  # complete Gamma for ns
  graphs_ic <-  list(aic = NA, bic = NA, mbic = NA)
  Gammas_ic <-  list(aic = NA, bic = NA, mbic = NA)

  if (reg_method == "ns") {
    for (l in seq_along(sel_methods)){
      est_graph <-  igraph::graph_from_adjacency_matrix(
        adj.ic.est[, ,l], mode = "undirected", diag = FALSE
      )

      if (complete_Gamma == FALSE) {
        Gamma_curr <- NULL
      } else {
        Gamma_curr <- try_complete_Gamma(
          est_graph, Gamma,
          key = "information criterion",
          val = sel_methods[l]
        )
      }

      graphs_ic[[l]] <- est_graph
      Gammas_ic[[l]] <- Gamma_curr
    }
  }

  return(list(
    graph = graphs,
    Gamma = Gammas,
    rholist = rhos,
    graph_ic = graphs_ic,
    Gamma_ic = Gammas_ic
  ))
}


try_complete_Gamma <- function(graph, Gamma, key, val){
  ## igraph numeric_matrix character double -> numeric_matrix | NA
  ## tries to call `complete_Gamma`, if it fails returns NULL
  if(!igraph::is.connected(graph)){
    message(paste0(
      'The estimated graph for ',
      key, ' = ', val,
      '} is not connected,  so it is not possible to complete `Gamma`.\n'
    ))
    return(NULL)
  }

  Gamma_comp <- complete_Gamma(graph = graph, Gamma = Gamma)
  graph_comp <- Gamma2graph(Gamma_comp)

  # Check if completed Gamma matches with given graph
  if (!graphs_equal(graph_comp, graph)) {
    message(paste0(
      'The completed Gamma for ',
      key, ' = ', val,
      ' does not match the estimated graph.\n'
    ))
  }

  # Return completed Gamma
  return(Gamma_comp)
}




#' Fitting extremal minimum spanning tree
#'
#' Fits an extremal minimum spanning tree, where the edge weights are:
#' - negative maximized log-likelihoods of the bivariate Huesler--Reiss distributions, 
#'   if `method = "ML"`. See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#' - empirical extremal variogram, if `method = "vario"`. See \insertCite{eng2020;textual}{graphicalExtremes} for details.
#' - empirical extremal correlation, if `method = "chi"`. See \insertCite{eng2020;textual}{graphicalExtremes} for details.
#'
#' @param data Numeric \nxd matrix, where `n` is the
#' number of observations and `d` is the dimension.
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#' @param method One of `"vario", "ML", "chi"`.
#' Default is `method = "vario"`.
#' @param cens Logical. This argument is considered only if `method = "ML"`.
#' If `TRUE`, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#'
#' @return List consisting of:
#' \item{`graph`}{An [`igraph::graph`] object. The fitted minimum spanning tree.}
#' \item{`Gamma`}{
#'   Numeric \dxd estimated variogram matrix \eGamma
#'   corresponding to the fitted minimum spanning tree.
#' }
#'
#' @examples
#' ## Fitting a 4-dimensional HR minimum spanning tree
#' my_graph <- igraph::graph_from_adjacency_matrix(
#'   rbind(
#'     c(0, 1, 0, 0),
#'     c(1, 0, 1, 1),
#'     c(0, 1, 0, 0),
#'     c(0, 1, 0, 0)
#'   ),
#'   mode = "undirected"
#' )
#' n <- 100
#' Gamma_vec <- c(.5, 1.4, .8)
#' complete_Gamma(Gamma = Gamma_vec, graph = my_graph) ## full Gamma matrix
#'
#' set.seed(123)
#' my_data <- rmpareto_tree(n, "HR", tree = my_graph, Gamma = Gamma_vec)
#' my_fit <- emst(my_data, p = NULL, method = "ML", cens = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
#'
emst <- function(
    data,
    p = NULL,
    method = c("vario", "ML", "chi"),
    cens = FALSE
){
  # Validate arguments
  method <- match.arg(method)

  # Check if you need to rescale data or not
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  # Estimate weight matrix
  if (method == "ML") {
    res <- ml_weight_matrix(data, cens = cens)
    weight_matrix <- res$llh_hr
    estimated_gamma <- res$est_gamma
  } else if (method == "vario") {
    weight_matrix <- emp_vario(data)
    estimated_gamma <- weight_matrix
  } else if (method == "chi") {
    weight_matrix <- - emp_chi(data)
    estimated_gamma <- chi2Gamma(- weight_matrix)
  }

  # Estimate tree
  graph.full <- igraph::make_full_graph(ncol(data))
  mst.tree <- igraph::mst(
    graph = graph.full,
    weights = weight_matrix[igraph::ends(graph.full, igraph::E(graph.full))],
    algorithm = "prim"
  )

  # Complete Gamma
  Gamma_comp <- complete_Gamma(estimated_gamma, mst.tree)

  # Return tree and completed Gamma
  return(list(
    graph = mst.tree,
    Gamma = Gamma_comp
  ))
}


#' Experimental: Fit graph using empirical Theta matrix
#'
#' Fits a graph to an empirical Gamma matrix by computing the corresponding
#' Theta matrix using [Gamma2Theta()] and greedily chooses `m` edges that
#' correspond to high absolute values in Theta.
#'
#' @param data The (standardized) data from which to compute Gamma
#' @param m The number of edges to add, defaults to the number of edges in a tree
#' @param Gamma_emp The empirical Gamma matrix
#' (can be `NULL` if `data` is given)
#'
#' @return A list containing an [igraph::graph] object and a fitted `Gamma` matrix
#'
#' @keywords internal
fit_graph_to_Theta <- function(data, m=NULL, Gamma_emp=NULL){

  if(is.null(Gamma_emp)){
    Gamma_emp <- emp_vario(data)
  }

  Theta_emp <- Gamma2Theta(Gamma_emp)

  d <- nrow(Gamma_emp)

  if(is.null(m)){
    m <- d-1
  }
  if(m < d-1){
    stop('m must be at least d-1!')
  }

  g_c <- igraph::make_full_graph(d)
  weights <- Theta_emp[igraph::as_edgelist(g_c)]
  g <- igraph::mst(g_c, weights)

  missing_edges <- igraph::as_edgelist(igraph::complementer(g))
  missing_weights <- Theta_emp[missing_edges]

  ord <- order(abs(missing_weights), decreasing = TRUE)

  new_edges <- missing_edges[ord[seq_len(m - (d-1))],]

  g <- igraph::add_edges(g, t(new_edges))

  Gamma_g <- complete_Gamma(Gamma_emp, g, N=100000, tol=1e-6, check_tol=100)

  return(list(
    graph = g,
    Gamma = Gamma_g
  ))
}