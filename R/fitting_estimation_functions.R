#' Fitting extremal graphical lasso
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


#' Learning extremal graph structure
#'
#' Fits an extremal graph structure using the neighborhood selection approach
#' (see \insertCite{meins2006;textual}{graphicalExtremes}) or graphical lasso
#' (see \insertCite{friedman2008;textual}{graphicalExtremes}).
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#'  number of observations and \eqn{d} is the dimension.
#'
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#'  it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#'  \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#'  to standardize the \code{data}.
#'
#' @param rholist Numeric vector of non-negative regularization parameters
#'  for the lasso.
#'  Default is `rholist = c(0.1, 0.15, 0.19, 0.205)`.
#'  For details see [glasso::glassopath].
#'
#' @param reg_method One of `"ns", "glasso"`, for neighborhood selection and
#'  graphical lasso, respectively.
#'  Default is `reg_method = "ns"`.
#'  For details see \insertCite{meins2006;textual}{graphicalExtremes},
#'  \insertCite{friedman2008;textual}{graphicalExtremes}.
#'
#' @param complete_Gamma Whether you want to try fto complete Gamma matrix.
#'  Default is `complete_Gamma = FALSE`.
#'
#' @return List made of:
#' \describe{
#'   \item{`graph`}{A list of [igraph::graph] objects representing the
#'   fitted graphs for each `rho` in `rholist`.}
#'   \item{`Gamma`}{A list of numeric \eqn{d\times d}{d x d} estimated
#'   variogram matrices \eqn{\Gamma} corresponding to the fitted graphs,
#'   for each `rho` in `rholist`. If `complete_Gamma = FALSE` or the
#'   underlying graph is not connected, it returns `NULL`.}
#'   \item{`rholist`}{The list of penalty coefficients.}
#'   \item{`graph_ic`}{A list of [igraph::graph] objects
#'   representing the optimal graph
#'   according to the `aic`, `bic`, and `mbic` information criteria.
#'   If `reg_method = "glasso"`, it returns a list of `NA`.}
#'   \item{`Gamma_ic`}{A list of numeric \eqn{d\times d}{d x d} estimated
#'   variogram matrices \eqn{\Gamma} corresponding
#'   to the `aic`, `bic`, and `mbic` information criteria.
#'   If `reg_method = "glasso"`, `complete_Gamma = FALSE`, or the underlying
#'   graph is not connected, it returns a list of `NA`.}
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
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  # Initialize variables
  Gamma <- emp_vario(data = data.std)
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
    }
    else if (reg_method == "ns") {
      idx_k <- which(data.std[, k] > 1)
      X <- log(data.std[idx_k, -k] / data.std[idx_k, k])
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
    }
    else {
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


#' Fitting extremal minimum spanning tree
#'
#' Fits an extremal minimum spanning tree, where the edge weights are:
#' \itemize{
#' \item negative maximized log-likelihoods of the bivariate Huesler--Reiss
#' distributions, if `method = "ML"`. See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#' \item empirical extremal variogram, if `method = "vario"`. See \insertCite{eng2020;textual}{graphicalExtremes} for details.
#' \item empirical extremal correlation, if `method = "chi"`. See \insertCite{eng2020;textual}{graphicalExtremes} for details.
#' }
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param method One of `"vario", "ML", "chi"`.
#' Default is `method = "vario"`.
#' @param cens Logical. This argument is considered only if `method = "ML"`.
#' If `TRUE`, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{graph}: An [igraph::graph()] object. The fitted minimum spanning tree.
#' \item \code{Gamma}: Numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}
#' corresponding to the fitted minimum spanning tree.
#' }
#'
#' @examples
#' ## Fitting a 4-dimensional HR minimum spanning tree
#'
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
#' my_data <- rmpareto_tree(n, "HR", tree = my_graph, par = Gamma_vec)
#' my_fit <- emst(my_data, p = NULL, method = "ML", cens = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
#'
emst <- function(data, p = NULL, method = c("vario", "ML", "chi"),
                 cens = FALSE) {

  # Validate arguments
  method <- match.arg(method)

  # Check if you need to rescale data or not
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  # Estimate weight matrix
  if (method == "ML") {
    res <- ml_weight_matrix(data = data.std, cens = cens)
    weight_matrix <- res$llh_hr
    estimated_gamma <- res$est_gamma

  } else if (method == "vario") {
    weight_matrix <- emp_vario(data = data.std)
    estimated_gamma <- weight_matrix

  } else if (method == "chi") {
    weight_matrix <- - emp_chi(data = data.std)
    estimated_gamma <- chi2Gamma(- weight_matrix)

  }

  # Estimate tree
  graph.full <- igraph::make_full_graph(ncol(data.std))
  mst.tree <- igraph::mst(
    graph = graph.full,
    weights = weight_matrix[igraph::ends(graph.full, igraph::E(graph.full))],
    algorithm = "prim"
  )

  # Set graphical parameters
  mst.tree <- mst.tree

  # Return tree and completed Gamma
  return(list(
    graph = mst.tree,
    Gamma = complete_Gamma(graph = mst.tree, Gamma = estimated_gamma)
  ))
}



#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on block graphs
#'
#' Fits the parameters of a multivariate Huesler--Reiss Pareto distribution using (censored) likelihood estimation.
#' See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#'
#' @param graph Undirected graph object from \code{igraph} package.
#'
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#'
#' @param method One of `"ML", "vario"`.
#' Default is `method = "ML"`.
#'
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{graph}: An [igraph::graph()] object.
#' \item \code{Gamma}: Numeric \eqn{d\times d}{d x d} estimated variogram matrix
#' \eqn{\Gamma}.
#' }
#'
#' @examples
#' ## Fitting a 4-dimensional HR distribution
#'
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
#' edges_to_add <- rbind(c(1, 3), c(1, 4), c(3, 4))
#'
#' set.seed(123)
#' my_data <- rmpareto_tree(n, "HR", tree = my_graph, par = Gamma_vec)
#' my_fit <- fmpareto_graph_HR(my_data,
#'   graph = my_graph, method = "ML",
#'   p = NULL, cens = FALSE
#' )
#' @references
#'  \insertAllCited{}
#' @export
fmpareto_graph_HR <- function(data, graph, p = NULL, method = c("ML", "vario"),
                              cens = FALSE){

  # Check arguments
  method <- match.arg(method)

  # Check graph
  d <- ncol(data)
  graph <- check_graph(graph, nVertices = d)

  # Infer graph type
  is_decomp <- igraph::is.chordal(graph)

  if (is_decomp$chordal) {

    max_clique <- 2 #!!!

    if (max_clique > 3) {
      method <- "vario"
      warning(paste0("The maximal clique size is larger than 3.",
                     " Forced to use empirical variogram."), call. = FALSE)
    }

    if (method == "ML") {
      fmpareto_graph_HR_decomposable(data, graph, p, cens)
    } else {
      fmpareto_graph_HR_general(data, graph, p) #!!! to optimize
    }

  } else {
    fmpareto_graph_HR_general(data, graph, p)
  }
}




#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on non-decomposable graphs
#'
#' Works by computing an empirical Gamma matrix and then fitting this to the given `graph`.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param graph Connected graph object from \code{igraph} package.
#' @param ... Further arguments passed to [complete_gamma_general()] if `graph`
#' is not decomposable
#'
#'
#' @export
fmpareto_graph_HR_general <- function(data, graph, p = NULL, ...) {
  # number of data columns (and expected vertices in graph)
  d <- NCOL(data)

  # make sure graph is valid
  graph <- check_graph(graph, graph_type = 'general', nVertices = d)

  # compute the empirical variogram
  Gamma_emp <- emp_vario(data, p = p)

  # fit the empirical variogram to the graph
  Gamma_graph <- complete_Gamma(Gamma_emp, graph, allowed_graph_type = 'general', ...)

  # check that completed Gamma matches the given graph
  completed_graph <- Gamma2graph(Gamma_graph, to_plot = FALSE)

  if (!(graphs_equal(completed_graph, graph))) {
    message(paste0("The completed Gamma", " does not match the given graph.\n"))
  }

  return(list(
    Gamma = Gamma_graph,
    graph = graph
  ))
}

#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on decomposable graphs
#'
#' Similar to [fmpareto_graph_HR()]. Differences are:
#' * Works with any decomposable graph, not just block graphs
#' * Does not support `edges_to_add`
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, \code{cens = FALSE}.
#' @param graph Decomposable graph object from \code{igraph} package.
#'
#' @export
fmpareto_graph_HR_decomposable <- function(data, graph, p = NULL, cens = FALSE) {
  # number of data columns (and expected vertices in graph)
  d <- NCOL(data)

  # make sure graph is valid
  graph <- check_graph(graph, graph_type = 'decomposable', nVertices = d)

  # rescale data if necessary:
  if(!is.null(p)){
    data <- data2mpareto(data, p)
  }

  # compute cliques:
  cliques <- get_cliques_and_separators(graph)

  # initialize variables:
  Ghat <- matrix(NA, d, d)
  fixParams <- logical(d)

  # loop through cliques and estimate Ghat:
  for(clique in cliques){
    # compute marginal pareto, on the nodes of the current clique
    data.cli <- mparetomargins(data = data, set_indices = clique)

    # find Ghat-entries that are already fixed by a previous clique/separator:
    G.cli <- Ghat[clique, clique]
    par.cli <- Gamma2par(G.cli)
    fixParams <- !is.na(par.cli)

    # get initial parameters, keeping the fixed ones:
    G.est0 <- emp_vario(data = data.cli)
    G.est <- replaceGammaSubMatrix(G.est0, G.cli)
    init.cli <- Gamma2par(G.est)

    # estimate parameters
    par.cli <- fmpareto_HR(
      data = data.cli,
      init = init.cli,
      fixParams = fixParams,
      cens = cens
    )$par

    # update Ghat
    G.cli <- par2Gamma(par.cli)
    Ghat[clique, clique] <- G.cli
  }

  # fill entries that don't correspond to edges:
  Ghat <- complete_Gamma(Ghat, allowed_graph_type = 'decomposable')

  # check that completed Gamma matches the given graph
  completed_graph <- Gamma2graph(Ghat, to_plot = FALSE)

  if (!(graphs_equal(completed_graph, graph))) {
    message(paste0("The completed Gamma", " does not match the given graph.\n"))
  }

  return(list(graph = graph, Gamma = Ghat))
}



#' Estimation of the variogram matrix \eqn{\Gamma} of the Huesler--Reiss distribution
#'
#' Estimates the variogram of the Huesler--Reiss distribution empirically.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param k Integer between 1 and \eqn{d}. Component of the multivariate
#' observations that is conditioned to be larger than the threshold \code{p}.
#' If \code{NULL} (default), then an average over all \code{k} is returned.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d}. The estimated
#' variogram of the Huesler--Reiss distribution.
#' @export
emp_vario <- function(data, k = NULL, p = NULL) {

  # helper ####
  G.fun <- function(i, data) {
    idx <- which(data[, i] > 1)
    if (length(idx) > 1) {
      xx <- Sigma2Gamma(stats::cov(log(data[idx, ])), full = TRUE)
    } else {
      xx <- matrix(NA, d, d)
    }
    return(xx)
  }

  # body ####
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  d <- ncol(data)
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  if (!is.null(k)) {
    G <- G.fun(k, data.std)

    if (any(is.na(G))) {
      warning(paste(
        "Produced NA matrix since there are no exceedances in the component k =",
        k
      ))
    }
  } else {

    # take the average
    row_averages <- rowMeans(sapply(1:d, FUN = function(i) {
      G.fun(i, data.std)
    }), na.rm = TRUE)
    G <- matrix(row_averages, nrow = d, ncol = d)
  }

  return(G)
}

#' @export
emp_vario_pairwise <- function(data, k = NULL, p = NULL, verbose = FALSE){
  vcat <- function(...){
    if(verbose){
      cat(...)
    }
  }
  d <- ncol(data)
  vario <- matrix(0, d, d)
  vario[] <- NA
  pct <- -Inf
  vcat('Computing emp_vario_pairwise, reporting progress:\n')
  for(i in seq_len(d)){
    vario[i,i] <- 0
    if(i + 1 > d){
      next
    }
    for(j in (i+1):d){
      data_ij <- data2mpareto(data[,c(i,j)], p, na.rm=TRUE)
      vario_ij <- emp_vario(data_ij, p=NULL)
      vario[i,j] <- vario_ij[1,2]
      vario[j,i] <- vario_ij[1,2]
      pct1 <- floor(100 * sum(!is.na(vario)) / length(vario))
      if(pct1 > pct){
        pct <- pct1
        vcat(pct, '% ', sep='')
      }
    }
  }
  cat('\nDone.\n')
  return(vario)
}


#' Empirical estimation of extremal correlation matrix \eqn{\chi}
#'
#' Estimates empirically the matrix of bivariate extremal correlation coefficients \eqn{\chi}.
#'
#' @inheritParams emp_chi_multdim
#'
#' @return Numeric matrix \eqn{d\times d}{d x d}. The matrix contains the
#' bivariate extremal coefficients \eqn{\chi_{ij}}, for \eqn{i, j = 1, ..., d}.
#' @examples
#' n <- 100
#' d <- 4
#' p <- .8
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = Gamma)
#' emp_chi(my_data, p)
#' @export
emp_chi <- function(data, p = NULL) {
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  n <- nrow(data.std)
  d <- ncol(data.std)

  ind <- data.std > 1
  ind_mat <- matrix(colSums(ind), byrow = TRUE, ncol = d, nrow = d)
  crossprod(ind, ind) / (1 / 2 * (ind_mat + t(ind_mat)))
}

#' @export
emp_chi_pairwise <- function(data, p = NULL, verbose=FALSE){
  vcat <- function(...){
    if(verbose){
      cat(...)
    }
  }
  d <- ncol(data)
  chi <- matrix(0, d, d)
  chi[] <- NA
  pct <- -Inf
  vcat('Computing emp_chi_pairwise, reporting progress:\n')
  for(i in seq_len(d)){
    chi[i,i] <- 1
    if(i + 1 > d){
      next
    }
    for(j in (i+1):d){
      data_ij <- data2mpareto(data[,c(i,j)], p, na.rm=TRUE)
      chi_ij <- emp_chi(data_ij, NULL)
      chi[i,j] <- chi_ij[1,2]
      chi[j,i] <- chi_ij[1,2]
      pct1 <- floor(100 * sum(!is.na(chi)) / length(chi))
      if(pct1 > pct){
        pct <- pct1
        vcat(pct, '% ', sep='')
      }
    }
  }
  cat('\nDone.\n')
  return(chi)
}



#' Empirical estimation of extremal correlation \eqn{\chi}
#'
#' Estimates the \eqn{d}-dimensional extremal correlation coefficient \eqn{\chi} empirically.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#'
#' @return Numeric. The empirical \eqn{d}-dimensional extremal correlation coefficient \eqn{\chi}
#' for the \code{data}.
#' @examples
#' n <- 100
#' d <- 2
#' p <- .8
#' G <- cbind(
#'   c(0, 1.5),
#'   c(1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = G)
#' emp_chi_multdim(my_data, p)
#' @export
emp_chi_multdim <- function(data, p = NULL) {
  if (!is.matrix(data)) {
    stop("The data should be a matrix")
  }
  if (ncol(data) <= 1) {
    stop("The data should be a matrix with at least two columns.")
  }

  d <- ncol(data)
  data <- stats::na.omit(data)

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }



  rowmin <- apply(data.std, 1, min)
  chi <- mean(sapply(1:ncol(data.std), FUN = function(i) {
    mean(rowmin[which(data.std[, i] > 1)] > 1)
  }), na.rm = TRUE)
  return(chi)
}



#' Compute Huesler--Reiss log-likelihood, AIC, and BIC
#'
#' Computes (censored) Huesler--Reiss log-likelihood, AIC, and BIC values.
#'
#' @param data Numeric matrix \eqn{n\times d}{n x d}. It contains
#' observations following a multivariate HR Pareto distribution.
#'
#' @param graph An [igraph::graph] object. The `graph` must be undirected and
#' connected
#'
#' @param Gamma Numeric matrix \eqn{n\times d}{n x d}.
#' It represents a variogram matrix \eqn{\Gamma}.
#'
#' @param cens Boolean. If true, then censored log-likelihood is computed.
#' By default, \code{cens = FALSE}.
#'
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale.
#'  Else, `p` is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#'
#' @return Numeric vector `c("loglik"=..., "aic"=..., "bic"=...)` with the evaluated
#' log-likelihood, AIC, and BIC values.
#'
#' @export
loglik_HR <- function(data, p = NULL, graph, Gamma, cens = FALSE){

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  loglik <- logLH_HR(
    data = data.std,
    Gamma = Gamma,
    cens = cens
  )

  n <- nrow(data)
  n_edges <-  igraph::ecount(graph)

  aic <- 2 * n_edges - 2 * loglik

  bic <- log(n) * n_edges - 2 * loglik

  c("loglik" = loglik, "aic" = aic, "bic" = bic)

}

# Todo:
# change log(n) above to log(n_edges)
# explain in help file how we compute bic, aic
