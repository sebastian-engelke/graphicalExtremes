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
#'
#' @return List with as many elements as entries in `rholist`. Each element of
#' the list contains:
#' \describe{
#'   \item{`rho`}{The penalty coefficient.}
#'   \item{`graph`}{An [igraph::graph] object representing the fitted graph.}
#'   \item{`Gamma`}{A numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}
#' corresponding to the fitted graph.}
#' }
#'
#'
#' @export
eglasso <- function(Gamma, rholist=c(0.1, 0.15, 0.19, 0.205),
                    reg_method =  c("mb", "glasso"),
                    eps=0.5){

  # Check args
  reg_method <- match.arg(reg_method)

  # Set main variables
  r <- length(rholist)
  d <- ncol(Gamma)
  null.vote <- array(
    0,
    dim=c(d, d, length(rholist))) # votes for EXCLUDING the edge

  for(k in 1:d){
    Sk <- stats::cov2cor(Gamma2Sigma(Gamma=Gamma, k=k))
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
      (abs(gl.tmp$wi)==0)  ## change back to == 0 .. <=1e-4

  }
  adj.est <- (null.vote/(ncol(null.vote)-2)) < .49

  output <- list()
  for(j in 1:r) {
    rho <- rholist[j]
    est_graph <- igraph::graph_from_adjacency_matrix(adj.est[,,j],
                                                     mode="undirected",
                                                     diag=FALSE)
    Gamma <- tryCatch({
      complete_Gamma(graph = est_graph, Gamma = Gamma)
    },
    error = function(e){
      if (e$message == "The given graph is not connected."){
        NA
      } else {
        stop(e)
      }
    })

    output[[j]] <- list(rho = rho, graph = est_graph, Gamma = Gamma)

  }

  return(output)
}



#' Fitting extremal minimum spanning tree
#'
#' Fits an extremal minimum spanning tree, where the edge weights are:
#' \itemize{
#' \item negative maximized log-likelihoods of the bivariate Huesler--Reiss
#' distributions, if `method = "ML"`. See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#' \item empirical extremal variogram, if `method = "vario"`.
#' \item empirical extremal correlation, if `method = "chi"`.
#' }
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param method One of `"ML", "vario", "chi"`.
#' @param cens Logical. This argument is considered only if `method = "ML"`.
#' If `TRUE`, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{tree}: Graph object from \code{igraph} package. The fitted minimum spanning tree.
#' \item \code{Gamma}: Numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}
#' corresponding to the fitted minimum spanning tree.
#' }
#'
#' @examples
#' ## Fitting a 4-dimensional HR MST tree
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
#' my_fit <- mst_HR(my_data, p = NULL, cens = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
#'
emst <- function(data, p = NULL, method = c("ML", "vario", "chi"),
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
    tree = mst.tree,
    Gamma = complete_Gamma(graph = mst.tree, Gamma = estimated_gamma)
  ))
}



#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on block graphs
#'
#' Fits the parameters of a multivariate Huesler--Reiss Pareto distribution using (censored) likelihood estimation.
#' Fitting is done separately on the cliques of the block graph. If  \code{edges_to_add}
#' are provided, then these edges are added in a greedy search to the original \code{graph},
#' such that in each step the likelihood is improved maximally and the new graph stays in the
#' class of block graphs. See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
#' it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
#' \code{p} is used as the probability in the function \code{\link{data2mpareto}}
#' to standardize the \code{data}.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, \code{cens = FALSE}.
#' @param graph Graph object from \code{igraph} package. The \code{graph} must be an undirected block graph, i.e., a decomposable, connected
#' graph with singleton separator sets.
#' @param edges_to_add Numeric matrix \eqn{m\times 2}{m x 2}, where \eqn{m} is
#' the number of edges that are tried to be added in the greedy search.
#' By default, \code{edges_to_add = NULL}.
#'
#' @return List consisting of:
#' \itemize{
#' \item \code{graph}: Graph object from \code{igraph} package. If \code{edges_to_add} are provided,
#' then this is a list of the resulting graphs in each step of the greedy search.
#' \item \code{Gamma}: Numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}.
#' If \code{edges_to_add} are provided,
#' then this is a list of the estimated variogram matrices in each step of the greedy search.
#' \item \code{AIC}: (only if \code{edges_to_add} are provided) List of AIC values of the fitted models
#' in each step of the greedy search.
#' \item \code{edges_added}: (only if \code{edges_to_add} are provided) Numeric matrix \eqn{m'\times 2}{m' x 2}, where
#' the \eqn{m'\leq m}{m'<=m} rows contain the edges that were added in the greedy search.
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
#'   graph = my_graph,
#'   p = NULL, cens = FALSE, edges_to_add = edges_to_add
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

  n <- nrow(data)
  d <- ncol(data)

  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }


  ind <- data.std > 1
  ind_mat <- matrix(colSums(ind), byrow = TRUE, ncol = d, nrow = d)
  crossprod(ind, ind) / (1 / 2 * (ind_mat + t(ind_mat)))
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
loglik_HR <- function(data, p = NULL, graph, Gamma, cens){

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
