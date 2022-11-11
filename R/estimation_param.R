
#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on block graphs
#'
#' Fits the parameters of a multivariate Huesler--Reiss Pareto distribution using (censored) likelihood estimation.
#' See \insertCite{eng2019;textual}{graphicalExtremes} for details.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#'
#' @param graph Undirected graph object from `igraph` package with `d` vertices.
#'
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#'
#' @param method One of `"ML", "vario"`.
#' Default is `method = "ML"`.
#'
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#'
#' @return List consisting of:
#' \item{`graph`}{An [igraph::graph()] object.}
#' \item{`Gamma`}{Numeric \eqn{d\times d}{d x d} estimated variogram matrix \eqn{\Gamma}.}
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
fmpareto_graph_HR_OLD <- function(
  data,
  graph,
  p = NULL,
  method = c("ML", "vario"),
  cens = FALSE
){
  # Check arguments
  method <- match.arg(method)

  # Check graph
  d <- ncol(data)
  graph <- check_graph(graph, nVertices = d)

  if(is_decomposable_graph(graph)) {

    max_clique <- 2 #!!!

    if (max_clique > 3) {
      method <- "vario"
      warning("The maximal clique size is larger than 3. Forced to use empirical variogram.")
    }

    if (method == "ML") {
      fmpareto_graph_HR_decomposable(data, graph, p, cens)
    } else {
      fmpareto_graph_HR_general(data, graph, p)
    }

  } else {
    fmpareto_graph_HR_general(data, graph, p)
  }
}

fmpareto_graph_HR <- function(
  data,
  graph,
  p = NULL,
  method = c('ML', 'vario'),
  handleCliques = c('full', 'average', 'sequential'),
  ...
){
  # match args
  method <- match.arg(method)
  handleCliques <- match.arg(handleCliques)

  # standardize data
  if(!is.null(p)){
    data <- data2mpareto(data, p)
  }
  
  # call other functions, depending on handleCliques and method
  if(handleCliques == 'full'){
    # Estimate the entire parameter matrix at once
    if(method == 'ML'){
      Gamma <- fmpareto_HR_MLE(
        data,
        graph = graph,
        ...
      )
    } else{ # method == 'vario'
      variogram <- emp_vario(data)
      Gamma <- complete_Gamma(variogram, graph)
    }
  } else if(handleCliques == 'sequential'){
    # Estimate one clique after the other, fixing previously estimated entries.
    # Works only with MLE on decomposable graphs
    if(method != 'ML' || !is_decomposable_graph(graph)){
      stop('Arguments handleCliques="sequential" only works with decomposable graphs and method="ML"!')
    }
    Gamma <- fmpareto_graph_HR_clique_sequential(
      data,
      graph,
      ...
    )
  } else { # handleCliques == 'average'
    # Estimate each clique separately, taking the average on separators
    Gamma <- fmpareto_graph_HR_clique_average(
      data,
      graph,
      method,
      ...
    )
    if(!is_sym_cnd(Gamma)){
      stop('Averaging on the separators did not yield a valid variogram matrix!')
    }
  }
  return(Gamma)
}

fmpareto_graph_HR_clique_average <- function(
  data,
  graph,
  method = c('ML', 'vario'),
  ...
){
  method <- match.arg(method)
  
  graph <- check_graph(graph)
  cliques <- igraph::max_cliques(graph)
  # Todo: use mclapply here
  if(method == 'vario'){
    subGammas <- lapply(cliques, function(cli){
      data.cli <- mparetomargins(data, cli)
      emp_vario(data.cli)
    })
  } else {
    subGammas <- lapply(cliques, function(cli){
      data.cli <- mparetomargins(data, cli)
      tmp <- fmpareto_HR_MLE(
        data.cli,
        ...
      )
      if(is.null(tmp$Gamma)){
        stop('MLE did not find a valid Gamma for clique: ', paste0(cli, collapse = ','))
      }
      return(tmp$Gamma)
    })
  }
  Gamma_partial <- combine_clique_estimates_by_averaging(cliques, subGammas)
  Gamma <- complete_Gamma(Gamma_partial, graph)
  return(Gamma)
}


combine_clique_estimates_by_averaging <- function(cliques, subGammas){
  d <- do.call(max, cliques)
  Gamma <- matrix(0, d, d)
  overlaps <- matrix(0, d, d)
  for(i in seq_along(cliques)){
    cli <- cliques[[i]]
    Gamma[cli, cli] <- Gamma[cli, cli] + subGammas[[i]]
    overlaps[cli, cli] <- overlaps[cli, cli] + 1
  }
  Gamma[overlaps == 0] <- NA
  Gamma <- Gamma / overlaps
  return(Gamma)
}

fmpareto_graph_HR_clique_sequential <- function(
  data,
  graph,
  ...
){
  ## TODO: use mclapply
  
  # check inputs
  d <- ncol(data)
  graph <- check_graph(graph, graph_type = 'decomposable', nVertices = d)
  
  # Compute cliques and (recursive) separators
  layers <- get_cliques_and_separators(graph, sortIntoLayers = TRUE)
  
  # Estimate entries corresponding to separators/cliques
  Ghat <- matrix(NA, d, d)
  for(cliques in layers){
    # Compute estimate for each new clique/separator
    subGammas <- lapply(cliques, function(cli){
      # get margins data
      data.cli <- mparetomargins(data, cli)
      
      # find (already) fixed entries
      G.cli <- Ghat[cli, cli]
      par.cli <- Gamma2par(G.cli)
      fixParams.cli <- !is.na(par.cli)
      
      # get initial parameters that agree with the fixed ones (heuristic, close to empirical variogram):
      G0 <- emp_vario(data.cli)
      G1 <- replaceGammaSubMatrix(G0, G.cli)
      init.cli <- Gamma2par(G1)
      
      # estimate parameters
      opt <- fmpareto_HR_MLE(
        data = data.cli,
        init = init.cli,
        fixParams = fixParams.cli,
        ...
      )
      if(!opt$convergence){
        stop('MLE did not converge for clique: ', cli)
      }
      return(opt$Gamma)
    })

    # Fill newly computed entries in Ghat
    for(i in seq_along(cliques)){
      cli <- cliques[[i]]
      Ghat[cli, cli] <- subGammas[[i]]
    }
  }
  
  # Complete non-edge entries
  G_comp <- complete_Gamma(Ghat, graph)
  return(G_comp)
}


#' Parameter fitting for multivariate Huesler--Reiss Pareto distributions on non-decomposable graphs
#'
#' Works by computing an empirical Gamma matrix and then fitting this to the given `graph`.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#' @param graph Connected graph object from `igraph` package.
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
    message(paste0("The completed Gamm does not match the given graph.\n"))
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
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#' @param cens Logical. If true, then censored likelihood contributions are used for
#' components below the threshold. By default, `cens = FALSE`.
#' @param graph Decomposable graph object from `igraph` package.
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

  # loop through cliques and estimate Ghat:
  for(clique in cliques){
    # compute marginal pareto, on the nodes of the current clique
    data.cli <- mparetomargins(data, set_indices = clique)

    # find Ghat-entries that are already fixed by a previous clique/separator:
    G.cli <- Ghat[clique, clique]
    par.cli <- Gamma2par(G.cli)
    fixParams <- !is.na(par.cli)

    # get initial parameters, keeping the fixed ones:
    G.est0 <- emp_vario(data.cli)
    G.est <- replaceGammaSubMatrix(G.est0, G.cli)
    init.cli <- Gamma2par(G.est)

    # estimate parameters
    par.cli <- fmpareto_HR_MLE(
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
#' observations that is conditioned to be larger than the threshold `p`.
#' If `NULL` (default), then an average over all `k` is returned.
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
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
    data <- data2mpareto(data, p)
  }

  if (!is.null(k)) {
    G <- G.fun(k, data)

    if (any(is.na(G))) {
      warning(paste(
        "Produced NA matrix since there are no exceedances in the component k =",
        k
      ))
    }
  } else {

    # take the average
    row_averages <- rowMeans(sapply(1:d, FUN = function(i) {
      G.fun(i, data)
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
    data <- data2mpareto(data, p)
  }

  n <- nrow(data)
  d <- ncol(data)

  ind <- data > 1
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
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale. Else,
#' `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#'
#' @return Numeric. The empirical \eqn{d}-dimensional extremal correlation coefficient \eqn{\chi}
#' for the `data`.
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
    data <- data2mpareto(data, p)
  }



  rowmin <- apply(data, 1, min)
  chi <- mean(sapply(1:ncol(data), FUN = function(i) {
    mean(rowmin[which(data[, i] > 1)] > 1)
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
#' By default, `cens = FALSE`.
#'
#' @param p Numeric between 0 and 1 or `NULL`. If `NULL` (default),
#' it is assumed that the `data` are already on multivariate Pareto scale.
#'  Else, `p` is used as the probability in the function [data2mpareto]
#' to standardize the `data`.
#'
#' @return Numeric vector `c("loglik"=..., "aic"=..., "bic"=...)` with the evaluated
#' log-likelihood, AIC, and BIC values.
#'
#' @export
loglik_HR <- function(data, p = NULL, graph, Gamma, cens = FALSE){
  # Todo:
  # change log(n) to log(n_edges)
  # explain in help file how we compute bic, aic
  if (!is.null(p)) {
    data <- data2mpareto(data, p)
  }

  loglik <- logLH_HR(
    data = data,
    Gamma = Gamma,
    cens = cens
  )

  n <- nrow(data)
  n_edges <-  igraph::ecount(graph)

  aic <- 2 * n_edges - 2 * loglik

  bic <- log(n) * n_edges - 2 * loglik

  return(c("loglik" = loglik, "aic" = aic, "bic" = bic))
}

