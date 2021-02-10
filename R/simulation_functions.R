#' Sampling of a multivariate Pareto distribution
#'
#' Simulates exact samples of a multivariate Pareto distribution.
#'
#' @param n Number of simulations.
#' @param model The parametric model type; one of:
#' \itemize{
#' \item \code{HR} (default),
#' \item \code{logistic},
#' \item \code{neglogistic},
#' \item \code{dirichlet}.
#' }
#' @param d Dimension of the multivariate Pareto
#' distribution.
#' @param par Respective parameter for the given \code{model}, that is,
#' \itemize{
#' \item \eqn{\Gamma}, numeric \eqn{d \times d}{d x d} variogram matrix,
#' if \code{model = HR}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, if \code{model = logistic}.
#' \item \eqn{\theta > 0}, if \code{model = neglogistic}.
#' \item \eqn{\alpha}, numeric vector of size \code{d} with positive entries,
#' if \code{model = dirichlet}.
#' }
#'
#' @return
#' Numeric matrix of size \eqn{n \times d}{n x d} of simulations of the
#' multivariate Pareto distribution.
#'
#' @details
#' The simulation follows the algorithm in \insertCite{eng2019;textual}{graphicalExtremes}.
#' For details on the parameters of the Huesler--Reiss, logistic
#' and negative logistic distributions see \insertCite{dom2016;textual}{graphicalExtremes}, and for the Dirichlet
#' distribution see \insertCite{coles1991modelling;textual}{graphicalExtremes}.
#'
#' @examples
#' ## A 4-dimensional HR distribution
#' n <- 10
#' d <- 4
#' G <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' rmpareto(n, "HR", d = d, par = G)
#'
#' ## A 3-dimensional logistic distribution
#' n <- 10
#' d <- 3
#' theta <- .6
#' rmpareto(n, "logistic", d, par = theta)
#'
#' ## A 5-dimensional negative logistic distribution
#' n <- 10
#' d <- 5
#' theta <- 1.5
#' rmpareto(n, "neglogistic", d, par = theta)
#'
#' ## A 4-dimensional Dirichlet distribution
#' n <- 10
#' d <- 4
#' alpha <- c(.8, 1, .5, 2)
#' rmpareto(n, "dirichlet", d, par = alpha)
#' @references
#'  \insertAllCited{}
#'
#' @export
rmpareto <- function(n,
                     model = c("HR", "logistic", "neglogistic", "dirichlet")[1],
                     d, par) {

  # methods
  model_nms <- c("HR", "logistic", "neglogistic", "dirichlet")

  # check arguments ####
  if (d != round(d) | d < 1) {
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1) {
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)) {
    stop(paste("The model must be one of ", paste(model_nms, collapse = ", "),
      ".",
      sep = ""
    ))
  }

  if (model == "HR") {
    if (!is.matrix(par)) {
      stop("The argument par must be a matrix, when model = HR.")
    }

    if (NROW(par) != d | NCOL(par) != d) {
      stop("The argument par must be a d x d matrix, when model = HR.")
    }
  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12) {
      stop(paste(
        "The argument par must be scalar between 1e-12 and 1 - 1e-12,",
        "when model = logistic."
      ))
    }
  } else if (model == "neglogistic") {
    if (par <= 1e-12) {
      stop(paste(
        "The argument par must be scalar greater than 1e-12,",
        "when model = neglogistic."
      ))
    }
  } else if (model == "dirichlet") {
    if (length(par) != d) {
      stop(paste(
        "The argument par must be a vector with d elements,",
        "when model = dirichlet."
      ))
    }

    if (any(par <= 1e-12)) {
      stop(paste(
        "The elements of par must be greater than 1e-12,",
        "when model = dirichlet."
      ))
    }
  }

  # prepare arguments ####
  if (model == "HR") {
    Gamma <- par

    # compute cholesky decomposition
    cov.mat <- Gamma2Sigma(Gamma, k = 1, full = FALSE)
    chol_mat <- matrix(0, d, d)

    result <- tryCatch(
      {
        chol_mat[-1, -1] <- chol(cov.mat)
      },
      error = function(e) {
        stop(paste(
          "The covariance matrix associated to Gamma",
          "cannot be factorized with Cholesky."
        ))
      }
    )

    # compute trend (matrix where each row is one variable)
    trend <- t(sapply(1:d, function(k) {
      sapply(1:d, function(j) {
        Gamma[j, k] / 2
      })
    }))
  } else if (model == "logistic") {
    theta <- par
  } else if (model == "neglogistic") {
    theta <- par
  } else if (model == "dirichlet") {
    alpha <- par
  }

  # function body ####
  counter <- 0
  res <- numeric(0)
  n.total <- 0
  while (n.total < n) {
    counter <- counter + 1
    shift <- sample(1:d, n, replace = TRUE)
    for (k in 1:d) {
      n.k <- sum(shift == k)

      if (n.k > 0) {
        proc <-
          switch(model,
            "HR" =
              simu_px_HR(
                n = n.k, idx = k, d = d, trend = trend[k, ],
                chol_mat = chol_mat
              ),
            "logistic" =
              simu_px_logistic(n = n.k, idx = k, d = d, theta = theta),
            "neglogistic" =
              simu_px_neglogistic(n = n.k, idx = k, d = d, theta = theta),
            "dirichlet" =
              simu_px_dirichlet(n = n.k, idx = k, d = d, alpha = alpha)
          )

        if (any(dim(proc) != c(n.k, d))) {
          stop("The generated sample has wrong size.")
        }

        proc <- proc / rowSums(proc) / (1 - stats::runif(NROW(proc)))
        idx.sim <- which(apply(proc, 1, max) > 1)
        res <- rbind(res, proc[idx.sim, ])
        n.total <- NROW(res)
      }
    }
  }

  return(res[sample(1:NROW(res), n, replace = FALSE), ])
}



#' Sampling of a multivariate Pareto distribution on a tree
#'
#' Simulates exact samples of a multivariate Pareto distribution that
#' is an extremal graphical model on a tree as defined in \insertCite{eng2019;textual}{graphicalExtremes}.
#'
#' @param n Number of simulations.
#' @param model The parametric model type; one of:
#' \itemize{
#' \item \code{HR} (default),
#' \item \code{logistic},
#' \item \code{dirichlet}.
#' }
#' @param tree Graph object from \code{igraph} package.
#' This object must be a tree, i.e., an
#' undirected graph that is connected and has no cycles.
#' @param par Respective parameter for the given \code{model}, that is,
#' \itemize{
#' \item \eqn{\Gamma}, numeric \eqn{d \times d}{d x d} variogram matrix,
#' where only the entries corresponding to the edges of the \code{tree} are used,
#' if \code{model = HR}. Alternatively, can be a vector of
#' length \eqn{d - 1} containing the entries of the variogram corresponding
#' to the edges of the given \code{tree}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, vector of length \eqn{d - 1}
#' containing the logistic parameters corresponding
#' to the edges of the given \code{tree}, if \code{model = logistic}.
#' \item a matrix of size \eqn{(d - 1) \times 2}{(d - 1) x 2}, where the rows
#' contain the parameters vectors \eqn{\alpha} of size 2 with positve entries
#' for each of the edges in \code{tree}, if \code{model = dirichlet}.
#' }
#'
#' @return
#' Numeric matrix of size \eqn{n \times d}{n x d} of simulations of the
#' multivariate Pareto distribution.
#'
#' @details
#' The simulation follows the algorithm in \insertCite{eng2019;textual}{graphicalExtremes}.
#' For details on the parameters of the Huesler--Reiss, logistic
#' and negative logistic distributions see \insertCite{dom2016;textual}{graphicalExtremes}, and for the Dirichlet
#' distribution see \insertCite{coles1991modelling;textual}{graphicalExtremes}.
#'
#' @examples
#' ## A 4-dimensional HR tree model
#'
#' my_tree <- igraph::graph_from_adjacency_matrix(rbind(
#'   c(0, 1, 0, 0),
#'   c(1, 0, 1, 1),
#'   c(0, 1, 0, 0),
#'   c(0, 1, 0, 0)
#' ),
#' mode = "undirected"
#' )
#' n <- 10
#' Gamma_vec <- c(.5, 1.4, .8)
#' set.seed(123)
#' rmpareto_tree(n, "HR", tree = my_tree, par = Gamma_vec)
#'
#' ## A 4-dimensional Dirichlet model with asymmetric edge distributions
#'
#' alpha <- cbind(c(.2, 1, .5), c(1.5, .6, .8))
#' rmpareto_tree(n, model = "dirichlet", tree = my_tree, par = alpha)
#' @references
#'  \insertAllCited{}
#'
#' @export
rmpareto_tree <- function(n, model = c("HR", "logistic", "dirichlet")[1],
                          tree, par) {

  # methods
  model_nms <- c("HR", "logistic", "dirichlet")

  # graph theory objects ####
  # check if it is directed
  if (igraph::is_directed(tree)) {
    warning("The given tree is directed. Converted to undirected.")
    tree <- igraph::as.undirected(tree)
  }

  # set graph theory objects
  adj <- as.matrix(igraph::as_adj(tree))
  d <- NROW(adj)
  e <- igraph::ecount(tree)
  ends.mat <- igraph::ends(tree, igraph::E(tree))

  # check if it is tree
  is_connected <- igraph::is_connected(tree)
  is_tree <- is_connected & (e == d - 1)

  if (!is_tree) {
    stop("The given graph is not a tree.")
  }


  # check arguments ####
  if (d != round(d) | d < 1) {
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1) {
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)) {
    stop(paste("The model must be one of ", paste(model_nms, collapse = ", "),
      ".",
      sep = ""
    ))
  }

  if (model == "HR") {
    if (!is.matrix(par)) {
      if (length(par) != (d - 1)) {
        stop(paste(
          "The argument par must be a d x d matrix,",
          "or a vector with d - 1 elements, when model = HR."
        ))
      }
    } else {
      if (NROW(par) != d | NCOL(par) != d) {
        stop(paste(
          "The argument par must be a d x d matrix,",
          "or a vector with d elements, when model = HR."
        ))
      }
    }
  } else if (model == "logistic") {
    if (any(par <= 1e-12) | any(par >= 1 - 1e-12)) {
      stop(paste(
        "The elements of par must be",
        "between 1e-12 and 1 - 1e-12,",
        "when model = logistic."
      ))
    }

    if (length(par) == 1) {
      par <- rep(par, d - 1)
      warning(paste(
        "The argument par was a scalar.",
        "Converted to a vector of size d - 1 with",
        "all same entries."
      ))
    }
    if (length(par) != d - 1) {
      stop(paste(
        "The argument par must have d - 1 elements,",
        "when model = logistic."
      ))
    }
  } else if (model == "dirichlet") {
    if (NROW(par) != d - 1 | NCOL(par) != 2) {
      stop(paste(
        "The argument par must be a (d-1) x 2 ,",
        "when model = dirichlet."
      ))
    }
    if (any(par <= 1e-12)) {
      stop(paste(
        "The elements of par must be greater than 1e-12,",
        "when model = dirichlet."
      ))
    }
  }

  # prepare arguments ####
  if (model == "HR") {
    if (is.matrix(par)) {
      par.vec <- par[ends.mat]
    } else {
      par.vec <- par
    }
  } else if (model == "logistic") {
    theta <- par
  } else if (model == "dirichlet") {
    alpha.mat <- par
  }

  # function body ####
  # Define a matrix A[[k]] choosing the paths from k to other vertices
  idx.e <- matrix(0, nrow = d, ncol = d)
  idx.e[ends.mat] <- 1:e
  idx.e <- idx.e + t(idx.e)

  # e.start[[k]][h] gives the index (1 or 2) of the starting node in the h edge
  # in the tree rooted at k
  A <- e.start <- e.end <- list()
  for (k in 1:d) {
    A[[k]] <- matrix(0, nrow = d, ncol = e)
    e.start[[k]] <- e.end[[k]] <- numeric(e)
    short.paths <- igraph::shortest_paths(tree, from = k, to = 1:d)
    for (h in 1:d) {
      path <- short.paths$vpath[[h]]
      idx.tmp <- idx.e[cbind(path[-length(path)], path[-1])]
      A[[k]][h, idx.tmp] <- 1
      e.start[[k]][idx.tmp] <-
        apply(ends.mat[idx.tmp, ] ==
          matrix(path[-length(path)], nrow = length(idx.tmp), ncol = 2),
        MARGIN = 1, FUN = function(x) which(x == TRUE)
        )
      e.end[[k]][idx.tmp] <-
        apply(ends.mat[idx.tmp, ] ==
          matrix(path[-1], nrow = length(idx.tmp), ncol = 2),
        MARGIN = 1, FUN = function(x) which(x == TRUE)
        )
    }
  }

  counter <- 0
  res <- numeric(0)
  n.total <- 0
  while (n.total < n) {
    counter <- counter + 1
    shift <- sample(1:d, n, replace = TRUE)
    for (k in 1:d) {
      n.k <- sum(shift == k)
      if (n.k > 0) {
        proc <-
          switch(model,
            "HR" =
              simu_px_tree_HR(n = n.k, Gamma_vec = par.vec, A = A[[k]]),
            "logistic" =
              simu_px_tree_logistic(
                n = n.k,
                theta = theta, A = A[[k]]
              ),
            "dirichlet" =
              simu_px_tree_dirichlet(
                n = n.k,
                alpha.start =
                  alpha.mat[cbind(1:e, e.start[[k]])],
                alpha.end =
                  alpha.mat[cbind(1:e, e.end[[k]])],
                A = A[[k]]
              )
          )

        if (any(dim(proc) != c(n.k, d))) {
          stop("The generated sample has wrong size.")
        }

        proc <- proc / rowSums(proc) / (1 - stats::runif(NROW(proc)))
        idx.sim <- which(apply(proc, 1, max) > 1)
        res <- rbind(res, proc[idx.sim, ])
        n.total <- NROW(res)
      }
    }
  }

  return(res[sample(1:NROW(res), n, replace = FALSE), ])
}



#' Sampling of a multivariate max-stable distribution
#'
#' Simulates exact samples of a multivariate max-stable distribution.
#'
#' @param n Number of simulations.
#' @param model The parametric model type; one of:
#' \itemize{
#' \item \code{HR} (default),
#' \item \code{logistic},
#' \item \code{neglogistic},
#' \item \code{dirichlet}.
#' }
#' @param d Dimension of the multivariate Pareto
#' distribution.
#' @param par Respective parameter for the given \code{model}, that is,
#' \itemize{
#' \item \eqn{\Gamma}, numeric \eqn{d \times d}{d x d} variogram matrix,
#' if \code{model = HR}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, if \code{model = logistic}.
#' \item \eqn{\theta > 0}, if \code{model = neglogistic}.
#' \item \eqn{\alpha}, numeric vector of size \code{d} with positive entries,
#' if \code{model = dirichlet}.
#' }
#'
#' @return
#' Numeric matrix of size \eqn{n \times d}{n x d} of simulations of the
#' multivariate max-stable distribution.
#'
#' @details
#' The simulation follows the extremal function algorithm in \insertCite{dom2016;textual}{graphicalExtremes}.
#' For details on the parameters of the Huesler--Reiss, logistic
#' and negative logistic distributions see \insertCite{dom2016;textual}{graphicalExtremes}, and for the Dirichlet
#' distribution see \insertCite{coles1991modelling;textual}{graphicalExtremes}.
#'
#' @examples
#' ## A 4-dimensional HR distribution
#' n <- 10
#' d <- 4
#' G <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' rmstable(n, "HR", d = d, par = G)
#'
#' ## A 3-dimensional logistic distribution
#' n <- 10
#' d <- 3
#' theta <- .6
#' rmstable(n, "logistic", d, par = theta)
#'
#' ## A 5-dimensional negative logistic distribution
#' n <- 10
#' d <- 5
#' theta <- 1.5
#' rmstable(n, "neglogistic", d, par = theta)
#'
#' ## A 4-dimensional Dirichlet distribution
#' n <- 10
#' d <- 4
#' alpha <- c(.8, 1, .5, 2)
#' rmstable(n, "dirichlet", d, par = alpha)
#' @references
#'  \insertAllCited{}
#'
#' @export
rmstable <- function(n,
                     model = c("HR", "logistic", "neglogistic", "dirichlet")[1],
                     d, par) {

  # methods
  model_nms <- c("HR", "logistic", "neglogistic", "dirichlet")

  # check arguments ####
  if (d != round(d) | d < 1) {
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1) {
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)) {
    stop(paste("The model must be one of ", paste(model_nms, collapse = ", "),
      ".",
      sep = ""
    ))
  }

  if (model == "HR") {
    if (!is.matrix(par)) {
      stop("The argument par must be a matrix, when model = HR.")
    }

    if (NROW(par) != d | NCOL(par) != d) {
      stop("The argument par must be a d x d matrix, when model = HR.")
    }
  } else if (model == "logistic") {
    if (length(par) != 1 | par <= 1e-12 | par >= 1 - 1e-12) {
      stop(paste(
        "The argument par must be scalar between 1e-12 and 1 - 1e-12,",
        "when model = logistic."
      ))
    }
  } else if (model == "neglogistic") {
    if (par <= 1e-12) {
      stop(paste(
        "The argument par must be scalar greater than 1e-12,",
        "when model = neglogistic."
      ))
    }
  } else if (model == "dirichlet") {
    if (length(par) != d) {
      stop(paste(
        "The argument par must be a vector with d elements,",
        "when model = dirichlet."
      ))
    }

    if (any(par <= 1e-12)) {
      stop(paste(
        "The elements of par must be greater than 1e-12,",
        "when model = dirichlet."
      ))
    }
  }

  # prepare arguments ####
  if (model == "HR") {
    Gamma <- par

    # compute cholesky decomposition
    cov.mat <- Gamma2Sigma(Gamma, k = 1, full = FALSE)
    chol_mat <- matrix(0, d, d)

    result <- tryCatch(
      {
        chol_mat[-1, -1] <- chol(cov.mat)
      },
      error = function(e) {
        stop(paste(
          "The covariance matrix associated to Gamma",
          "cannot be factorized with Cholesky."
        ))
      }
    )

    # compute trend (matrix where each row is one variable)
    trend <- t(sapply(1:d, function(k) {
      sapply(1:d, function(j) {
        Gamma[j, k] / 2
      })
    }))
  } else if (model == "logistic") {
    theta <- par
  } else if (model == "neglogistic") {
    theta <- par
  } else if (model == "dirichlet") {
    alpha <- par
  }

  # function body ####
  counter <- rep(0, times = n)
  res <- matrix(0, nrow = n, ncol = d)
  for (k in 1:d) {
    poisson <- stats::rexp(n)

    while (any(1 / poisson > res[, k])) {
      ind <- (1 / poisson > res[, k])
      n.ind <- sum(ind)
      idx <- (1:n)[ind]
      counter[ind] <- counter[ind] + 1
      proc <-
        switch(model,
          "HR" =
            simu_px_HR(
              n = n.ind, idx = k, d = d, trend = trend[k, ],
              chol_mat = chol_mat
            ),
          "logistic" =
            simu_px_logistic(n = n.ind, idx = k, d = d, theta = theta),
          "neglogistic" =
            simu_px_neglogistic(n = n.ind, idx = k, d = d, theta = theta),
          "dirichlet" =
            simu_px_dirichlet(n = n.ind, idx = k, d = d, alpha = alpha)
        )

      if (any(dim(proc) != c(n.ind, d))) {
        stop("The generated sample has wrong size.")
      }

      if (k == 1) {
        ind.upd <- rep(TRUE, times = n.ind)
      } else {
        ind.upd <- sapply(1:n.ind, function(i) {
          all(1 / poisson[idx[i]] * proc[i, 1:(k - 1)] <=
            res[idx[i], 1:(k - 1)])
        })
      }
      if (any(ind.upd)) {
        idx.upd <- idx[ind.upd]
        res[idx.upd, ] <- pmax(res[idx.upd, ], 1 / poisson[idx.upd] *
          proc[ind.upd, ])
      }
      poisson[ind] <- poisson[ind] + stats::rexp(n.ind)
    }
  }

  return(res[sample(1:NROW(res), n, replace = FALSE), ])
}



#' Sampling of a multivariate max-stable distribution on a tree
#'
#' Simulates exact samples of a multivariate max-stable distribution that
#' is an extremal graphical model on a tree as defined in \insertCite{eng2019;textual}{graphicalExtremes}.
#'
#' @param n Number of simulations.
#' @param model The parametric model type; one of:
#' \itemize{
#' \item \code{HR} (default),
#' \item \code{logistic},
#' \item \code{dirichlet}.
#' }
#' @param tree Graph object from \code{igraph} package.
#' This object must be a tree, i.e., an
#' undirected graph that is connected and has no cycles.
#' @param par Respective parameter for the given \code{model}, that is,
#' \itemize{
#' \item \eqn{\Gamma}, numeric \eqn{d \times d}{d x d} variogram matrix,
#' where only the entries corresponding to the edges of the \code{tree} are used,
#' if \code{model = HR}. Alternatively, can be a vector of
#' length \eqn{d - 1} containing the entries of the variogram corresponding
#' to the edges of the given \code{tree}.
#' \item \eqn{\theta \in (0, 1)}{0 < \theta < 1}, vector of length \eqn{d - 1}
#' containing the logistic parameters corresponding
#' to the edges of the given \code{tree}, if \code{model = logistic}.
#' \item a matrix of size \eqn{(d - 1) \times 2}{(d - 1) x 2}, where the rows
#' contain the parameter vectors \eqn{\alpha} of size 2 with positve entries
#' for each of the edges in \code{tree}, if \code{model = dirichlet}.
#' }
#'
#' @return
#' Numeric matrix of size \eqn{n \times d}{n x d} of simulations of the
#' multivariate max-stable distribution.
#'
#' @details
#' The simulation follows a combination of the extremal function algorithm in \insertCite{dom2016;textual}{graphicalExtremes}
#' and the theory in \insertCite{eng2019;textual}{graphicalExtremes} to sample from a single extremal function.
#' For details on the parameters of the Huesler--Reiss, logistic
#' and negative logistic distributions see \insertCite{dom2016;textual}{graphicalExtremes}, and for the Dirichlet
#' distribution see \insertCite{coles1991modelling;textual}{graphicalExtremes}.
#'
#' @examples
#' ## A 4-dimensional HR tree model
#'
#' my_tree <- igraph::graph_from_adjacency_matrix(rbind(
#'   c(0, 1, 0, 0),
#'   c(1, 0, 1, 1),
#'   c(0, 1, 0, 0),
#'   c(0, 1, 0, 0)
#' ),
#' mode = "undirected"
#' )
#' n <- 10
#' Gamma_vec <- c(.5, 1.4, .8)
#' rmstable_tree(n, "HR", tree = my_tree, par = Gamma_vec)
#'
#' ## A 4-dimensional Dirichlet model with asymmetric edge distributions
#'
#' alpha <- cbind(c(.2, 1, .5), c(1.5, .6, .8))
#' rmstable_tree(n, model = "dirichlet", tree = my_tree, par = alpha)
#' @references
#'  \insertAllCited{}
#'
#' @export
rmstable_tree <- function(n, model = c("HR", "logistic", "dirichlet")[1],
                          tree, par) {

  # methods
  model_nms <- c("HR", "logistic", "dirichlet")

  # graph theory objects ####
  # check if it is directed
  if (igraph::is_directed(tree)) {
    warning("The given tree is directed. Converted to undirected.")
    tree <- igraph::as.undirected(tree)
  }

  # set graph theory objects
  adj <- as.matrix(igraph::as_adj(tree))
  d <- NROW(adj)
  e <- igraph::ecount(tree)
  ends.mat <- igraph::ends(tree, igraph::E(tree))

  # check if it is tree
  is_connected <- igraph::is_connected(tree)
  is_tree <- is_connected & (e == d - 1)

  if (!is_tree) {
    stop("The given graph is not a tree.")
  }


  # check arguments ####
  if (d != round(d) | d < 1) {
    stop("The argument d must be a positive integer.")
  }

  if (n != round(n) | n < 1) {
    stop("The argument n must be a positive integer.")
  }

  if (!(model %in% model_nms)) {
    stop(paste("The model must be one of ", paste(model_nms, collapse = ", "),
      ".",
      sep = ""
    ))
  }

  if (model == "HR") {
    if (!is.matrix(par)) {
      if (length(par) != (d - 1)) {
        stop(paste(
          "The argument par must be a d x d matrix,",
          "or a vector with d - 1 elements, when model = HR."
        ))
      }
    } else {
      if (NROW(par) != d | NCOL(par) != d) {
        stop(paste(
          "The argument par must be a d x d matrix,",
          "or a vector with d elements, when model = HR."
        ))
      }
    }
  } else if (model == "logistic") {
    if (any(par <= 1e-12) | any(par >= 1 - 1e-12)) {
      stop(paste(
        "The elements of par must be",
        "between 1e-12 and 1 - 1e-12,",
        "when model = logistic."
      ))
    }

    if (length(par) == 1) {
      par <- rep(par, d - 1)
      warning(paste(
        "The argument par was a scalar.",
        "Converted to a vector of size d - 1 with",
        "all same entries."
      ))
    }
    if (length(par) != d - 1) {
      stop(paste(
        "The argument par must have d - 1 elements,",
        "when model = logistic."
      ))
    }
  } else if (model == "dirichlet") {
    if (NROW(par) != d - 1 | NCOL(par) != 2) {
      stop(paste(
        "The argument par must be a (d-1) x 2 ,",
        "when model = dirichlet."
      ))
    }
    if (any(par <= 1e-12)) {
      stop(paste(
        "The elements of par must be greater than 1e-12,",
        "when model = dirichlet."
      ))
    }
  }

  # prepare arguments ####
  if (model == "HR") {
    if (is.matrix(par)) {
      par.vec <- par[ends.mat]
    } else {
      par.vec <- par
    }
  } else if (model == "logistic") {
    theta <- par
  } else if (model == "dirichlet") {
    alpha.mat <- par
  }

  # function body ####
  # Define a matrix A[[k]] choosing the paths from k to other vertices
  idx.e <- matrix(0, nrow = d, ncol = d)
  idx.e[ends.mat] <- 1:e
  idx.e <- idx.e + t(idx.e)

  # e.start[[k]][h] gives the index (1 or 2) of the starting node in the h edge
  # in the tree rooted at k
  A <- e.start <- e.end <- list()
  for (k in 1:d) {
    A[[k]] <- matrix(0, nrow = d, ncol = e)
    e.start[[k]] <- e.end[[k]] <- numeric(e)
    short.paths <- igraph::shortest_paths(tree, from = k, to = 1:d)
    for (h in 1:d) {
      path <- short.paths$vpath[[h]]
      idx.tmp <- idx.e[cbind(path[-length(path)], path[-1])]
      A[[k]][h, idx.tmp] <- 1
      e.start[[k]][idx.tmp] <-
        apply(ends.mat[idx.tmp, ] ==
          matrix(path[-length(path)], nrow = length(idx.tmp), ncol = 2),
        MARGIN = 1, FUN = function(x) which(x == TRUE)
        )
      e.end[[k]][idx.tmp] <-
        apply(ends.mat[idx.tmp, ] ==
          matrix(path[-1], nrow = length(idx.tmp), ncol = 2),
        MARGIN = 1, FUN = function(x) which(x == TRUE)
        )
    }
  }


  counter <- rep(0, times = n)
  res <- matrix(0, nrow = n, ncol = d)
  for (k in 1:d) {
    poisson <- stats::rexp(n)

    while (any(1 / poisson > res[, k])) {
      ind <- (1 / poisson > res[, k])
      n.ind <- sum(ind)
      idx <- (1:n)[ind]
      counter[ind] <- counter[ind] + 1
      proc <-
        switch(model,
          "HR" =
            simu_px_tree_HR(n = n.ind, Gamma_vec = par.vec, A = A[[k]]),
          "logistic" =
            simu_px_tree_logistic(
              n = n.ind,
              theta = theta, A = A[[k]]
            ),
          "dirichlet" =
            simu_px_tree_dirichlet(
              n = n.ind,
              alpha.start =
                alpha.mat[cbind(1:e, e.start[[k]])],
              alpha.end =
                alpha.mat[cbind(1:e, e.end[[k]])],
              A = A[[k]]
            )
        )

      if (any(dim(proc) != c(n.ind, d))) {
        stop("The generated sample has wrong size.")
      }

      if (k == 1) {
        ind.upd <- rep(TRUE, times = n.ind)
      } else {
        ind.upd <- sapply(1:n.ind, function(i) {
          all(1 / poisson[idx[i]] * proc[i, 1:(k - 1)] <=
            res[idx[i], 1:(k - 1)])
        })
      }
      if (any(ind.upd)) {
        idx.upd <- idx[ind.upd]
        res[idx.upd, ] <- pmax(res[idx.upd, ], 1 / poisson[idx.upd] *
          proc[ind.upd, ])
      }
      poisson[ind] <- poisson[ind] + stats::rexp(n.ind)
    }
  }

  return(res[sample(1:NROW(res), n, replace = FALSE), ])
}
