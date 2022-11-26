
mst_HR <- function(data, p = NULL, cens = FALSE) {

  # check if you need to rescale data or not
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  # set some params
  n <- nrow(data.std)
  d <- ncol(data.std)
  graph.full <- igraph::make_full_graph(d)

  # compute weight matrix
  G.emp <- emp_vario(data = data.std)
  res <- which(upper.tri(matrix(nrow = d, ncol = d)), arr.ind = TRUE)
  if (cens) {
    bivLLH <- apply(res[, 1:2], 1, function(x) {
      fmpareto_obj <- fmpareto_HR_MLE(data = data.std[, x],
                                  init = G.emp[x[1], x[2]],
                                  cens = cens)
      par.est <- fmpareto_obj$par
      llh_hr <- -(fmpareto_obj$nllik
                  - 2 * (sum(log(data.std[which(data.std[, x[1]] > 1), x[1]]))
                         + sum(log(data.std[which(data.std[, x[2]] > 1), x[2]]))))
      c(par = par.est, llh_hr = llh_hr)
    })
  }


  if (!cens) {
    bivLLH <- apply(res[, 1:2], 1, function(x) {
      par.est <- fmpareto_HR_MLE(
        data = data.std[, x], init = G.emp[x[1], x[2]],
        cens = cens
      )$par

      llh_hr <- logLH_HR(data = data.std[, x], Gamma = par2Gamma(par.est)) +
        2 * (sum(log(data.std[, x[1]])) + sum(log(data.std[, x[2]])))

      c(par = par.est, llh_hr = llh_hr)
    })
  }

  bivLLH.mat <- par2Gamma(bivLLH["llh_hr", ])

  # Estimated tree
  mst.tree <- igraph::mst(
    graph = graph.full, weights =
      -bivLLH.mat[igraph::ends(
        graph.full,
        igraph::E(graph.full)
      )],
    algorithm = "prim"
  )

  # set graphical parameters
  mst.tree <- mst.tree

  # Estimated Gamma
  est_Gamma <- par2Gamma(bivLLH["par", ])

  # return tree
  return(list(
    tree = mst.tree,
    Gamma = complete_Gamma(graph = mst.tree, Gamma = est_Gamma)
  ))
}

fmpareto_graph_HR_add_edges <- function(data, graph, p = NULL, cens = FALSE,
                                        edges_to_add = NULL) {

  # set up main variables
  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # check if it is directed
  if (igraph::is_directed(graph)) {
    warning("The given graph is directed. Converted to undirected.")
    graph <- igraph::as.undirected(graph)
  }

  # check if it is connected
  is_connected <- igraph::is_connected(graph)

  if (!is_connected) {
    stop("The given graph is not connected.")
  }

  # check if graph is decomposable
  is_decomposable <- igraph::is_chordal(graph)$chordal
  if (!is_decomposable) {
    stop("The given graph is not decomposable (i.e., chordal).")
  }

  # check if it is block graph
  cli <- igraph::max_cliques(graph)
  ncli <- length(cli)
  min_sep <- 0

  for (i in 1:ncli) {
    cli1 <- cli[[i]]
    for (j in 1:ncli) {
      if (j <= i) {
        next
      }
      cli2 <- cli[[j]]

      min_sep <- max(min_sep, length(intersect(cli1, cli2)))

      if (min_sep > 1) {
        break
      }
    }
  }

  if (min_sep > 1) {
    stop("The given graph is not a block graph.")
  }

  # check if the number of nodes in the graph matches the number
  # of variables in the data matrix
  nnodes <- igraph::vcount(graph)
  if (nnodes != NCOL(data)) {
    stop(paste(
      "The number of nodes in the graph doesn't match with the number",
      "of variables (i.e., columns) in the data matrix."
    ))
  }

  # check if you need to rescale data or not
  if (!is.null(p)) {
    data.std <- data2mpareto(data, p)
  } else {
    data.std <- data
  }

  l <- 1

  graph.cur <- list()
  graph.cur[[l]] <- graph
  Ghat <- list()
  Ghat[[l]] <- matrix(NA, nrow = nnodes, ncol = nnodes)

  # loop through all cliques
  for (i in 1:ncli) {
    # pick the curren cliques
    cli.idx <- cli[[i]]
    # how many nodes in the current cliques?
    cli.len <- length(cli.idx)
    # compute marginal pareto, on the nodes of the current clique
    data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)

    G.est <- emp_vario(data = data.cli)
    init <- Gamma2par(G.est)
    Ghat[[l]][cli.idx, cli.idx] <- fmpareto_HR_MLE(
      data = data.cli,
      init = init, cens = cens
    )$Gamma
  }

  Ghat[[l]] <- complete_Gamma(graph = graph.cur[[l]], Gamma = Ghat[[l]])

  # if you want to add some edges
  if (!is.null(edges_to_add)) {
    # check if edges_to_add is vector
    if (is.vector(edges_to_add)) edges_to_add <- t(as.matrix(edges_to_add))

    # check if any proposed edge is already in the given graph
    adj_mat <- igraph::as_adjacency_matrix(graph, sparse = FALSE) > 0

    m <- nrow(edges_to_add)
    check_new_edges <- 0
    for (k in 1:m) {
      current_edge <- edges_to_add[k, ]

      is_already_edge <- adj_mat[current_edge[1], current_edge[2]] |
        adj_mat[current_edge[2], current_edge[1]]

      if (is_already_edge) {
        break
      }
    }

    if (is_already_edge) {
      stop(paste(
        "The argument edges_to_add cannot contain edges already",
        "present in the given graph."
      ))
    }


    stop.flag <- FALSE
    AIC <- 2 * igraph::ecount(graph.cur[[l]]) - 2 * logLH_HR(
      data = data.std,
      Gamma = Ghat[[l]], cens = cens
    )
    edges_added <- c()

    while (length(edges_to_add) != 0 & stop.flag == FALSE) {
      m <- nrow(edges_to_add)
      AIC.tmp <- rep(NA, times = m)
      Ghat.tmp <- list()

      # go through proposed edges one after the other while retaining a block
      # graph
      # m number of proposed edges
      for (k in 1:m) {
        # current temporary graph
        Ghat.tmp[[k]] <- Ghat[[l]]
        # add the current proposed edge to the graph
        graph.tmp <- igraph::add_edges(
          graph = graph.cur[[l]],
          edges = edges_to_add[k, ]
        )

        # if the obtained graph is decomposable
        if (igraph::is_chordal(graph.tmp)$chordal) {
          # find list of max cliques
          cli <- igraph::max_cliques(graph.tmp)
          # find in which clique the new proposed edge is. It can be in at most
          # one clique, otherwise, the original graph were not decomposable.
          intersections <-
            sapply(cli, FUN = function(x) length(intersect(x, edges_to_add[k, ])) == 2)
          ii <- which(intersections == TRUE)


          # only in the clique itself the separator can be of size > 1
          if (sum(sapply(cli, FUN = function(x) {
            length(intersect(x, cli[[ii]])) > 1
          })) == 1) {
            cat("\nTry edge", edges_to_add[k, ])
            cli.idx <- cli[[ii]]
            cli.len <- length(cli.idx)
            data.cli <- mparetomargins(data = data.std, set_indices = cli.idx)

            G.est <- emp_vario(data = data.cli)
            init <- Gamma2par(G.est)
            Ghat.tmp[[k]][cli.idx, cli.idx] <- fmpareto_HR_MLE(
              data = data.cli,
              init = init, cens = cens
            )$Gamma
            Ghat.tmp[[k]] <- complete_Gamma(graph = graph.tmp, Gamma = Ghat.tmp[[k]])
            AIC.tmp[k] <- 2 * igraph::ecount(graph.tmp) -
              2 * logLH_HR(data = data.std, Gamma = Ghat.tmp[[k]], cens = cens)
          }
        }
      }
      if (!all(is.na(AIC.tmp))) {
        add.idx <- which(AIC.tmp == min(AIC.tmp, na.rm = TRUE))
        cat("\nAdded edge ", edges_to_add[add.idx, ])
        l <- l + 1
        graph.cur[[l]] <-
          igraph::add_edges(graph = graph.cur[[l - 1]], edges = edges_to_add[add.idx, ])
        graph.cur[[l]] <- graph.cur[[l]]
        Ghat[[l]] <- Ghat.tmp[[add.idx]]
        AIC <- c(AIC, AIC.tmp[add.idx])
        edges_added <- rbind(edges_added, t(as.matrix(edges_to_add[add.idx, ])))
        edges_to_add <- edges_to_add[-add.idx, ]
      }
      if (all(is.na(AIC.tmp))) stop.flag <- TRUE
    }
    return(list(
      graph = graph.cur,
      Gamma = Ghat, AIC = AIC, edges_added = edges_added
    ))
  }

  return(list(graph = graph, Gamma = Ghat[[1]]))
}

emp_chi_deprecated <- function(data, p) {
  d <- ncol(data)
  res <- as.matrix(expand.grid(1:d, 1:d))
  res <- res[res[, 1] > res[, 2], , drop = FALSE]
  chi <- apply(res, 1, function(x) {
    emp_chi_multdim(cbind(data[, x[1]], data[, x[2]]), p = p)
  })
  chi.mat <- matrix(NA, ncol = d, nrow = d)
  chi.mat[res] <- chi
  chi.mat[res[, 2:1, drop = FALSE]] <- chi
  diag(chi.mat) <- 1

  return(chi.mat)
}



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
#' The corresponding \eTheta matrix produced by [Gamma2Theta] has values
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



#' DEPRECATED: Fitting extremal graphical lasso
#'
#' Fits an extremal minimum spanning tree
#'
#' @param Gamma Numeric \nxd matrix.
#' It represents a variogram matrix \eGamma.
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
#'   \item{`Gamma`}{A list of numeric \dxd estimated
#'   variogram matrices \eGamma corresponding to the fitted graphs,
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

