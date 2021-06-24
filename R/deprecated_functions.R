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
      fmpareto_obj <- fmpareto_HR(data = data.std[, x],
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
      par.est <- fmpareto_HR(
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
    Ghat[[l]][cli.idx, cli.idx] <- fmpareto_HR(
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
            Ghat.tmp[[k]][cli.idx, cli.idx] <- fmpareto_HR(
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
