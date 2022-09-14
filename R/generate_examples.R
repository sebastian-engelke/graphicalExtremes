#' Generate random Hüsler--Reiss Models
#'
#' Generates a random connected graph and Gamma matrix with conditional independence
#' structure corresponding to that graph.
#'
#' @param d Number of vertices in the graph
#' @param graph_type `"tree"`, `"block"`, `"decomposable"`, or `"general"`
#' @param ... Further arguments passed to functions generating the graph and Gamma matrix
#'
#' @family Example generations
#'
#' @examples
#' set.seed(1)
#' d <- 12
#'
#' generate_random_model(d, 'tree')
#' generate_random_model(d, 'block')
#' generate_random_model(d, 'decomposable')
#' generate_random_model(d, 'general')
#'
#' @export
generate_random_model <- function(d, graph_type='tree', ...){
  graph <- if(graph_type == 'tree'){
    generate_random_tree(d)
  } else if(graph_type == 'block'){
    generate_random_chordal_graph(d, ..., block_graph = TRUE)
  } else if(graph_type == 'decomposable'){
    generate_random_chordal_graph(d, ...)
  } else if(graph_type == 'general'){
    generate_random_connected_graph(d, ...)
  } else{
    stop('Invalid graph_type!')
  }

  Gamma <- generate_random_graphical_Gamma(graph, ...)

  return(list(
    graph = graph,
    Gamma = Gamma
  ))
}

#' Generate a random Gamma matrix for a given graph
#'
#' Generates a valid Gamma matrix with conditional independence structure
#' specified by a graph
#'
#' @param graph An [igraph::graph] object
#' @param ... Furhter arguments passed to [generate_random_spd_matrix()]
#' @family Example generations
generate_random_graphical_Gamma <- function(graph, ...){
  d <- igraph::vcount(graph)
  cliques <- igraph::maximal.cliques(graph)
  P <- matrix(0, d, d)
  for(cli in cliques){
    d_cli <- length(cli)
    P_cli <- generate_random_spd_matrix(d_cli, ...)
    ID <- diag(d_cli) - matrix(1/d_cli, d_cli, d_cli)
    P_cli <- ID %*% P_cli %*% ID
    P[cli, cli] <- P[cli, cli] + P_cli
  }
  Gamma <- Theta2Gamma(P)
  return(Gamma)
}


#' Generate a random Gamma matrix containing only integers
#'
#' Generates a random variogram Matrix by producing a
#' \eqn{(d-1) \times (d-1)}{(d-1)x(d-1)} matrix `B` with random
#' integer entries between `-b` and `b`, computing `S = B %*% t(B)`,
#' and passing this `S` to [Sigma2Gamma()].
#' This process is repeated with an increasing `b` until a valid Gamma matrix
#' is produced.
#'
#' @param d Number of rows/columns in the output matrix
#' @param b Initial `b` used in the algorithm described above
#' @param b_step By how much `b` is increased in each iteration
#'
#' @return A \eqn{d \times d}{d x d} variogram matrix with integer entries
#'
#' @family Example generations
#'
#' @examples
#'
#' generate_random_integer_Gamma(5, 2, 0.1)
#'
#' @export
generate_random_integer_Gamma <- function(d, b=2, b_step=1){
  d1 <- d - 1
  if(b_step <= 0){
    stop('Make sure that b_step > 0!')
  }
  repeat {
    B <- floor(b * (stats::runif(d1**2)*2 - 1))
    B <- matrix(B, d1, d1)
    S <- B %*% t(B)
    if(matrixcalc::is.positive.definite(S)){
      break
    }
    b <- b+b_step
  }
  G <- Sigma2Gamma(S, k=1)
  return(G)
}

#' Generate a random symmetric positive definite matrix
#'
#' Generates a random \eqn{d \times d}{dxd} symmetric positive definite matrix.
#' This is done by generating a random \eqn{d \times d}{dxd} matrix `B`,
#' then computing `B %*% t(B)`,
#' and then normalizing the matrix to approximately single digit entries.
#'
#' @param d Number of rows/columns
#' @param bMin Minimum value of entries in `B`
#' @param bMax Maximum value of entries in `B`
#' @param ... Ignored, only allowed for compatibility
#' @family Example generations
generate_random_spd_matrix <- function(d, bMin=-10, bMax=10, ...){
  M <- 0
  while(!matrixcalc::is.symmetric.matrix(M)) {
    B <- matrix(bMin + stats::runif(d**2) * (bMax-bMin), d, d)
    M <- B %*% t(B)
  }
  m <- max(floor(log(det(M), 10) / d), 0)
  M <- M * 10**(-m)
  if(!matrixcalc::is.positive.definite(M)){
    stop('Failed to produce an SPD matrix!')
  }
  return(M)
}

#' Generate a random chordal graph
#'
#' Generates a random chordal graph by starting with a (small) complete graph
#' and then adding new cliques until the specified size is reached.
#' The sizes of cliques and separators can be specified.
#'
#' @param d Number of vertices in the graph
#' @param cMin Minimal size of cliques (last clique might be smaller if necessary)
#' @param cMax Maximal size of cliques
#' @param sMin Minimal size of separators
#' @param sMax Maximal size of separators
#' @param block_graph Force `sMin == sMax == 1` to produce a block graph
#' @param ... Ignored, only allowed for compatibility
#'
#' @return An [igraph::graph] object
#' @family Example generations
generate_random_chordal_graph <- function(d, cMin=2, cMax=6, sMin=1, sMax=4, block_graph=FALSE, ...){
  if(block_graph){
    sMin <- 1
    sMax <- 1
  }
  if(cMax < cMin || sMax < sMin || cMin < sMin || cMax <= sMin || cMin > d){
    stop('Inconsistent parameters')
  }
  c0 <- floor(stats::runif(1, max(cMin, sMin+1), min(d, cMax)))
  g <- igraph::make_full_graph(c0, directed = FALSE)
  missingVertices <- d - igraph::vcount(g)
  while(missingVertices > 0){
    cliques <- igraph::maximal.cliques(g)
    cliqueSizes <- sapply(cliques, length)
    cM <- max(cliqueSizes)
    c1 <- floor(stats::runif(1, cMin, cMax))
    s1 <- floor(stats::runif(1, sMin, sMax))
    dVertices <- c1 - s1
    if(dVertices <= 0 || s1 >= cM){
      next
    } else if(dVertices > missingVertices){
      dVertices <- missingVertices
      c1 <- s1 + dVertices
    }
    indLargeCliques <- which(cliqueSizes > s1)
    if(length(indLargeCliques) == 1){
      ind <- indLargeCliques[1]
    } else{
      ind <- sample(indLargeCliques, 1)
    }
    cli <- cliques[[ind]]
    sep <- as.vector(sample(cli, s1))
    oldSize <- igraph::vcount(g)
    newVertices <- oldSize + (1:dVertices)
    g <- igraph::add_vertices(g, dVertices)
    A <- igraph::as_adjacency_matrix(g, sparse = FALSE)
    A[c(sep, newVertices), c(sep, newVertices)] <- 1
    diag(A) <- 0
    g <- igraph::graph_from_adjacency_matrix(A, mode = 'undirected')
    missingVertices <- d - igraph::vcount(g)
  }
  return(g)
}


#' Generate a random connected graph
#'
#' Generates a random connected graph.
#' First tries to generate an Erdoes-Renyi graph, if that fails, falls back
#' to producing a tree and adding random edges to that tree.
#'
#' @param d Number of vertices in the graph
#' @param m Number of edges in the graph (specify this or `p`)
#' @param p Probability of each edge being in the graph (specify this or `m`)
#' @param maxTries Maximum number of tries to produce a connected Eroes-Renyi graph
#' @param ... Ignored, only allowed for compatibility
#'
#' @return An [igraph::graph] object
#'
#' @family Example generations
generate_random_connected_graph <- function(d, m=NULL, p=2/(d+1), maxTries=1000, ...){
  # Try producing an Erdoesz-Renyi graph
  # Usually works for small d / large m:
  if(!is.null(m)){
    if(m < d-1){
      stop('m must be at least d-1!')
    } else if(m == d-1){
      return(generate_random_tree(d))
    } else{
      for(i in seq_len(maxTries)){
        g <- igraph::sample_gnm(d, m)
        if(igraph::is.connected(g)){
          return(g)
        }
      }
    }
  } else{
    for(i in seq_len(maxTries)){
      g <- igraph::sample_gnp(d, p)
      if(igraph::is.connected(g)){
        return(g)
      }
    }
  }

  # Fall back to making a tree and adding m-(d-1) edges:
  if(is.null(m)){
    m <- p*d
  }
  m <- fitInInterval(m, d-1, d*(d-1) / 2)
  g <- generate_random_tree(d)

  for(i in seq_len(m - (d-1))){
    edges <- igraph::as_edgelist(igraph::complementer(g))
    r <- sample(nrow(edges), 1)
    g <- igraph::add_edges(g, edges[r,])
  }

  return(g)
}

#' Generate a random tree
#'
#' Generates a random tree from a random Pruefer sequence
#'
#' @param d Number of vertices in the graph
#'
#' @return An [igraph::graph] object
#' @family Example generations
generate_random_tree <- function(d){
  pruefer <- floor(stats::runif(d-2, 1, d-1))
  pruefer_to_graph(pruefer)
}

#' Convert a Pruefer sequence to a graph
pruefer_to_graph <- function(pruefer){
  d <- length(pruefer) + 2
  adj <- matrix(0, d, d)
  vertices <- 1:d

  while(length(pruefer)>0){
    p <- pruefer[1]
    v <- min(setdiff(vertices, pruefer))
    adj[v, p] <- 1
    vertices <- setdiff(vertices, v)
    pruefer <- pruefer[-1]
  }
  adj[vertices[1], vertices[2]] <- 1
  adj <- adj + t(adj)

  g <- igraph::graph_from_adjacency_matrix(adj, mode='undirected')
  return(g)
}



