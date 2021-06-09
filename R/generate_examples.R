



generate_random_model <- function(d, graph_type='tree', ...){
  graph <- 0
  if(graph_type == 'tree'){
    graph <- generate_random_tree(d)
  } else if(graph_type == 'decomposable'){
    graph <- generate_random_chordal_graph(d, ...)
  } else {
    graph <- generate_random_connected_graph(d, ...)
  }
  
  Gamma <- generate_random_graphical_Gamma(graph)
  
  return(list(
    graph = graph,
    Gamma = Gamma
  ))
}

generate_random_graphical_Gamma <- function(graph){
  d <- igraph::vcount(graph)
  cliques <- igraph::maximal.cliques(graph)
  P <- matrix(0, d, d)
  for(cli in cliques){
    d_cli <- length(cli)
    P_cli <- generate_random_spd_matrix(d_cli)
    ID <- diag(d_cli) - matrix(1/d_cli, d_cli, d_cli)
    P_cli <- ID %*% P_cli %*% ID
    P[cli, cli] <- P[cli, cli] + P_cli
  }
  Gamma <- Theta2Gamma(P)
  return(Gamma)
}


generate_random_spd_matrix <- function(d, bMin=-10, bMax=10, tol=1e-6){
  B <- bMin + runif(d**2) * (bMax-bMin)
  B <- matrix(B, d, d)
  M <- B %*% t(B)
  m <- max(floor(log(det(M), 10) / d), 0)
  M <- M * 10**(-m)
  if(!matrixcalc::is.positive.definite(M)){
    stop('Failed to produce an SPD matrix!')
  }
  return(M)
}

generate_random_chordal_graph <- function(d, cMin=2, cMax=6, sMin=1, sMax=4){
  if(cMax < cMin || sMax < sMin || cMin < sMin || cMax <= sMin || cMin > d){
    stop('Inconsistent parameters')
  }
  c0 <- floor(runif(1, max(cMin, sMin+1), min(d, cMax)))
  g <- igraph::make_full_graph(c0, directed = FALSE)
  missingVertices <- d - igraph::vcount(g)
  while(missingVertices > 0){
    cliques <- igraph::maximal.cliques(g)
    cliqueSizes <- sapply(cliques, length)
    cM <- max(cliqueSizes)
    c1 <- floor(runif(1, cMin, cMax))
    s1 <- floor(runif(1, sMin, sMax))
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


generate_random_connected_graph <- function(d, m=2*d, p=NULL, maxTries=1000){
  # Try producing a Erdoesz-Renyi graph
  # Usually works for small d / large m
  if(is.null(p)){
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

  if(is.null(m)){
    m <- p*d
  }
  m <- max(m, d-1)
  m <- min(m, d * (d-1) / 2)
  g <- generate_random_tree(d)
  
  for(i in seq_len(m - (d-1))){
    edges <- igraph::as_edgelist(igraph::complementer(g))
    r <- sample(nrow(edges), 1)
    g <- igraph::add_edges(g, edges[r,])
  }

  if(!igraph::is.connected(g)){
    stop('Failed to generate a connected graph!')
  }
  return(g)
}

generate_random_tree <- function(d){
  pruefer <- floor(runif(d-2, 1, d-1))
  pruefer_to_graph(pruefer)
}

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



