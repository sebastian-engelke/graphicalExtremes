

# Runs an igraph command (or any other expression) with
# igraph_options(add.params = FALSE)
# and resets this option afterwards
without_igraph_params <- function(expr){
  tmp <- igraph_options(add.params = FALSE)
  ret <- eval(expr, parent.frame())
  igraph_options(tmp)
  return(ret)
}

makeVertexNames <- function(g){
  paste0('v', seq_along(igraph::V(g)))
}

decompose_graph <- function(g){
  # todo: handle user-defined names with duplicates?
  if(is.null(igraph::V(g)$names) || any(duplicated(igraph::V(g)$names))){
    igraph::V(g)$names <- paste0('v', seq_along(igraph::V(g)))
  }

  graphs <- list(g)
  for(cli in cliques){
  }
}

vertexNamesToIds <- function(vNames, g){
  gNames <- igraph::V(g)$names
  match(vNames, gNames)
}

split_graph_once <- function(g){
  # compute cliques and sort by size:
  cliques <- igraph::cliques(g)
  ord <- order(sapply(cliques(g), length))
  cliques <- cliques[ord]

  for(cli in cliques){
    if(igraph::is.separator(g, cli)){
    }
  }
}

split_graph <- function(g, sep){
  inds <- seq_len(igraph::vcount(g))
  inds <- inds[-sep]
  g2 <- igraph::delete.vertices(g, sep)
  comps <- getConnectedComponents(g2)
  parts <- lapply(comps, c, sep)
}

getConnectedComponents <- function(g){
  tmp <- igraph::components(g)
  components <- lapply(seq_len(tmp$no), function(i){
    which(tmp$membership == i)
  })
  return(components)
}



