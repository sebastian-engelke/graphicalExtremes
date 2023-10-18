
#' Check input graph
#'
#' Checks that the input graph is a valid graph for an extremal graphical model.
#' If necessary, converts the graph to an undirected graph.
#' Removes vertex labels if present.
#'
#' @param graph An \[`igraph::graph`\] object.
#' @param graph_type `"general"`, `"decomposable"`, `"block"`, `"tree"`. The required type of graph.
#' @param check_connected Whether to check if the graph is connected.
#' @param nVertices The number of vertices required in the graph.
#'
#' @return The given `graph`, if necessary converted to undirected.
#' If the graph is not valid an error is thrown.
#'
#' @family Input checks
#' @keywords internal
check_graph <- function(
  graph,
  graph_type = c('general', 'decomposable', 'block', 'tree'),
  check_connected = TRUE,
  nVertices = NULL
){
  graph_type <- match.arg(graph_type)

  # check that graph is actually a graph object
  if(!igraph::is.igraph(graph)){
    stop('The given object is not an igraph-graph.')
  }

  # set up main variables
  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # check graph size
  if(!is.null(nVertices) && d != nVertices){
    stop('The given graph does not have the correct number of vertices.')
  }

  # check if it is directed
  if (igraph::is_directed(graph)) {
    graph <- igraph::as.undirected(graph)
  }

  if(!is.null(igraph::vertex_attr(graph)[['name']])){
    graph <- igraph::remove.vertex.attribute(graph, 'name')
  }

  # check if it is connected
  is_connected <- igraph::is_connected(graph)
  if (!is_connected && check_connected) {
    stop("The given graph is not connected.")
  }

  if(graph_type == 'tree'){
    if(!is_tree_graph(graph)){
      stop("The given graph is not a tree.")
    }
  } else if(graph_type == 'block'){
    if(!is_block_graph(graph)){
      stop("The given graph is not a block graph.")
    }
  } else if(graph_type == 'decomposable'){
    if(!is_decomposable_graph(graph)) {
      stop("The given graph is not decomposable.")
    }
  } else if(graph_type != 'general'){
    stop("Not a valid graph_type.")
  }

  return(graph)
}

#' Check input graph and Gamma matrix
#'
#' Checks and converts the (incomplete) Gamma matrix and graph given for a
#' HR graphical model.
#'
#' @param Gamma A Gamma matrix or vector of entries corresponding to the edges
#' of `graph`
#' @param graph A graph object or `NULL` if the graph structure is specified by
#' `NA` in the Gamma matrix
#' @param graph_type Passed to [check_graph()].
#'
#' @return A list consisting of
#' \item{`Gamma`}{The Gamma matrix given as input or implied by the input}
#' \item{`graph`}{The graph given as input or implied by the input}
#' Throws an error if the input is not valid.
#'
#' @family Input checks
#' @keywords internal
check_Gamma_and_graph <- function(Gamma, graph = NULL, graph_type = 'general'){
  # make graph from Gamma if necessary
  if (is.null(graph) && is.matrix(Gamma)) {
    graph <- partialMatrixToGraph(Gamma)
  } else if (is.null(graph)) {
    stop("Supply a graph or a valid Gamma matrix")
  }

  # check graph
  graph <- check_graph(graph, graph_type)

  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # transform Gamma if needed
  if (is.vector(Gamma)) {
    if (length(Gamma) != e) {
      stop(paste(
        "The argument Gamma must be a symmetric d x d matrix,",
        "or a vector with as many entries as the number of edges",
        "in the graph."
      ))
    }
    G <- matrix(NA, d, d)
    edgeList <- igraph::as_edgelist(graph)
    diag(G) <- 0 # diagonal = 0
    G[edgeList] <- Gamma # upper tri
    G[edgeList[,c(2,1),drop=FALSE]] <- Gamma # lower tri
    Gamma <- G
  }

  # check that Gamma is d x d:
  if (!is_symmetric_matrix(Gamma)) {
    stop(paste(
      "The argument Gamma must be a symmetric d x d matrix,",
      "or a vector with as many entries as the number of edges",
      "in the graph."
    ))
  }

  # return Gamma and graph
  return(list(
    Gamma = Gamma,
    graph = graph
  ))
}


#' Parameter matrix checks
#' 
#' Checks wheter the matrix given is a valid Huesler-Reiss parameter matrix
#' in the role of \eGamma, \eTheta, or \eSigma, respectively.
#' 
#' @inheritParams sharedParamsMatrixTransformations
#' 
#' @return The input matrix, passed through [`ensure_matrix_symmetry_and_truncate_zeros`],
#' if it is ok, or throws an error, otherwise.
#' 
#' @rdname checkGamma
#' @export
checkGamma <- function(
  Gamma,
  alert=NULL,
  tol=get_small_tol(),
  returnBoolean=FALSE
){
  alert <- get_alert_function(alert)
  if(!is_symmetric_matrix(Gamma)){
    alert('Gamma must be a symmetric matrix!')
    if(returnBoolean) return(FALSE)
  }
  maxAbsDiag <- max_without_warning(abs(diag(Gamma)))
  if(maxAbsDiag > tol){
    alert('Gamma must have 0 diagonal (err=', maxAbsDiag, ')')
    if(returnBoolean) return(FALSE)
  }
  diag(Gamma) <- 0
  smallestEntry <- min_without_warning(upper.tri.val(Gamma))
  if(smallestEntry < -tol){
    alert('Gamma must have only non-negative entries (err=', smallestEntry, ')')
    if(returnBoolean) return(FALSE)
  }
  if(!is_cnd(Gamma)){
    alert('Gamma is not a valid variogram matrix!')
    if(returnBoolean) return(FALSE)
  }
  if(returnBoolean) return(TRUE)
  return(ensure_matrix_symmetry_and_truncate_zeros(Gamma, tol))
}
is_sym_cnd <- function(M, tol=get_small_tol){
  checkGamma(M, alert=FALSE, tol=tol, returnBoolean = TRUE)
}

#' @rdname checkGamma
#' @export
checkSigmaTheta <- function(
  M,
  k,
  full,
  matrixName='Sigma',
  tol=get_small_tol(),
  alert=NULL,
  returnBoolean=FALSE
){
  alert <- get_alert_function(alert)

  if(!is_symmetric_matrix(M)){
    alert(matrixName, ' must be a symmetric matrix!')
    if(returnBoolean) return(FALSE)
  }
  
  # Compute d, check 1<=k<=d, return early for empty matrix
  d <- computeD(M, k, full)
  if(!is.null(k) && (k < 0 || k > d)){
    alert('k (', k, ') must be in 1, ..., d (', d, ').')
    if(returnBoolean) return(FALSE)
  }
  if(d == 0 || (d == 1 && !full)){
    # Empty matrix is ok
    if(returnBoolean) return(TRUE)
    return(M)
  }
  
  # Handle positive semi-definite case with zero rowsums
  if(is.null(k)){
    maxAbsRowsum <- get_max_abs_rowsum(M)
    if(maxAbsRowsum > tol){
      alert(matrixName, ' must have zero row sums (err=', maxAbsRowsum, ').')
      if(returnBoolean) return(FALSE)
    }
    if(d == 1){
      # 1d Sigma/Theta is ok, if it is 0 (i.e. maxRowSum<=tol)
      M[1,1] <- 0
      if(returnBoolean) return(TRUE)
      return(M)
    }
    ev2 <- get_critical_ev_Sigma_Theta(M)
    if(ev2 <= 0){
      alert(matrixName, ' must be pos. semi-def. (eigenvalue ', ev2, ' should be >0).')
      if(returnBoolean) return(FALSE)
    }
    if(returnBoolean) return(TRUE)
    return(ensure_matrix_symmetry_and_truncate_zeros(M))
  }

  # k not NULL, full==FALSE -> only check pos. def.
  if(!full){
    ev1 <- get_smallest_ev(M)
    if(ev1 <= 0){
      alert(matrixName, ' without k-th row/column must be pos. def. (', ev1, ' should be >0).')
      if(returnBoolean) return(FALSE)
    }
    if(returnBoolean) return(TRUE)
    return(ensure_matrix_symmetry_and_truncate_zeros(M))
  }
  
  # Remaining case: k not NULL, full==TRUE -> check that kth row/column are zero, and rest is pos. def.
  maxZeroEntry <- max_without_warning(abs(M[,k]))
  if(maxZeroEntry > tol){
    alert('The k-th row/column of ', matrixName, ' must be zero (err=', maxZeroEntry, ').')
    if(returnBoolean) return(FALSE)
  }
  # Make sure k-th row/col are exactly zero
  M[k,] <- 0
  M[,k] <- 0
  
  # Zero matrix is ok
  if(d == 1){
    if(returnBoolean) return(TRUE)
    return(M)
  }
  
  # Check that rest is pos. def.
  ev1 <- get_smallest_ev(M[-k,-k])
  if(ev1 <= 0){
    alert(matrixName, ' without k-th row/column must be pos. def. (', ev1, ' should be >0).')
    if(returnBoolean) return(FALSE)
  }
  if(returnBoolean) return(TRUE)
  return(ensure_matrix_symmetry_and_truncate_zeros(M))
}
#' @rdname checkGamma
#' @export
checkTheta <- function(
  Theta,
  k=NULL,
  full=FALSE,
  tol=get_small_tol(),
  alert=NULL,
  returnBoolean=FALSE
){
  checkSigmaTheta(Theta, k, full, 'Theta', tol, alert, returnBoolean)
}
#' @rdname checkGamma
#' @export
checkSigma <- function(
  Sigma,
  k=NULL,
  full=FALSE,
  tol=get_small_tol(),
  alert=NULL,
  returnBoolean=FALSE
){
  checkSigmaTheta(Sigma, k, full, 'Sigma', tol, alert, returnBoolean)
}
#' @rdname checkGamma
#' @param M Numeric matrix, \eGamma, \eSigma, or \eTheta.
#' @export
checkMatrix <- function(
  M,
  name=c('Gamma', 'Sigma', 'Theta')[1],
  k=NULL,
  full=FALSE,
  tol=get_small_tol(),
  alert=NULL,
  returnBoolean=FALSE
){
  if(name == 'Gamma'){
    checkGamma(M, alert=alert, returnBoolean=returnBoolean)
  } else if(name == 'Sigma'){
    checkSigma(M, k, full, alert=alert, returnBoolean=returnBoolean)
  } else if(name == 'Theta'){
    checkTheta(M, k, full, alert=alert, returnBoolean=returnBoolean)
  } else{
    stop('Invalid matrix name!')
  }
}

computeD <- function(M, k=NULL, full=FALSE){
  # M is full matrix -> return number of rows/columns
  if(is.null(k) || full){
    return(nrow(M))
  }
  # M is d-1 x d-1 -> return 1 + number of rows/columns
  return(nrow(M) + 1)
}



#' Ensure numerical matrix symmetry/zero values
#' 
#' Ensures the symmetry of a square matrix by averaging it with its transpose.
#' Optionally verifies that the matrix was close to symmetric before.
#' 
#' @param M Numeric square matrix.
#' @param checkTol Positive scalar. If the maximum absolute difference between `M`
#' and `t(M)` is larger, show a warning.
#' 
#' @return The adjusted value of `M`.
#' 
#' @rdname ensure_matrix_symmetry_and_truncate_zeros
#' @export
ensure_matrix_symmetry <- function(M, checkTol=Inf){
  if(checkTol < Inf && max_without_warning(abs(M - t(M))) > checkTol){
    warning('Matrix not symmetric (up to tolerance: ', checkTol, ')!')
  }
  (M + t(M))/2
}
#' @description Makes sure zeros are "numerically zero", by truncating all small values.
#' @param tol All entries with absolute value below this value are truncated to zero.
#' 
#' @rdname ensure_matrix_symmetry_and_truncate_zeros
#' @export
truncate_zeros <- function(M, tol=get_small_tol()){
  M[abs(M) <= tol] <- 0
  M
}
#' @rdname ensure_matrix_symmetry_and_truncate_zeros
#' @export
ensure_matrix_symmetry_and_truncate_zeros <- function(M, tol=get_small_tol(), checkTol=Inf){
  M <- ensure_matrix_symmetry(M, checkTol)
  truncate_zeros(M, tol)
}


is_symmetric_matrix <- function(M, tol=get_small_tol()){
  (
    is.matrix(M) && ncol(M) == nrow(M)
    && all(is.na(M) == t(is.na(M)))
    && max_without_warning(M - t(M), na.rm = TRUE) < tol
  )
}

is_cnd <- function(M){
  Sk <- Gamma2Sigma(M, k=1, full=FALSE, check=FALSE)
  eig <- eigen(Sk, symmetric = TRUE, only.values = TRUE)$values
  return(eig[length(eig)] > 0)
}

check_pos_def <- function(M, alert=get_alert_function()){
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  smallestEV <- eig[length(eig)]
  if(smallestEV <= 0){
    alert('Matrix is not positive definite (', smallestEV, ' should be positive)')
  }
  return(M)
}

get_max_abs_rowsum <- function(M){
  return(max_without_warning(abs(rowSums(M))))
}

is_valid_Theta <- function(M, tol=get_small_tol()){
  # Must be symmetric
  if(!is_symmetric_matrix(M, tol)){
    return(FALSE)
  }
  d <- nrow(M)
  # Empty matrix is valid
  if(d == 0){
    return(TRUE)
  }
  # Check that rowsums are zero
  if(any(abs(rowSums(M)) > tol)){
    return(FALSE)
  }
  # 1x1 matrix is just 0 => ok
  if(d == 1){
    return(TRUE)
  }
  # Check that there are d-1 positive and one zero eigenvalue
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if(abs(eig[d]) > tol || abs(eig[d-1]) <= 0){
    return(FALSE)
  }
  return(TRUE)
}

min_without_warning <- function(..., na.rm = FALSE){
  if(length(c(...)) == 0){
    return(Inf)
  }
  min(..., na.rm = na.rm)
}

max_without_warning <- function(..., na.rm = FALSE){
  if(length(c(...)) == 0){
    return(-Inf)
  }
  max(..., na.rm = na.rm)
}


is_pos_def <- function(M){
  # Assumes that input is symmetric!
  # If M is symmetric, pos.def. is equivalent to all positive eigenvalues
  get_smallest_ev(M) > 0
}

get_smallest_ev <- function(M){
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  return(eig[length(eig)])
}

get_second_smallest_ev <- function(M){
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  d <- length(eig)
  return(eig[d-1])
}

get_critical_ev_Sigma_Theta <- function(M){
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  eig <- eig[-which.min(abs(eig))]
  return(min(eig))
}
