

#' @param Sigma Numeric \dxd or \d1xd1 covariance matrix.
#' @param Gamma Numeric \dxd variogram matrix.
#' @param Theta Numeric \dxd or \d1xd1 precision matrix.
#' @param k `NULL` if the input/output matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param k1 `NULL` if the input matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param k2 `NULL` if the output matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param full Logical. If `TRUE` and `!is.null(k)`,
#' the input/output matrix is a \dxd matrix with the kth row filled with zeros.
#' @param full1 Logical. If `TRUE` and `!is.null(k1)`,
#' the input is a \dxd matrix with the kth row filled with zeros.
#' @param full2 Logical. If `TRUE` and `!is.null(k2)`,
#' the output is a \dxd matrix with the kth row filled with zeros.
#' @param M Partial matrix with `NA` entries indicating missing edges.
#' @param tol Numeric scalar. Entries in the precision matrix with absolute value
#' smaller than this are considered to be zero.
#' @param check Whether to check the inputs and ensure numerical symmetry of outputs.
#' 
#' @name sharedParamsMatrixTransformations
NULL

#' Convert matrix to graph
#' 
#' Creates a graph object representing the graph structure implied by a parameter matrix.
#' 
#' @inheritParams sharedParamsMatrixTransformations
#' 
#' @return An [`igraph::graph`] object.
#' 
#' @rdname matrix2graph
#' @export
Gamma2graph <- function(Gamma, tol=get_large_tol(), check=TRUE){
    Theta <- Gamma2Theta(Gamma, check=check)
    Theta2graph(Theta, tol=tol, check=FALSE)
}
#' @rdname matrix2graph
#' @export
Sigma2graph <- function(Sigma, tol=get_large_tol(), k=NULL, full=FALSE, check=TRUE){
    Theta <- Sigma2Theta(Sigma, k1=k, full1=full, check=check)
    Theta2graph(Theta)
}
#' @rdname matrix2graph
#' @export
Theta2graph <- function(Theta, tol=get_large_tol(), k=NULL, full=FALSE, check=TRUE){
    Theta <- Theta2Theta(Theta, k1=k, full1=full, check=check)
    A <- 1*(abs(Theta) > tol)
    graph <- igraph::graph_from_adjacency_matrix(
        A,
        mode = 'undirected',
        diag = FALSE
    )
    return(graph)
}
#' @rdname matrix2graph
#' @export
partialMatrixToGraph <- function(M){
    A <- !is.na(Matrix)
    graph <- igraph::graph_from_adjacency_matrix(
        A,
        mode = 'undirected',
        diag = FALSE
    )
    return(graph)
}


#' Conversion between Hüsler-Reiss parameter matrices
#' 
#' Converts between different matrices that parametrize the same
#' Hüsler-Reiss distribution:
#' \eGamma, \eSigma, \eTheta, \eSigmaK, \eThetaK.
#' The \d1xd1 matrices \eSigmaK and \eThetaK can also be given/returned
#' as \dxd matrices with the kth row and column filled with zeros.
#' 
#' @inheritParams sharedParamsMatrixTransformations
#' 
#' @details
#' If `k`, `k1`, or `k2` is `NULL`, the corresponding `full*` argument is ignored.
#' 
#' @return
#' The desired parameter matrix corresponding to the specified inputs.
#' 
#' @rdname parameterMatrixConversion
#' @export
Theta2Theta <- function(Theta, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    if(check){
        Theta <- checkTheta(Theta, k1, full1)
    }

    # Compute full Theta matrix
    if(full1 || is.null(k1)){
        ThetaFull <- Theta
    } else{
        d <- computeD(Theta, k1, full1)
        ThetaFull <- matrix(0, d, d)
        ThetaFull[-k1,-k1] <- Theta
    }
    if(!is.null(k1)){
        # Fill missing row/column s.t. rowSums/colSums are zero
        ThetaFull[,k1] <- -rowSums(ThetaFull)
        ThetaFull[k1,] <- ThetaFull[,k1]
        ThetaFull[k1,k1] <- -sum(ThetaFull[,k1])
    }
    
    # Compute return matrix
    if(check){
        ThetaFull <- ensure_matrix_symmetry_and_truncate_zeros(ThetaFull)
    }
    if(is.null(k2)){
        return(ThetaFull)
    }
    if(!full2){
        return(ThetaFull[-k2,-k2])
    }
    ThetaFull[k2,] <- 0
    ThetaFull[,k2] <- 0
    return(ThetaFull)
}

#' @rdname parameterMatrixConversion
#' @export
Sigma2Sigma <- function(Sigma, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    if(check){
        Sigma <- checkSigma(Sigma)
    }

    d <- computeD(Sigma, k1, full1)

    # Compute full Sigma matrix
    if(full1 || is.null(k1)){
        SigmaFull <- Sigma
    } else{
        SigmaFull <- matrix(0, d, d)
        SigmaFull[-k1, -k1] <- Sigma
    }
    
    # Compute return matrix
    if(identical(k1, k2)){
        Sigma2 <- SigmaFull
    } else if(is.null(k2)){
        ID <- diag(d) - matrix(1/d, d, d)
        Sigma2 <- ID %*% SigmaFull %*% ID
    } else {
        Sigma2 <- (
            SigmaFull
            - matrix(SigmaFull[,k2], d, d, byrow = FALSE)
            - matrix(SigmaFull[,k2], d, d, byrow = TRUE)
            + SigmaFull[k2,k2]
        )

        if(!full2){
            Sigma2 <- Sigma2[-k2, -k2]
        }
    }

    # If requested, ensure numerical symmetry of return matrix
    if(check){
        Sigma2 <- ensure_matrix_symmetry_and_truncate_zeros(Sigma2)
    }
    return(Sigma2)
}

#' @rdname parameterMatrixConversion
#' @export
Gamma2Sigma <- function(Gamma, k=NULL, full=FALSE, check=TRUE){
    if(check){
        Gamma <- checkGamma(Gamma)
    }

    d <- nrow(Gamma)
    if(is.null(k)){
        ID <- diag(d) - matrix(1/d, d, d)
        Sigma <- ID %*% (-1/2 * Gamma) %*% ID
    } else{
        Sigma <- Sigma2Sigma(-1/2 * Gamma, k2=k, full2=full, check=FALSE)
    }

    if(check){
        Sigma <- ensure_matrix_symmetry_and_truncate_zeros(Sigma)
    }
    return(Sigma)
}

#' @rdname parameterMatrixConversion
#' @export
Gamma2Theta <- function(Gamma, k=NULL, full=FALSE, check=TRUE){
    Sigma <- Gamma2Sigma(Gamma, check=check)
    Theta <- Sigma2Theta(Sigma, k2=k, full2=full, check=FALSE)
    
    if(check){
        Theta <- ensure_matrix_symmetry_and_truncate_zeros(Theta)
    }
    
    return(Theta)
}

#' @rdname parameterMatrixConversion
#' @export
Sigma2Gamma <- function(Sigma, k=NULL, full=FALSE, check=TRUE){
    # The below transformation works for all dxd matrices that correspond to Gamma,
    # not just the one Sigma matrix -> only transform if not full
    if(!is.null(k) && !full){
        Sigma <- Sigma2Sigma(Sigma, k1=k, full1=full, check=check)
    } else if(check){
        Sigma <- checkSigma(Sigma, k, full)
    }

    d <- nrow(Sigma)
    oneVec <- makeOneVec(d)
    Gamma <- oneVec %*% t(diag(Sigma)) + diag(Sigma) %*% t(oneVec) - 2 * Sigma

    if(check){
        Gamma <- ensure_matrix_symmetry_and_truncate_zeros(Gamma)
    }
    return(Gamma)
}

#' @rdname parameterMatrixConversion
#' @export
Theta2Gamma <- function(Theta, k=NULL, full=FALSE, check=TRUE){
    Sigma <- Theta2Sigma(Theta, k1=k, full1=full, check=check)
    Gamma <- Sigma2Gamma(Sigma, check=FALSE)
    if(check){
        Gamma <- ensure_matrix_symmetry_and_truncate_zeros(Gamma)
    }
    return(Gamma)
}

#' @rdname parameterMatrixConversion
#' @export
Sigma2Theta <- function(Sigma, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    # Convert input to Sigma (from Sigma^k) if necessary
    Sigma <- Sigma2Sigma(Sigma, k1=k1, full1=full1, check=check)
    
    Theta <- corpcor::pseudoinverse(Sigma)

    # Convert output to Theta^k if necessary
    if(!is.null(k2)){
        Theta <- Theta2Theta(Theta, k2=k2, full2=full2, check=FALSE)
    }

    if(check){
        Theta <- ensure_matrix_symmetry_and_truncate_zeros(Theta)
    }
    return(Theta)
}

#' @rdname parameterMatrixConversion
#' @export
Theta2Sigma <- function(Theta, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    # Convert input to Theta (from Theta^k) if necessary
    Theta <- Theta2Theta(Theta, k1=k1, full1=full1, check=check)

    Sigma <- corpcor::pseudoinverse(Theta)
    
    # Convert output to Sigma^k if necessary
    if(!is.null(k2)){
        Sigma <- Sigma2Sigma(Sigma, k2=k2, full2=full2, check=FALSE)
    }
    
    if(check){
        Sigma <- ensure_matrix_symmetry_and_truncate_zeros(Sigma)
    }
    return(Sigma)
}


# par2Matrix <- NULL
# par2Gamma <- NULL
# par2Sigma <- NULL
# matrix2par <- NULL

# #' @param chi 
# chi2Gamma <- NULL
# Gamma2chi <- NULL
# Gamma2chi_3D <- NULL



checkGamma <- function(Gamma){
    return(Gamma)
}
checkSigma <- function(Sigma, k, full){
    return(Sigma)
}
checkTheta <- function(Theta, k, full){
    return(Theta)
}

computeD <- function(M, k=NULL, full=FALSE){
    # M is full matrix -> return number of rows/columns
    if(is.null(k) || full){
        return(nrow(M))
    }
    # M is d-1 x d-1 -> return 1 + number of rows/columns
    return(nrow(M) + 1)
}

