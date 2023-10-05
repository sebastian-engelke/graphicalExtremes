

#' @param Sigma Numeric \dxd or \d1xd1 covariance matrix.
#' @param Gamma Numeric \dxd variogram matrix.
#' @param Theta Numeric \dxd or \d1xd1 precision matrix.
#' @param k `NULL` or Integer between one and `d` (dimension of the corresponding matrix).
#' @param full Logical. If `TRUE` and `k!=NULL`, include the krh row, filled with zeros.
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
#' @rdname matrix2graph
#' @export
# Gamma2graph <- function(Gamma, tol=get_large_tol()){}
# Theta2graph <- NULL
# Sigma2graph <- NULL
# partialMatrixToGraph <- NULL
NULL


#' Convert between representations of \eTheta, \eThetaK / \eSigma, \eSigmaK
#' 
#' Converts between matrices \eTheta and \eThetaK,
#' possibly for different values of k,
#' and possibly filling the kth row/colum with zeros.
#' Same for \eSigma, \eSigmaK.
#' 
#' @param k1 `NULL` if the input matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param k2 `NULL` if the output matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param full1 Logical. If `TRUE` and `!is.null(k1)`, the input is a \dxd matrix with the kth row filled with zeros.
#' @param full2 Logical. If `TRUE` and `!is.null(k2)`, the output is a \dxd matrix with the kth row filled with zeros.
#' 
#' @inheritParams sharedParamsMatrixTransformations
#' 
#' @details
#' If `k1` or `k2` is `NULL`, the corresponding `full_` argument is ignored.
#' 
#' @return
#' The precision matrix \eTheta or \eThetaK, corresponding to the specified arguments.
#' 
#' @rdname Theta2Theta
#' @export
Theta2Theta <- function(Theta, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    # TODO: check input

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
    # strict symmetry is ensured by assignments above if input is symmetric (!)
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

#' @rdname Theta2Theta
#' @return
#' The covariance matrix \eSigma or \eSigmaK, corresponding to the specified arguments.
#' @export
Sigma2Sigma <- function(Sigma, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    # TODO: check input

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
        ProjK <- makeProjK(d, k2)
        Sigma2 <- t(ProjK) %*% SigmaFull %*% ProjK
        if(!full2){
            Sigma2 <- Sigma2[-k2, -k2]
        }
    }

    # If requested, ensure numerical symmetry of return matrix
    if(check){
        Sigma2 <- ensure_matrix_symmetry(Sigma2)
    }
    return(Sigma2)
}


computeD <- function(M, k=NULL, full=FALSE){
    # M is full matrix -> return number of rows/columns
    if(is.null(k) || full){
        return(nrow(M))
    }
    # M is d-1 x d-1 -> return 1 + number of rows/columns
    return(nrow(M) + 1)
}

Gamma2Sigma <- function(Gamma, k=NULL, full=FALSE, check=TRUE){
    if(check){
        Gamma <- checkGamma(Gamma)
    }
    
    # -1/2 Gamma is a valid Sigma matrix, if positive semi-definiteness is not checked
    Sigma <- Sigma2Sigma(-1/2 * Gamma, k2=k, full2=full, check=FALSE)

    if(check){
        Sigma <- ensure_matrix_symmetry(Sigma)
    }
    return(Sigma)
}
Gamma2Theta <- function(Gamma, k=NULL, full=FALSE, check=TRUE){
    if(check){
        Gamma <- checkGamma(Gamma)
    }
    
    Sigma <- Gamma2Sigma(Gamma, check=FALSE)
    Theta <- Sigma2Theta(Sigma, k2=k, full2=full, check=FALSE)
    
    if(check){
        Theta <- ensure_matrix_symmetry(Theta)
    }
    
    return(Theta)
}
Sigma2Gamma <- function(Sigma, k=NULL, full=FALSE, check=TRUE){
    if(check){
        Sigma <- checkSigma(Sigma, k, full)
    }
    
    if(!is.null(k) && !full){
        # The below transformation works for all dxd matrices that correspond to Gamma,
        # not just the one Sigma matrix -> only transform if not full
        Sigma <- Sigma2Sigma(Sigma, k1=k, full1=full, check=FALSE)
    }
    d <- nrow(Sigma)
    oneVec <- makeOneVec(d)
    Gamma <- oneVec %*% t(diag(Sigma)) + diag(Sigma) %*% t(oneVec) - 2 * Sigma

    if(check){
        Gamma <- ensure_matrix_symmetry(Gamma)
    }
    return(Gamma)
}
# Theta2Gamma <- NULL
#' @export
Sigma2Theta <- function(Sigma, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    if(check){
        Sigma <- checkSigma(Sigma, k1, full1)
    }

    # Convert input to Sigma (from Sigma^k) if necessary
    if(!is.null(k1)){
        Sigma <- Sigma2Sigma(Sigma, k1=k1, full1=full1, check=FALSE)
    }
    
    Theta <- corpcor::pseudoinverse(Sigma)

    # Convert output to Theta^k if necessary
    if(!is.null(k2)){
        Theta <- Theta2Theta(Theta, k2=k2, full2=full2, check=FALSE)
    }

    if(check){
        Theta <- ensure_matrix_symmetry(Theta)
    }
    return(Theta)
}
#' @export
Theta2Sigma <- function(Theta, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
    if(check){
        Theta <- checkTheta(Theta, k1, full1)
    }
    
    # Convert input to Theta (from Theta^k) if necessary
    if(!is.null(k1)){
        Theta <- Theta2Theta(Theta, k1=k1, full1=full1, check=FALSE)
    }

    Sigma <- corpcor::pseudoinverse(Theta)
    
    # Convert output to Sigma^k if necessary
    if(!is.null(k2)){
        Sigma <- Sigma2Sigma(Sigma, k2=k2, full2=full2, check=FALSE)
    }
    
    if(check){
        Sigma <- ensure_matrix_symmetry(Sigma)
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

