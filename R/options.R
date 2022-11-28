

#' Number of cores to be used in parallel computations
#' 
#' Helper function that returns the number of cores to be used in parallel computations.
#' Will always be 1 on windows. On other systems, this value can be set using
#' `setOption('graphicalExtremes.mc.cores', ...)`.
#' 
#' @param overwrite Use this value (if it is valid and not on windows)
#' @return An integer to be used as number of cores
get_mc_cores <- function(overwrite = NULL){
    # Always 1 on windows
    if(.Platform$OS.typ == 'windows'){
        return(1L)
    }
    # Use overwrite if specified
    if(is.numeric(overwrite) && length(overwrite) >= 1){
        nc_overwrite <- fitInInterval(overwrite[1], 1, Inf)
        return(as.integer(nc_overwrite))
    }
    # Try to use package option
    nc_option <- getOption('graphicalExtremes.mc.cores')
    if(is.numeric(nc_option) && length(nc_option) >= 1){
        nc_option <- fitInInterval(nc_option[1], 1, Inf)
        return(as.integer(nc_option))
    }
    # Try to detect
    nc_detected <- parallel::detectCores()
    if(!is.na(nc_detected)){
        return(nc_detected)
    }
    # Fall back to 1
    return(1L)
}


#' Tolerance to be used instead of machine precision
#' 
#' Helper function that returns the tolerance to be used e.g. when checking 
#' matrices for symmetry and definiteness.
#' Can be set using `setOption('graphicalExtremes.tol', ...)`.
#' 
#' In general, this value is used only as a "permissive" tolerance, in the sense
#' that if a value has to be positive, it is compared to actual zero, but if
#' it has to be zero, its absolute value is compared to this tolerance.
#' 
#' @param overwrite Use this value if specified
#' @return A non-negative numerical scalar
get_tol <- function(overwrite = NULL){
    # Use overwrite if specified
    if(is.numeric(overwrite) && length(overwrite) >= 1){
        tol_overwrite <- fitInInterval(overwrite[1], 0, Inf)
        return(tol_overwrite)
    }
    # Try package option
    tol_option <- getOption('graphicalExtremes.tol')
    if(is.numeric(tol_option) && length(tol_option) >= 1){
        tol_option <- fitInInterval(tol_option[1], 0, Inf)
        return(tol_option)
    }
    # Fall back to some default value
    return(1e-12)
}

