
# Helper function to format (coordinate) degrees used as axis labels
formatDegrees <- function(decDeg, direction = 'NS', latex=TRUE){
  dir <- 1 + (decDeg < 0)
  dirStrings <- substring(direction, dir, dir)
  if(latex){
    degString <- '^{\\circ}'
    delim <- '$'
    dirStrings <- paste0('\\mathrm{', dirStrings, '}')
  } else{
    degString <- '\u00B0'
    delim <- ''
  }
  decDeg <- abs(decDeg)
  isNa <- is.na(decDeg)
  degStrings <- rep('', length(decDeg))
  degStrings[!isNa] <- paste0(
    delim,
    formatDegrees2(decDeg[!isNa], dirStrings[!isNa], degString),
    delim
  )
  return(degStrings)
}
formatDegrees2 <- function(decDeg, dirStrings, degString){
  dms <- decDegToDegMinSec(decDeg)
  x <- paste0(dms[,1], degString, dirStrings)
  if(!any(duplicated(x))){
    return(x)
  }
  # print(x)
  x <- paste0(dms[,1], degString, ' ', sprintf('%02d', dms[,2]), "'", dirStrings)
  if(!any(duplicated(x))){
    return(x)
  }
  # print(x)
  x <- paste0(dms[,1], degString, ' ', sprintf('%02d', dms[,2]), "' ", dms[,3], '"', dirStrings)
  return(x)
}

decDegToDegMinSec <- function(decDeg, asString = FALSE){
  deg <- floor(decDeg)
  decMin <- (decDeg - deg) * 60
  min <- floor(decMin)
  decSec <- (decMin - min) * 60
  if(asString){
    return(paste(deg, min, decSec))
  }
  return(cbind(deg, min, decSec))
}

#' Compute plot limits
#' 
#' Helper function to compute the axis limits of a plot
#' with given x, y data and optionally a fixed x-y-ratio and
#' correcting the latitude/longitude scale at different latitudes
#' @keywords internal
computeLimits <- function(xData, yData, xyRatio=1, convertLatLong=TRUE, stretch = 1){
  xRange <- range(xData)
  yRange <- range(yData)

  xMid <- mean(xRange)
  yMid <- mean(yRange)

  xRadius <- diff(xRange) / 2
  yRadius <- diff(yRange) / 2

  if(is.null(xyRatio)){
    # Just don't scale
    xScale <- 1
    yScale <- 1
    xyRatio0 <- NULL
  } else{
    xyRatio0 <- xyRatio
    # Use lat/long to account for spherical coords:
    if(convertLatLong){
      xyRatio <- xyRatio / cos(pi/180 * yMid)
    }

    # Scale x, y according to xyRatio
    xScale <- xyRatio / (xRadius / yRadius)
    yScale <- 1/xScale
    
    xScale <- max(1, xScale)
    yScale <- max(1, yScale)
  }

  return(list(
    xlim = xMid + xRadius * xScale * c(-1, 1) * stretch,
    ylim = yMid + yRadius * yScale * c(-1, 1) * stretch,
    xyRatio = xyRatio0
  ))
}
