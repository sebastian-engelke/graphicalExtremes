

#' @export
plotFlights <- function(
  airportIndices = NULL,
  airports_sel = NULL,
  connections_sel = NULL,
  graph = NULL,
  plotAirports = TRUE,
  plotConnections = TRUE,
  returnGGPlot = FALSE,
  useAirportNFlights = FALSE,
  useConnectionNFlights = FALSE,
  minNFlights = 0,
  map = 'state',
  vertexColors = NULL,
  vertexShapes = NULL,
  xyRatio = NULL,
  clipMap = FALSE,
  useLatex = FALSE,
  edgeAlpha = 0.2
) {
  # Make sure ggplot2 is installed
  ggplotAvailable <- require('ggplot2')
  if(!ggplotAvailable){
    stop('ggplot2 needs to be installed')
  }

  ## Fill unspecified inputs
  # Use all airports, conenctions if not specified
  if(is.null(airports_sel)){
    airports_sel <- flights$airports
  }
  # Use airportIndices from graph, if not specified:
  if(!is.null(graph)){
    vNames <- igraph::V(graph)$name
    if(!is.null(vNames) && is.null(airportIndices)){
      airportIndices <- vNames
    }
  }
  # Set map to NULL if not specified:
  if(is.null(map) || is.na(map) || nchar(map) == 0){
    map <- NULL
  }
  # Make selection of airports:
  if(!is.null(airportIndices)){
    airports_sel <- airports_sel[airportIndices,]
  }
  IATAS <- airports_sel[,'IATA']
  
  ## If necessary, compute connections and flight counts
  computeConnections <- is.null(connections_sel) && plotConnections
  computeNFlights <- is.null(airports_sel$nFlights) && useAirportNFlights && plotAirports
  if(computeConnections || computeNFlights){
    # used multiple times below
    nFlightMat <- apply(flights$flightCounts[IATAS,IATAS,], c(1,2), sum)
  }
  # Make connections list if not given
  if(computeConnections){
    connections_sel <- flightCountMatrixToConnectionList(nFlightMat)
  }
  # Compute flight counts per airport
  if(computeNFlights){
    # add arriving and departing flights
    airports_sel$nFlights <- rowSums(nFlightMat) + colSums(nFlightMat)
  }

  # Make sure number of graph vertices and selected airports match:
  if(!is.null(graph)){
    nVertices <- igraph::vcount(graph)
    nAirports <- nrow(airports_sel)
    if(nVertices != nAirports){
      stop(sprintf(
        'Number of vertices (%d) and number of selected airports (%d) are different!',
        nVertices,
        nAirports
      ))
    }
  }

  # Add vertex colors to data frame if specified
  aesVertexColor <- NULL
  if(!is.null(vertexColors)){
    aesVertexColor <- 'vertexColor'
    airports_sel$vertexColor <- vertexColors[IATAS]
  }

  # Add vertex shapes to data frame if specified
  aesVertexShape <- NULL
  if(!is.null(vertexShapes)){
    aesVertexShape <- 'vertexShape'
    airports_sel$vertexShape <- as.character(vertexShapes[IATAS])
  }

  # Prepare connections plotting:
  if(plotConnections){
    # Make selection of connections:
    if(is.null(graph)){
      # Select all connections, that are between selected airports
      # and with >= minNFlights flights:
      ind <- (
        (
          connections_sel$departureAirport %in% IATAS
          | connections_sel$arrivalAirport %in% IATAS
        )
        & connections_sel$nFlights >= minNFlights
      )
      connections_sel <- connections_sel[ind,]
    } else{
      # Convert graph to connections list
      m <- igraph::get.edgelist(graph, names=FALSE)
      connections_graph <- data.frame(matrix(IATAS[m], ncol = 2))
      airportColNames <- c('departureAirport', 'arrivalAirport')
      colnames(connections_graph) <- airportColNames

      # Read nFlights per connection from connections_sel
      rownames(connections_sel) <- paste0(
        connections_sel$departureAirport, '_', connections_sel$arrivalAirport 
      )
      rownames(connections_graph) <- paste0(
        connections_graph$departureAirport, '_', connections_graph$arrivalAirport 
      )
      if(useConnectionNFlights){
        connections_graph$nFlights <- 0
        connections_graph$nFlights <- connections_sel[
          rownames(connections_graph),
          'nFlights'
        ]
        connections_graph$nFlights[is.na(connections_graph$nFlights)] <- 1
        # for(i in seq_len(nrow(connections_graph))){
        #   ap1 <- connections_graph$departureAirport[i]
        #   ap2 <- connections_graph$arrivalAirport[i]
        #   connections_graph$nFlights[i] <- sum(connections_sel$nFlights[
        #     connections_sel$departureAirport %in% c(ap1, ap2)
        #     & connections_sel$arrivalAirport %in% c(ap1, ap2)
        #   ], na.rm = TRUE)
        # }
      }
      connections_sel <- connections_graph
    }


    # Add coordinates to selected connections:
    connections_sel$x0 <- airports_sel[connections_sel$departureAirport, 'Longitude']
    connections_sel$y0 <- airports_sel[connections_sel$departureAirport, 'Latitude']
    connections_sel$x1 <- airports_sel[connections_sel$arrivalAirport, 'Longitude']
    connections_sel$y1 <- airports_sel[connections_sel$arrivalAirport, 'Latitude']
  }
  
  # Specify whether to size vertices/edges by nFlights:
  aesSizeNodes <- NULL
  aesSizeEdges <- NULL
  if(useAirportNFlights){
    aesSizeNodes <- 'nFlights'
  }
  if(useConnectionNFlights){
    aesSizeEdges <- 'nFlights'
  }

  # Main plot object:
  ggp <- (
    ggplot2::ggplot()
    + ggplot2::xlab(NULL)
    + ggplot2::ylab(NULL)
    + ggplot2::scale_x_continuous(labels = function(x) formatDegrees(x, 'EW', useLatex))
    + ggplot2::scale_y_continuous(labels = function(x) formatDegrees(x, 'NS', useLatex))
  )
  
  # Plot US map in background:
  if(!is.null(map)){
    dmap <- ggplot2::map_data(map)
    ggp <- ggp + ggplot2::geom_polygon(
      data = dmap,
      ggplot2::aes(x = long, y = lat, group = group),
      color = "grey65",
      fill = "#f9f9f9",
      size = 0.2
    )
  }
  
  # Manually set axes limits (clips map, sets aspect ratio):
  # TODO: consider different crs from
  # https://ggplot2.tidyverse.org/reference/ggsf.html
  if(!is.null(xyRatio) || (clipMap && !is.null(map))){
    xData <- airports_sel$Longitude
    yData <- airports_sel$Latitude
    if(!clipMap && !is.null(map)){
      m <- map_data(map)
      xData <- c(xData, range(m$long))
      yData <- c(yData, range(m$lat))
    }
    limits <- computeLimits(
      xData,
      yData,
      xyRatio = xyRatio
    )
    ggp <- ggp + ggplot2::coord_cartesian(
      xlim = limits$xlim,
      ylim = limits$ylim
    ) + theme(
      aspect.ratio = limits$xyRatio
    )
  }

  # Plot airports:
  if(plotAirports){
    ggp <- ggp + ggplot2::geom_point(
      data = airports_sel,
      ggplot2::aes_string(
        x = 'Longitude',
        y = 'Latitude',
        size = aesSizeNodes,
        col = aesVertexColor,
        shape = aesVertexShape
      ),
      na.rm = TRUE,
      alpha = 1
    )
  }
  
  # Plot connections:
  if(plotConnections){
    ggp <- ggp + ggplot2::geom_segment(
      data = connections_sel,
      ggplot2::aes_string(
        x = 'x0',
        xend = 'x1',
        y = 'y0',
        yend = 'y1',
        size = aesSizeEdges
      ),
      alpha = edgeAlpha
    )
  }
  
  # Return ggplot object if specified:
  if(returnGGPlot){
    return(ggp)
  }

  # Call plot:
  plot(ggp)
  return(invisible(NULL))
}


#' Convert flight counts to connection list
#' 
#' Convert a numeric matrix containing flight counts between ariports to a data
#' frame containing a list of connections.
#' 
#' @param nFlightsPerConnection A square, numeric matrix with identical column- and row-names.
#' Each entry represents the number of flights from the airport indexing the row to
#' the airport indexing the column in some arbitrary time period.
#' @param directed Logical scalar. Whether flights A->B and B->A should be considered separately.
#' 
#' @return A data frame with columns `departureAirport`, `arrivalAirport`, `nFlights`.
#' Each row represents one connection with >=1 flights in the input matrix.
#' 
#' @export
flightCountMatrixToConnectionList <- function(nFlightsPerConnection, directed=TRUE){
  # order rows/columns:
  perm <- order(colnames(nFlightsPerConnection))
  nFlightsPerConnection <- nFlightsPerConnection[perm, perm]
  allAirports <- colnames(nFlightsPerConnection)
  
  if(directed){
    mask <- matrix(TRUE, length(allAirports), length(allAirports))
  } else{
    nFlightsPerConnection <- nFlightsPerConnection + t(nFlightsPerConnection)
    mask <- upper.tri(nFlightsPerConnection)
  }

  ind <- which(
    (nFlightsPerConnection > 0) & mask,
    arr.ind = TRUE
  )
  flightCounts <- nFlightsPerConnection[ind]
  dep <- allAirports[ind[,1]]
  arr <- allAirports[ind[,2]]
  df <- data.frame(
    departureAirport = dep,
    arrivalAirport = arr,
    nFlights = flightCounts
  )
  return(df)
}

#' Helper function to format (coordinate) degrees used as axis labels
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
  dmsString <- measurements::conv_unit(decDeg, from='dec_deg', to='deg_min_sec')
  dms <- do.call(rbind, strsplit(dmsString, ' '))
  x <- paste0(dms[,1], degString, dirStrings)
  if(!any(duplicated(x))){
    return(x)
  }
  print(x)
  x <- paste0(dms[,1], degString, ' ', dms[,2], "'", dirStrings)
  if(!any(duplicated(x))){
    return(x)
  }
  print(x)
  x <- paste0(dms[,1], degString, ' ', dms[,2], "' ", dms[,3], '"', dirStrings)
  return(x)
}

#' Helper function to compute the axis limits of a plot
#' with given x, y data and optionally a fixed x-y-ratio and
#' correcting the latitude/longitude scale at different latitudes
computeLimits <- function(xData, yData, xyRatio=1, convertLatLong=TRUE){
  xRange <- range(xData)
  yRange <- range(yData)
  
  if(is.null(xyRatio)){
    return(list(
      xlim = xRange,
      ylim = yRange
    ))
  }

  xMid <- mean(xRange)
  yMid <- mean(yRange)

  xRadius <- diff(xRange) / 2
  yRadius <- diff(yRange) / 2

  # Use lat/long to account for spherical coords:
  xyRatio0 <- xyRatio
  if(convertLatLong){
    xyRatio <- xyRatio / cos(pi/180 * yMid)
  }

  xScale <- xyRatio / (xRadius / yRadius)
  yScale <- 1/xScale
  
  xScale <- max(1, xScale)
  yScale <- max(1, yScale)

  return(list(
    xlim = xMid + xRadius * xScale * c(-1, 1),
    ylim = yMid + yRadius * yScale * c(-1, 1),
    xyRatio = xyRatio0
  ))
}
