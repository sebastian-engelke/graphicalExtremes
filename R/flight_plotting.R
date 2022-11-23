
# library(dplyr)
# library(ggplot2)
# library(magrittr)

source('plotThinning.R')

if(is.null(globalenv()[['DO_PLOT']])){
  DO_PLOT <- TRUE
}

### ONLY WORKS WHEN SOURCED AFTER/FROM `prep.R`
airports <- airports
connections <- connections
nDays <- length(dates)

theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank(),
                  strip.background = element_rect(fill = "white")))

my_palette <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9",
                "#CC79A7")


nFlightMatrixToConnectionList <- function(nFlightsPerConnection, directed=TRUE){
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



# Plot flight connection
plotFlights <- function(
  airportIndices = NULL,
  airports_sel = NULL,
  connections_sel = NULL,
  useNFlights = c(FALSE, FALSE),
  minNFlights = 0,
  graph = NULL,
  map = 'state',
  plotAirports = TRUE,
  plotConnections = TRUE,
  returnGGPlot = FALSE,
  vertexColors = NULL,
  vertexShapes = NULL,
  xyRatio = NULL,
  clipMap = FALSE,
  useLatex = FALSE,
  edgeAlpha = 0.2
) {

  # Use all airports, conenctions if not specified
  if(is.null(airports_sel)){
    airports_sel <- flights$airports
  }
  if(is.null(connections_sel)){
    nFlightMat <- apply(flights$flightCounts, c(1,2), sum)
    connections_sel <- nFlightMatrixToConnectionList(nFlightMat)
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
  
  if(!is.null(vertexColors)){
    airports_sel[,'vertexColor'] <- vertexColors[IATAS]
    aesVertexColor <- 'vertexColor'
  } else{
    airports_sel[,'vertexColor'] <- NA
    aesVertexColor <- NULL
  }
  
  if(!is.null(vertexShapes)){
    airports_sel[,'vertexShape'] <- as.character(vertexShapes[IATAS])
    aesVertexShape <- 'vertexShape'
  } else{
    airports_sel[,'vertexShape'] <- NA
    aesVertexShape <- NULL
  }
  
  airports_sel %<>% mutate(dailyFlights = nFlights / nDays)

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

  # Prepare connections plotting:
  if(plotConnections){
    # Make selection of connections:
    if(is.null(graph)){
      connections_sel <- (
        connections_sel
        %>% filter(departureAirport %in% IATAS)
        %>% filter(arrivalAirport %in% IATAS)
        %>% mutate(dailyFlights = nFlights / nDays)
      )
    } else{
      m <- igraph::get.edgelist(graph, names=FALSE)
      connections_graph <- data.frame(matrix(IATAS[m], ncol = 2))
      airportColNames <- c('departureAirport', 'arrivalAirport')
      colnames(connections_graph) <- airportColNames
      connections_sel <- (
        left_join(connections_graph, connections_sel, by=airportColNames)
        %>% replace_na(list(nFlights = 0))
        %>% mutate(dailyFlights = nFlights / nDays)
      )
    }

    # Add coordinates to selected connections:
    connections_sel %<>% mutate(
      x0 = airports_sel[departureAirport, 'Longitude'],
      y0 =  airports_sel[departureAirport, 'Latitude'],
      x1 = airports_sel[arrivalAirport, 'Longitude'],
      y1 =  airports_sel[arrivalAirport, 'Latitude']
    )

    # Filter out connections with <minNFlights flights:
    connections_sel %<>% filter(nFlights >= minNFlights)
  }
  
  # Specify whether to size vertices/edges by nFlights:
  if(length(useNFlights) < 2){
    useNFlights <- c(useNFlights, useNFlights)
  }
  aesSizeNodes <- NULL
  aesSizeEdges <- NULL
  if(useNFlights[1]){
    aesSizeNodes <- 'dailyFlights'
  }
  if(useNFlights[2]){
    aesSizeEdges <- 'dailyFlights'
  }

  # Main plot object:
  ggp <- (
    ggplot()
    + xlab(NULL)
    + ylab(NULL)
    + scale_x_continuous(labels = function(x) formatDegrees(x, 'EW', useLatex))
    + scale_y_continuous(labels = function(x) formatDegrees(x, 'NS', useLatex))
  )
  
  # Legend title
  if(any(useNFlights)){
    ggp <- ggp + labs(size = 'Daily Flights')
  }

  # Plot US map in background:
  if(!is.null(map)){
    dmap <- map_data(map)
    ggp <- ggp + geom_polygon(
      data = dmap,
      aes(x = long, y = lat, group = group),
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
    ggp <- ggp + coord_cartesian(
      xlim = limits$xlim,
      ylim = limits$ylim
    ) + theme(
      aspect.ratio = limits$xyRatio
    )
  }

  # Plot airports:
  if(plotAirports){
    ggp <- ggp + geom_point(
      data = airports_sel,
      aes_string(
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
    ggp <- ggp + geom_segment(
      data = connections_sel,
      aes_string(
        x = 'x0',
        xend = 'x1',
        y = 'y0',
        yend = 'y1',
        size = aesSizeEdges
      ),
      alpha = edgeAlpha
    )
  }
  
  # ggp <- ggp + geom_sf() + coord_sf(
  #   label_axes = '--EN',
  #   label_graticule = 'EN'
  # )
  
  # ggp <- ggp + theme(axis.ticks.x = )
  
  # Return ggplot object if specified:
  if(returnGGPlot){
    return(ggp)
  }

  # Actual plot:
  plot(ggp)
  return(invisible(NULL))
}

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
