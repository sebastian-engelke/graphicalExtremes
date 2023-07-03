

#' @export
plotDanube <- function(
  stationIndices = NULL,
  graph = NULL,
  plotStations = TRUE,
  plotConnections = TRUE,
  labelStations = FALSE,
  returnGGPlot = FALSE,
  useStationVolume = FALSE,
  useConnectionVolume = FALSE,
  mapCountries = c('Germany'),
  vertexColors = NULL,
  vertexShapes = NULL,
  edgeColors = NULL,
  xyRatio = NULL,
  clipMap = FALSE,
  useLatex = FALSE,
  edgeAlpha = 0.2
) {
  # Make sure ggplot2 is installed
  ggplotAvailable <- requireNamespace('ggplot2')
  if(!ggplotAvailable){
    stop('ggplot2 needs to be installed')
  }
  
  # This makes sure `flights` is available, even if the package was not called with `library()`
  danube <- graphicalExtremes::danube

  ## Fill unspecified inputs
  # Use all airports, connections if not specified
  if(is.null(stationIndices)){
    stationIndices <- seq_len(nrow(danube$info))
  } else if(is.logical(stationIndices)){
    stationIndices <- which(stationIndices)
  }
  # Make selection of airports:
  stations <- danube$info[stationIndices,]

  # Set map to NULL if not specified:
  if(is.null(mapCountries) || is.na(mapCountries) || length(mapCountries) == 0 || nchar(mapCountries) == 0){
    map <- NULL
  }

  # If clipMap is numeric, it is used to zoom in/out of the clipped map:
  stretchMap <- 1
  if(is.numeric(clipMap)){
    stretchMap <- fitInInterval(1 * clipMap, 0, Inf)
    clipMap <- (clipMap > 0)
  }
  
  # Make sure number of graph vertices and selected airports match:
  if(!is.null(graph)){
    nVertices <- igraph::vcount(graph)
    nStations <- nrow(stations)
    if(nVertices != nStations){
      stop(sprintf(
        'Number of vertices (%d) and number of selected stations (%d) are different!',
        nVertices,
        nStations
      ))
    }
  }

  # Add vertex colors to data frame if specified
  aesVertexColor <- NULL
  if(!is.null(vertexColors)){
    aesVertexColor <- 'vertexColor'
    stations$vertexColor <- vertexColors
  }

  # Add vertex shapes to data frame if specified
  aesVertexShape <- NULL
  if(!is.null(vertexShapes)){
    aesVertexShape <- 'vertexShape'
    stations$vertexShape <- as.character(vertexShapes)
  }

  # Prepare connections plotting:
  if(plotConnections){
    # Select all physical flows between selected stations
    if(is.null(graph)){
      ind <- (
        danube$flow_edges[,1] %in% stationIndices
        | danube$flow_edges[,2] %in% stationIndices
      )
      edgeList <- danube$flow_edges[ind,]
    } else{
      edgeList <- igraph::as_edgelist(graph, names=FALSE)
    }

    connections_sel <- data.frame(
      x0 = danube$info$Long[edgeList[,1]],
      x1 = danube$info$Long[edgeList[,2]],
      y0 = danube$info$Lat[edgeList[,1]],
      y1 = danube$info$Lat[edgeList[,2]],
      AveVol = danube$info$AveVol[edgeList[,1]]
    )
  }
  
  # Handle edge coloring
  aesEdgeColor <- NULL
  if(plotConnections && !is.null(edgeColors)){
    aesEdgeColor <- 'edgeColors'
    if(is.matrix(edgeColors)){
      # Assume square matrix, matching number of selected stations
      # Read colors from matrix and add to connections_sel
      connections_sel$edgeColors <- edgeColors[edgeList]
    } else if(is.vector(edgeColors)){
      # Assume the vector matches the order of connections
      connections_sel$edgeColors <- edgeColors
    } else{
      stop('Argument `edgeColors` must be a vector with one entry per edge, or a matrix.')
    }
  }
  
  # Specify whether to size vertices/edges by nFlights:
  aesSizeNodes <- NULL
  aesSizeEdges <- NULL
  if(useStationVolume){
    aesSizeNodes <- 'AveVol'
  }
  if(useConnectionVolume){
    aesSizeEdges <- 'AveVol'
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
  if(!is.null(mapCountries)){
    dmap <- ggplot2::map_data('world', mapCountries)
    ggp <- ggp + ggplot2::geom_polygon(
      data = dmap,
      ggplot2::aes_string(x = 'long', y = 'lat', group = 'group'),
      color = "grey65",
      fill = "#f9f9f9",
      size = 0.2
    )
  }
  
  # Manually set axes limits (clips map, sets aspect ratio):
  # Note: might be improved using a different crs from
  # https://ggplot2.tidyverse.org/reference/ggsf.html
  if(!is.null(xyRatio) || (clipMap && !is.null(mapCountries))){
    xData <- stations$Long
    yData <- stations$Lat
    if(!clipMap && !is.null(mapCountries)){
      m <- ggplot2::map_data('world', mapCountries)
      xData <- c(xData, range(m$long))
      yData <- c(yData, range(m$lat))
    }
    limits <- computeLimits(
      xData,
      yData,
      xyRatio = xyRatio,
      stretch = stretchMap
    )
    ggp <- ggp + ggplot2::coord_cartesian(
      xlim = limits$xlim,
      ylim = limits$ylim
    )
    if(!is.null(limits$xyRatio)){
      ggp <- ggp + ggplot2::theme(
        aspect.ratio = 1/limits$xyRatio
      )
    }
  }

  # Plot airports:
  if(plotStations){
    ggp <- ggp + ggplot2::geom_point(
      data = stations,
      ggplot2::aes_string(
        x = 'Long',
        y = 'Lat',
        size = aesSizeNodes,
        col = aesVertexColor,
        shape = aesVertexShape
      ),
      na.rm = TRUE,
      alpha = 1
    )
  }
  if(plotStations && labelStations){
    justs <- c(
      'bl',
      'bl',
      'bl',
      'bl',
      'tl',
      'tl',
      'br', #7
      'br',
      'br',
      'br',
      'cr',
      'cr', #12
      'tl',
      'cr',
      'tl',
      'tl',
      'cl',
      'cl',
      'tl', #19
      'cl',
      'cl',
      'tl', #22
      'cr',
      'br',
      'bl',
      'bl',
      'bl', #27
      'cl',
      'tl',
      'cl',
      'cl' #31
    )[stationIndices]
    stations$vjust <- sapply(justs, function(x){
      s <- substr(x, 1, 1)
      c(t = 'top', b = 'bottom', c = 'center')[s]
    })
    stations$hjust <- sapply(justs, function(x){
      s <- substr(x, 2, 2)
      c(l = 'left', r = 'right', c = 'center')[s]
    })
    stations$label <- paste0('  ', as.character(stationIndices), ' ')
    ggp <- ggp + ggplot2::geom_text(
      data = stations,
      ggplot2::aes_string(
        x = 'Long',
        y = 'Lat',
        label = 'label',
        hjust = 'hjust',
        vjust = 'vjust'
      )#,
      # hjust = 'outward',
      # vjust = 'outward'
      # position = 'top-left'
      # hjust = 'center',
      # nudge_x = 0.1
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
        col = aesEdgeColor,
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



plotDanubePng <- function(){
  pngAvailable <- requireNamespace('png')
  if(!pngAvailable){
    stop('R-package "png" needs to be installed')
  }
  print(system.file('danube.png', 'graphicalExtremes'))
  img <- png::readPNG(system.file('danube.png'))
  plot(0,type='n',axes=FALSE,ann=FALSE, xlim=c(0,1), ylim=c(0,1))
  rasterImage(img, 0, 0, 1, 1)
}









plotDanube0 <- function(){
  # Make sure ggplot2 is installed
  ggplotAvailable <- requireNamespace('ggplot2')
  if(!ggplotAvailable){
    stop('ggplot2 needs to be installed')
  }
  
  # This makes sure `danube` is available, even if the package was not called with `library()`
  danube <- graphicalExtremes::danube
  
  ggp <- ggplot2::ggplot()
  
  stretchMap <- 1.2
  xyRatio <- 1
  clipMap <- FALSE
  drawMap <- TRUE

  if(drawMap){
    dmap <- ggplot2::map_data('world', 'Germany')
    ggp <- ggp + ggplot2::geom_polygon(
      data = dmap,
      ggplot2::aes_string(x = 'long', y = 'lat', group = 'group'),
      color = "grey65",
      fill = "#f9f9f9",
      size = 0.2
    )
  }

  # Manually set axes limits (clips map, sets aspect ratio):
  # Note: might be improved using a different crs from
  # https://ggplot2.tidyverse.org/reference/ggsf.html
  if(!is.null(xyRatio) || (clipMap && drawMap)){
    xData <- danube$info$Long
    yData <- danube$info$Lat
    if(!clipMap && drawMap){
      m <- ggplot2::map_data('world', 'Germany')
      xData <- c(xData, range(m$long))
      yData <- c(yData, range(m$lat))
    }
    limits <- computeLimits(
      xData,
      yData,
      xyRatio = xyRatio,
      stretch = stretchMap
    )
    ggp <- ggp + ggplot2::coord_cartesian(
      xlim = limits$xlim,
      ylim = limits$ylim
    )
    if(!is.null(limits$xyRatio)){
      ggp <- ggp + ggplot2::theme(
        aspect.ratio = 1/limits$xyRatio
      )
    }
  }
  
  ggp <- ggp + ggplot2::geom_point(
    data = danube$info,
    ggplot2::aes_string(
      x = 'Long',
      y = 'Lat'
    )
  )
  
  edges <- data.frame(
    x0 = danube$info$Long[danube$flow_edges[,1]],
    x1 = danube$info$Long[danube$flow_edges[,2]],
    y0 = danube$info$Lat[danube$flow_edges[,1]],
    y1 = danube$info$Lat[danube$flow_edges[,2]]
  )

  ggp <- ggp + ggplot2::geom_segment(
    data = edges,
    ggplot2::aes_string(
      x = 'x0',
      xend = 'x1',
      y = 'y0',
      yend = 'y1'
    )
  )
  
  return(ggp)
}



