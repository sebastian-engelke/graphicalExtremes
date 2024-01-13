


#' Plot flight data
#'
#' Plotting function to visualize the flight connections from the [`flights`] dataset.
#' This function requires the package `ggplot2` to be installed.
#'
#' @param airportIndices The indices of the airports (w.r.t. `airports_sel`) to include.
#' @param airports_sel The airports to plot. Might be further subset by arguments `airportIndices`, `graph`.
#' If `NULL`, then [`flights`]`$airports` will be used.
#' @param connections_sel A three columns data frame as output by [flightCountMatrixToConnectionList()].
#' If `NULL`, then [`flights`]`$nFlights` will be used to construct one.
#' @param graph An optional [`igraph::graph`] object, containing a flight graph to plot.
#' Vertices should either match the selected airports in number and order,
#' or be named with the corresponding IATA codes of the airports they represent.
#' @param plotAirports Logical. Whether to plot the airports specified.
#' @param plotConnections Logical. Whether to plot the connections specified.
#' @param labelAirports Logical. Whether to show the IATA code next to each plotted airport.
#' @param returnGGPlot If `TRUE`, a [`ggplot2::ggplot`] object is returned and not plotted immediately.
#' @param useAirportNFlights Logical. Whether to vary the size of the circles representing airports in the plot,
#' according to the number of flights at that airport.
#' @param useConnectionNFlights Logical. Whether to vary the size of the edges representing connections in the plot,
#' according to the number of flights on that connection.
#' @param minNFlights Numeric scalar. Only plot connections with at least this many flights.
#' @param map String or [`data.frame`] or `NULL`. What map to use as the background image.
#' Strings are passed to [`ggplot2::map_data()`], data frames are assumed to be the output of [`ggplot2::map_data()`].
#' @param vertexColors Optional vector, named with IATA codes, to be used as colors for the vertices/airports.
#' @param vertexShapes Optional vector, named with IATA codes, to be used as shapes for the vertices/airports. Is coerced to `character`.
#' @param edgeColors Optional vector or symmetric matrix (character or numeric), to be used as colors for edges/connections.
#' If this is a vector, its entries must match the plotted connections (in the order specified in `connections_sel` or implied by [`igraph::get.edgelist`]).
#' If this is a matrix, its row/column names must be IATA codes, or its rows/columns match the plotted airports (in number and order).
#' @param xyRatio Approximate X-Y-ratio (w.r.t. distance on the ground) of the area shown in the plot.
#' @param clipMap Logical or numeric scalar. Whether to ignore the map image when determining the axis limits of the plot.
#' If it is a positive scalar, the plot limits are extended by that factor.
#' @param useLatex Whether to format numbers etc. as latex code (useful when plotting to tikz).
#' @param edgeAlpha Numeric scalar between 0 and 1. The alpha value to be used when plotting edges/connections.
#'
#' @return If `returnGGPlot` is `TRUE`, a [`ggplot2::ggplot`] object, otherwise `NULL`.
#' @examples
#' # Plot all airports in the dataset
#' plotFlights(plotConnections = FALSE, map = 'world')
#' 
#' # Plot a selection of airports
#' plotFlights(c('JFK', 'SFO', 'LAX'), useConnectionNFlights = TRUE, useAirportNFlights = TRUE)
#' 
#' # Plot airports with a custom connections graph
#' IATAs <- c('ACV', 'BFL', 'EUG', 'SFO', 'MRY')
#' graph <- igraph::make_full_graph(length(IATAs))
#' plotFlights(IATAs, graph=graph, clipMap = 1.5)
#' 
#' @family flightData
#' @seealso [`plotDanube`]
#' 
#' @export
plotFlights <- function(
  airportIndices = NULL,
  airports_sel = NULL,
  connections_sel = NULL,
  graph = NULL,
  plotAirports = TRUE,
  plotConnections = TRUE,
  labelAirports = FALSE,
  returnGGPlot = FALSE,
  useAirportNFlights = FALSE,
  useConnectionNFlights = FALSE,
  minNFlights = 0,
  map = 'state',
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
  flights <- graphicalExtremes::flights

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
  if(is.null(map) || is.na(map) || identical(map, '')){
    map <- NULL
  }
  # Make selection of airports:
  if(!is.null(airportIndices)){
    airports_sel <- airports_sel[airportIndices,]
  }
  IATAS <- airports_sel[,'IATA']
  # If clipMap is numeric, it is used to zoom in/out of the clipped map:
  stretchMap <- 1
  if(is.numeric(clipMap)){
    stretchMap <- fitInInterval(1 * clipMap, 0, Inf)
    clipMap <- (clipMap > 0)
  }
  
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
      }
      connections_sel <- connections_graph
    }


    # Add coordinates to selected connections:
    connections_sel$x0 <- airports_sel[connections_sel$departureAirport, 'Longitude']
    connections_sel$y0 <- airports_sel[connections_sel$departureAirport, 'Latitude']
    connections_sel$x1 <- airports_sel[connections_sel$arrivalAirport, 'Longitude']
    connections_sel$y1 <- airports_sel[connections_sel$arrivalAirport, 'Latitude']
  }
  
  # Handle edge coloring
  aesEdgeColor <- NULL
  if(plotConnections && !is.null(edgeColors)){
    aesEdgeColor <- 'edgeColors'
    if(is.matrix(edgeColors)){
      # Handle matrix with edge colors as entries (ignore non-edge entries)
      if(is.null(dimnames(edgeColors))){
        # No dimnames -> assume square matrix, matching IATAs
        dimnames(edgeColors) <- list(IATAS, IATAS)
      }
      # Read colors from matrix and add to connections_sel
      ind <- cbind(connections_sel$departureAirport, connections_sel$arrivalAirport)
      connections_sel$edgeColors <- edgeColors[ind]
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
    + ggplot2::theme(legend.position = 'none')
  )
  
  # Plot US map in background:
  if(!is.null(map)){
    if(is.character(map)){
      dmap <- ggplot2::map_data(map)
    } else if(is.data.frame(map)){
      dmap <- map
    } else{
      stop('Argument `map` has to be a string, data.frame, or NULL.')
    }
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
  if(!is.null(xyRatio) || (clipMap && !is.null(map))){
    xData <- airports_sel$Longitude
    yData <- airports_sel$Latitude
    if(!clipMap && !is.null(map)){
      m <- dmap
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
  if(plotAirports && labelAirports){
    ggp <- ggp + ggplot2::geom_text(
      data = airports_sel,
      ggplot2::aes_string(
        x = 'Longitude',
        y = 'Latitude',
        label = 'IATA'
      ),
      hjust = 'left',
      nudge_x = 1/2
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
  graphics::plot(ggp)
  return(invisible(NULL))
}


#' Get filtered flight delays
#' 
#' Get filtered flight delay data, containing only a selection of dates and airports.
#' Currently, all possible selections correspond to the case study in \localCiteT{hen2022}.
#' 
#' @param what Whether to get the array of delays (numerical),
#' or just the vector of airport codes (`"IATAs"`, strings)
#' or dates (as strings). Specify exactly one.
#' @param airportFilter Which airports to include. Specify exactly one. See details below.
#' @param dateFilter Which dates to include. Specify exactly one. See details below.
#' @param delayFilter Which kinds of delays to include. Specify one or more.
#' Possible values are `"arrivals"`, `"departures"`, and `"totals"` (computed as sum of arrival and departure delays).
#' 
#' @details
#' The provided lists of airports and dates correspond to the ones used in
#' the case study of \localCiteT{hen2022}.
#' The argument `airportFilter="tcCluster"` corresponds to the airports in the analyzed "Texas Cluster",
#' `airportFilter="tcAll"` corresponds to all airports used in the previous clustering step,
#' `airportFilter="all"` corresponds to all airports in the dataset.
#' 
#' Similarly, `dateFilter="tcTrain"` selects the dates from the training set,
#' `dateFilter="tcTest"` the ones from the test/validation set.
#' To get the union of these sets, specify `dateFilter="tcAll"`.
#' To get all dates in the dataset (possibly more than for "tcAll"),
#' specify `dateFilter="all"`.
#' 
#' @return
#' If `what="IATAs"` or `what="dates"`, a character vector.
#' If required, it can be converted to [`Date`] objects using [`as.Date()`].
#' 
#' If `what="delays"`, a three-dimensional array or two-dimensional matrix,
#' with dimensions corresponding to dates, airports, and delay types.
#' 
#' 
#' @references \insertAllCited{}
#' 
#' @family flightData
#' 
#' @export
getFlightDelayData <- function(
  what = c('delays', 'IATAs', 'dates'),
  airportFilter = c('all', 'tcCluster', 'tcAll'),
  dateFilter = c('all', 'tcTrain', 'tcTest', 'tcAll'),
  delayFilter = c('totals', 'arrivals', 'departures')[1]
){
  # Filenames of internal data files
  TC_IATAS_ALL <- 'Texas_cluster_IATAS_all.rds'
  TC_IATAS_CLUSTER <- 'Texas_cluster_IATAS_cluster.rds'
  TC_DATES_TRAIN <- 'Texas_cluster_days_train_data.rds'
  TC_DATES_TEST <- 'Texas_cluster_days_test_data.rds'

  # Check args
  what <- match.arg(what)
  airportFilter <- match.arg(airportFilter)
  dateFilter <- match.arg(dateFilter)
  delayFilter <- match.arg(delayFilter, c('totals', 'arrivals', 'departures'), several.ok = TRUE)
  
  # Compute filtered IATA list (if necessary)
  if(what == 'dates' || airportFilter == 'all'){
    IATAs <- dimnames(graphicalExtremes::flights$delays)[[2]]
  } else if(airportFilter == 'tcCluster'){
    IATAs <- getPackageData(TC_IATAS_CLUSTER)
  } else if(airportFilter == 'tcAll'){
    IATAs <- getPackageData(TC_IATAS_ALL)
  }

  # Compute filtered date list (if necessary)
  allDates <- dimnames(graphicalExtremes::flights$delays)[[1]]
  if(what == 'IATAs' || dateFilter == 'all'){
    dates <- allDates
  } else{
    dates <- c()
    if(dateFilter == 'tcTrain' || dateFilter == 'tcAll'){
      dates <- c(dates, getPackageData(TC_DATES_TRAIN))
    }
    if(dateFilter == 'tcTest' || dateFilter == 'tcAll'){
      dates <- c(dates, getPackageData(TC_DATES_TEST))
    }
    if(length(setdiff(dates, allDates)) > 0){
      stop(
        'The installed version of the package does not contain the full `flights` dataset. ',
        'Make sure to install from GitHub to include the full dataset.'
      )
    }
  }

  # Compute dataset to return
  if(what == 'dates'){
    return(dates)
  }
  if(what == 'IATAs'){
    return(IATAs)
  }
  # else: what == 'delays'
  filteredDelays <- graphicalExtremes::flights$delays[dates, IATAs, ]
  ret <- array(
    0,
    c(length(dates), length(IATAs), length(delayFilter)),
    list(dates, IATAs, delayFilter)
  )
  for(i in seq_along(delayFilter)){
    filter <- delayFilter[i]
    if(filter == 'totals'){
      ret[,,i] <- apply(filteredDelays, c(1,2), sum)      
    } else {
      ret[,,i] <- filteredDelays[,,filter]
    }
  }
  ret <- drop(ret)
  return(ret)
}



#' Convert flight counts to connection list
#' 
#' Convert a numeric matrix containing flight counts between airports to a data
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
#' @examples
#' flightCountMatrixToConnectionList(flights$flightCounts[1:100, 1:100, 1])
#' 
#' @family flightData
#' @export
flightCountMatrixToConnectionList <- function(nFlightsPerConnection, directed=TRUE){
  # sum up over years if necessary
  if(length(dim(nFlightsPerConnection))){
    nFlightsPerConnection <- apply(nFlightsPerConnection, c(1,2), sum)
  }
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
