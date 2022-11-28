#' Upper Danube basin dataset
#'
#' A dataset containing river discharge data for tributaries of the Danube.
#' 
#' @format A named `list` with four entries
#' \describe{
#'  \item{`data`}{A numeric matrix, containing pre-processed discharge data for each gauging station}
#'  \item{`dailyData`}{A numeric matrix, containing daily (raw) discharge data for each gauging station}
#'  \item{`info`}{A data frame, containing information about each gauging station}
#'  \item{`flow_edges`}{
#'    A two-column numeric matrix. Each row contains the indices (in `info`)
#'    of a pair of gauging stations that are directly connected by a river.
#'  }
#' }
#' 
#' @details
#' To obtain the matrix `data`, the daily discharge data from the summer months of
#' 1960 to 2010, given in `dailyData`, was declustered, yielding between seven and ten observations per year.
#' Each row corresponds to one observation from this declustered time series,
#' the non-unique rownames indicate which year an observation is from.
#' Each column corresponds to one of the gauging stations,
#' with column indices in `data` corresponding to row indices in `info`.
#' See \insertCite{asadi2015}{graphicalExtremes} for details on the preprocessing and declustering.
#' 
#' `info` is a data frame containing the following information for
#' each of the gauging stations or its corresponding catchment area:
#' \describe{
#'  \item{`RivNames`}{Name of the river at the gauging station}
#'  \item{`Lat`, `Long`}{Coordinates of the gauging station}
#'  \item{`Lat_Center`, `Long_Center`}{Coordinates of the center of the catchment corresponding to the gauging station}
#'  \item{`Alt`}{Mean altitude of the catchment}
#'  \item{`Area`}{Area of the catchment corresponding to the gauging station}
#'  \item{`Slope`}{Mean slope of the catchment}
#'  \item{`PlotCoordX`, `PlotCoordY`}{
#'    X-Y-coordinates which can be used to arrange the gauging stations when plotting a flow graph.
#'  }
#' }
#' 
#' @examples # TODO
#' 
#' @seealso [`flights`]
#' 
#' @references
#' \insertAllCited{}
#'
#' @source Bavarian Environmental Agency <https://www.gkd.bayern.de>.
#'
"danube"


#' DEPRECATED: Flights delay data
#'
#' A dataset containing daily total delays (Jan 1, 2015 -- Dec 31, 2015)
#' of Southwest Airlines at different airports in California, Nevada,
#' Arizona, Utah, and Texas.
#' 
#' This dataset was called `flights` in earlier development versions of this package
#' and is is superseded by the new dataset now called [`flights`].
#' 
#' @seealso[`flights`], [`danube`]
#'
#' @source U.S. Department of Transportation's (DOT) Bureau of Transportation Statistics <https://www.bts.gov/>.
"flights_old"

#' Flights delay data
#' 
#' A dataset containing daily total delays of major airlines in the USA.
#' The raw data was obtained from the U.S.
#' [Bureau of Transportation Statistics](https://www.bts.dot.gov/),
#' and pre-processed as described in
#' \insertCite{hen2022;textual}{graphicalExtremes}.
#' 
#' @format A named `list` with three entries:
#' \describe{
#'  \item{`airports`}{A `data.frame`, containing information about US airports}
#'  \item{`delays`}{A numeric matrix, containing daily aggregated delays at the airports in the dataset}
#'  \item{`flightCounts`}{
#'    A numeric array, containing yearly flight numbers between airports in the dataset
#'  }
#' }
#' 
#' @details
#' `flightCounts` is a three-dimensional array, containing the number of flights in the dataset
#' between each pair of airports, aggregated on a yearly basis. 
#' Each entry is the total number of flights between the departure airport (row)
#' and destination airport (column) in a given year (dimension 3).
#' This array does not contain any `NA`s, even if an airport did not operate
#' at all in a given year, which is simply indicated by zeros.
#' 
#' `delays` is a three-dimensional array containing daily total positive delays,
#' in minutes, of incoming and outgoing flights respectively.
#' Each column corresponds to an airport in the dataset and each row corresponds
#' to a day. The third dimension has length two, `'arrivals'` containing delays of
#' incoming flights and `'departures'` containing delays of outgoing flights.
#' Zeros indicate that there were flights arriving/departing at that airport
#' on a given day, but none of them had delays. `NA`s indicate that there were
#' no flights arriving/departing at that airport on that day at all.
#' 
#' `airports` is a data frame containing the following information about a number of US airports.
#' Some entries are missing, which is indicated by `NA`s.
#' \describe{
#'  \item{`IATA`}{3-letter IATA code}
#'  \item{`Name`}{name of the airport}
#'  \item{`City`}{main city served by the airport}
#'  \item{`Country`}{country or territory where the airport is located (mostly `"United States"`)}
#'  \item{`ICAO`}{4-letter ICAO code}
#'  \item{`Latitude`}{latitude of the airport, in decimal degrees}
#'  \item{`Longitude`}{longitude of the airport, in decimal degrees}
#'  \item{`Altitude`}{altitude of the airport, in feet}
#'  \item{`Timezone`}{timezone of the airport, in hours offset from UTC}
#'  \item{`DST`}{Daylight savings time used at the airport. 'A'=US/Canada, 'N'=None.}
#'  \item{`Timezone2`}{name of the timezone of the airport}
#' }
#' 
#' @references \insertAllCited{}
#' 
#' @examples 
#' # Get total number of flights in the dataset:
#' totalFlightCounts <- apply(flights$flightCounts, c(1,2), sum)
#' 
#' # Get number of flights for specific years in the dataset:
#' flightCounts_08_09 <- apply(flights$flightCounts[,,c('2008', '2009')], c(1,2), sum)
#' 
#' # Get list of connections:
#' connections <- flightCountMatrixToConnectionList(flights$flightCounts)
#' connections_08 <- flightCountMatrixToConnectionList(flights$flightCounts[,,'2008'])
#' 
#' # Get total delays (arriving + departing):
#' totalDelays <- apply(flights$delays, c(1,2), sum)
#' 
#' @seealso [`danube`], [`flightCountMatrixToConnectionList`], [`plotFlights`]
#' 
#' @source
#' Raw delays data:
#' - <https://www.bts.dot.gov/browse-statistical-products-and-data/bts-publications/airline-service-quality-performance-234-time>
#' 
#' Fields/Forms used in the raw data:
#' - <https://esubmit.rita.dot.gov/ViewReports.aspx>
#' - <https://esubmit.rita.dot.gov/On-Time-Form1.aspx>
#' - <https://esubmit.rita.dot.gov/On-Time-Form3A.aspx>
#' 
#' Airports (includes license information):
#' - <https://openflights.org/data.html>
#' 
"flights"
