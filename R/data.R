#' Upper Danube basin dataset
#'
#' A dataset containing river discharge data for tributaries of Danube.
#'
#' @source Bavarian Environmental Agency <http://www.gkd.bayern.de>.
#'
"danube"


#' DEPRECATED: Flights delay data
#'
#' A dataset containing daily total delays (Jan 1, 2015 -- Dec 31, 2015)
#' of Southwest Airlines at different airports in California, Nevada,
#' Arizona, Utah, and Texas.
#' 
#' This dataset was called `flights` in earlier development versions of this package
#' and is is superseded by the new dataset now called `flights`.
#'
#' @source U.S. Department of Transportation's (DOT) Bureau of Transportation Statistics <https://www.bts.gov/>.
"flights_old"

#' Flights delay data
#' 
#' A dataset containing daily total delays of major airlines in the USA.
#' TODO:
#' - add license(s)?
#' - add interface functions?
#'   - nFlights -> connection-list?
#' - add plot function
#' 
#' @format A named `list` with three entries:
#' \describe{
#'  \item{`airports`}{A `data.frame` containing information about US airports}
#'  \item{`delays`}{Daily aggregated delays at US airports}
#'  \item{`flightCounts`}{An array containing yearly number of flights between US airports}
#' }
#' 
#' @details
#' `flightCounts` is a three-dimensional array, containing the number of flights in the dataset
#' between each pair of airports, aggregated on a yearly basis. 
#' Each entry is the total number of flights between the departure airport (row)
#' and destination airport (column) in a given year (dimension 3).
#' 
#' `delays` is a three-dimensional array containing daily total positive delays,
#' in minutes, of incoming and outgoing flights respectively.
#' Each column corresponds to an airport in the dataset and each row corresponds
#' to a day. The third dimension has length two, `'arrivals'` containing delays of
#' incoming flights and `'departures'` containing delays of outgoing flights.
#' 
#' `airports` is a data frame containing the following information about a number of US airports:
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
#' @examples 
#' # Get total number of flights in the dataset:
#' totalFlightCounts <- apply(flights$flightCounts, c(1,2), sum)
#' 
#' # Get number of flights for specific years in the dataset:
#' flightCounts_08_09 <- apply(flights$flightCounts[,,c('2008', '2009')], c(1,2), sum)
#' 
#' # Get total delays (arriving + departing):
#' totalDelays <- apply(flights$delays, c(1,2), sum)
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
#' Airports:
#' - <https://openflights.org/data.html>
#' 
"flights"
