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
#' TODO: Details, add licenses!
#' TODO: convert Timezone to numeric, DST: "U"->NA
#' TODO: see paper
#' 
#' @format A named `list` with four entries:
#' `airports` is a `data.frame` containing the following information about a number of US airports:
#' \describe{
#'  \item{`IATA`}{3-letter IATA code}
#'  \item{`Name`}{name of the airport}
#'  \item{`City`}{main city served by the airport}
#'  \item{`Country`}{country or territory where the airport is located (mostly "United States")}
#'  \item{`ICAO`}{4-letter ICAO code}
#'  \item{`Latitude`}{latitude of the airport, in decimal degrees}
#'  \item{`Longitude`}{longitude of the airport, in decimal degrees}
#'  \item{`Altitude`}{altitude of the airport, in feet}
#'  \item{`Timezone`}{timezone of the airport, in hours offset from UTC}
#'  \item{`DST`}{
#'   Daylight savings time used at the airport.
#'   Either 'A' (US/Canada), 'N' (None) or NA.
#'  }
#'  \item{`Timezone2`}{name of the timezone of the airport}
#' }
#' 
#' @source
#' Data (.asc files in .zip files):
#' - https://www.bts.dot.gov/browse-statistical-products-and-data/bts-publications/airline-service-quality-performance-234-time
#' 
#' Fields/Forms:
#' - https://esubmit.rita.dot.gov/ViewReports.aspx
#'   - https://esubmit.rita.dot.gov/On-Time-Form1.aspx
#'   - https://esubmit.rita.dot.gov/On-Time-Form2A.aspx
#'   - https://esubmit.rita.dot.gov/On-Time-Form2B.aspx
#'   - https://esubmit.rita.dot.gov/On-Time-Form3A.aspx
#'   - https://esubmit.rita.dot.gov/On-Time-Form3B.aspx
#' 
#' Airports:
#' - https://openflights.org/data.html
#' 
"flights"
