% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{danube}
\alias{danube}
\title{Upper Danube basin dataset}
\format{
A named \code{list} with four entries
\describe{
\item{\code{data_clustered}}{A numeric matrix, containing pre-processed discharge data for each gauging station}
\item{\code{data_raw}}{A numeric matrix, containing daily (raw) discharge data for each gauging station}
\item{\code{info}}{A data frame, containing information about each gauging station}
\item{\code{flow_edges}}{
A two-column numeric matrix. Each row contains the indices (in \code{info})
of a pair of gauging stations that are directly connected by a river.
}
}
}
\source{
Bavarian Environmental Agency \url{https://www.gkd.bayern.de}.
}
\usage{
danube
}
\description{
A dataset containing river discharge data for tributaries of the Danube.
}
\details{
To obtain the matrix \code{data_clustered}, the daily discharge data from the summer months of
1960 to 2010, given in \code{data_raw}, was declustered, yielding between seven and ten observations per year.
Each row corresponds to one observation from this declustered time series,
the \emph{non-unique rownames} indicate which year an observation is from.
Each column corresponds to one of the gauging stations,
with column indices in \code{data_raw}/\code{data_clustered} corresponding to row indices in \code{info}.
See \insertCite{asadi2015}{graphicalExtremes} for details on the preprocessing and declustering.

\code{info} is a data frame containing the following information for
each of the gauging stations or its corresponding catchment area:
\describe{
\item{\code{RivNames}}{Name of the river at the gauging station}
\item{\code{Lat}, \code{Long}}{Coordinates of the gauging station}
\item{\code{Lat_Center}, \code{Long_Center}}{Coordinates of the center of the catchment corresponding to the gauging station}
\item{\code{Alt}}{Mean altitude of the catchment}
\item{\code{Area}}{Area of the catchment corresponding to the gauging station}
\item{\code{Slope}}{Mean slope of the catchment}
\item{\code{PlotCoordX}, \code{PlotCoordY}}{
X-Y-coordinates which can be used to arrange the gauging stations when plotting a flow graph.
}
}
}
\examples{
dim(danube$data_clustered)
colnames(danube$info)

}
\references{
\insertAllCited{}
}
\seealso{
Other danubeData: 
\code{\link{getDanubeFlowGraph}()},
\code{\link{plotDanube}()}

Other datasets:
\code{\link{flights}}
}
\concept{danubeData}
\concept{datasets}
\keyword{datasets}
