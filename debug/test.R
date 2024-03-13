
library(graphicalExtremes)

FILE0 <- 'plotFlights0.pdf'
FILE1 <- 'plotFlights1.pdf'

IATAs <- getFlightDelayData("IATAs", "tcCluster")

pdf(FILE0)

MIN_N_FLIGHTS <- 20

gg0 <- plotFlights(
  IATAs,
  minNFlights = MIN_N_FLIGHTS,
  useAirportNFlights = TRUE,
  useConnectionNFlights = FALSE,
  returnGGPlot = TRUE,
  clipMap = 1.3,
  xyRatio = 1
)
plot(gg0)
dev.off()

pdf(FILE1)
g <- getFlightGraph(IATAs, minNFlights = MIN_N_FLIGHTS)
gg1 <- plotFlights(
  IATAs,
  graph = g,
  useAirportNFlights = TRUE,
  useConnectionNFlights = FALSE,
  returnGGPlot = TRUE,
  clipMap = 1.3,
  xyRatio = 1
)
plot(gg1)
dev.off()

silentRdiff <- function(from, to){
  capture.output(ret <- suppressWarnings(tools::Rdiff(from, to, Log=TRUE)))
  ret$status != 0L
}


pdfsDifferent <- silentRdiff(FILE0, FILE1)
if(pdfsDifferent){
  cat('PDF plots are different (', FILE0, ' != ', FILE1, ').\n', sep='')
} else{
  cat('PDF plots are the same (', FILE0, ' == ', FILE1, ').\n', sep='')
}
