
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}
library(igraph)

library(ggplot2)

# Specify years
yearNames <- as.character(seq(2010, 2020))
minNFlights <- length(yearNames) * 1000

# Compute departures + arrivals per airport
flightsPerConnection <- apply(flights$flightCounts[, , yearNames], c(1, 2), sum)
flightsPerAirport <- rowSums(flightsPerConnection) + colSums(flightsPerConnection)

# Select airports (more than minNFlights and rough westcoast coordinates)
ind <- (
  flightsPerAirport >= minNFlights
  & flights$airports$Longitude < -119
  & flights$airports$Longitude > -130
  & flights$airports$Latitude > 25
  & flights$airports$Latitude < 50
)
IATAs <- flights$airports$IATA[ind]

# Plot airports + connections with at least monthly flights
minNConnections <- length(yearNames) * 12


# Compute undirected flights per connection
flightsPerConnectionUD <- flightsPerConnection + t(flightsPerConnection)
# Consider only connections between selected airports
flightsPerConnectionUD <- flightsPerConnectionUD[IATAs, IATAs]

# Make flight graph
A <- 1 * (flightsPerConnectionUD > minNConnections)




flight_graph <- graph_from_adjacency_matrix(A, diag = FALSE, mode = "undirected")


# We use only departure delays, at selected airports, within the selected years
dates <- as.Date(dimnames(flights$delays)[[1]])
indDates <- format(dates, "%Y") %in% yearNames
mat <- flights$delays[indDates, IATAs, "departures"]
# We remove all rows containing NAs
rowHasNA <- apply(is.na(mat), 1, any)
mat <- mat[!rowHasNA, ]

p <- .7


Gamma <- emp_vario(mat, p = p)
rholist = seq(1e-4, 0.10, length.out = 10)
flights_eglasso_fit <- eglearn(mat, p = p, rholist = rholist, complete_Gamma = TRUE)

