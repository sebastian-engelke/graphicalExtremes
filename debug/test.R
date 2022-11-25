
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}



# plotFlights(plotConnections = FALSE, map='world')

ind <- 20:30

g <- igraph::make_ring(length(ind))

plotFlights(ind, graph=g)
plotFlights(ind, graph=g, useConnectionNFlights=TRUE, useAirportNFlights = TRUE)
plotFlights(ind, useConnectionNFlights=TRUE, useAirportNFlights = TRUE)

R -e "? graphicalExtremes::danube; Sys.sleep(3)"
R -e "? graphicalExtremes::flights; Sys.sleep(3)"
