


vCols <- as.character(seq_len(nrow(danube$info)))
eCols <- as.character(seq_len(nrow(danube$info)-1))

ind <- danube$info$AveVol > mean(danube$info$AveVol)
vCols <- as.character(danube$info$AveVol > mean(danube$info$AveVol))


ggp <- plotDanube(
    returnGGPlot = TRUE,
    # stationIndices = ind,
    xyRatio = 1,
    clipMap = 1.2,
    plotStations = TRUE,
    # useStationVolume = TRUE,
    # useConnectionVolume = TRUE,
    # vertexColors = vCols,
    labelStations = TRUE,
    # edgeColors = eCols,
    # edgeAlpha = 1,
    mapCountries = c('Germany', 'Austria', 'Czech Republic')
)
# plot(ggp)

# plotDanubePng()

plotDanube2()



# airports <- graphicalExtremes::flights$airports
# ind <- seq_along(airports$IATA)
# # Only mainland US:
# ind <- (
#   ind
#   & airports$Latitude > 24
#   & airports$Latitude < 50
#   & airports$Longitude < -60
#   & airports$Longitude > -130
# )
# # Remove NAs
# ind[is.na(ind)] <- FALSE

# plotFlights(airportIndices = ind, plotConnections = FALSE, xyRatio = 1)
