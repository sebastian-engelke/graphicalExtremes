

dd <- getFlightDelayData(
    IATAfilter = 'tcCluster',
    dateFilter = 'tcTrain',
    delayFilter = c('totals')
)

print(dim(dd))

print(str(dd))
