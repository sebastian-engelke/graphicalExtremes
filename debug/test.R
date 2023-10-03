
d <- igraph::vcount(getDanubeFlowGraph())

cc <- sample(5, d-1, replace=TRUE)

# plotDanube(edgeColors = cc)

# plotDanube2(stationIndices = c(10:20, 1, 25, ))
