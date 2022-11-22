
## Try to figure out what the coordinates in danube data are

if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}
# library(igraph)
# library(tictoc)
library(implicitExpansion)


g <- igraph::graph_from_edgelist(danube$flow_edges)

loc <- danube$coords_to_plot

loc <- as.matrix(danube$info[,c('Long', 'Lat')])

loc1 <- as.matrix(danube$info[,c('Long', 'Lat')])
loc2 <- as.matrix(danube$info[,c('Long_Center', 'Lat_Center')])

locAll <- rbind(loc1, loc2)

xLim <- range(locAll[,1])
yLim <- range(locAll[,2])

edgeData <- lapply(seq_len(nrow(danube$flow_edges)), function(i){
    edge <- danube$flow_edges[i,]
    loc1[c(edge, NA),]
})
locEdges <- do.call(rbind, edgeData)

lineData <- lapply(seq_len(nrow(loc1)), function(i){
    rbind(
        loc1[i,],
        loc2[i,],
        NA
    )
})
locLines <- do.call(rbind, lineData)

plot(loc1, xlim = xLim, ylim = yLim)
lines(locEdges)
points(loc2, col='red')
lines(locLines, col='red', lty=2)

# dev.new()
plot(g, layout=loc1)

rn <- rownames(danube$data)
colSums(rn %m==% t(unique(rn)))
