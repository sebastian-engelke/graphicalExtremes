
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}

library(ggplot2)



g <- igraph::graph_from_edgelist(danube$flow_edges)
loc <- as.matrix(danube$info[,c('PlotCoordX', 'PlotCoordY')])
plot(g, layout = loc)

loc <- as.matrix(danube$info[,c('Lat', 'Long')])

aveVol <- colMeans(danube$dailyData, na.rm = TRUE)

info <- danube$info
info$aveVol <- aveVol

# plot(g, layout = loc)

map <- map_data('world', region='germany')

makeGgp <- function(){
    ggp <- ggplot()

    ggp <- ggp + geom_polygon(
        data = map,
        aes(x = long, y = lat, group = group),
        color = "grey65",
        fill = "#f9f9f9",
    )

    ggp <- ggp + geom_point(
        data = info,
        aes(
            size = aveVol,
            x = Long,
            y = Lat
        )
    )

    fe <- danube$flow_edges
    lat <- danube$info$Lat
    long <- danube$info$Long

    df <- data.frame(
        x0 = long[fe[,1]],
        x1 = long[fe[,2]],
        y0 = lat[fe[,1]],
        y1 = lat[fe[,2]],
        size = aveVol[fe[,1]]
    )
    ggp <- ggp + geom_segment(
        data = df,
        aes(
            # size = size,
            x = x0,
            xend = x1,
            y = y0,
            yend = y1
        )
        # arrow = arrow()
        # arrow = arrow(length = unit(0.01, 'npc'))
    )

    ggp
}

plot(makeGgp())
