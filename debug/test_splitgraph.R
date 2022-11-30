
devtools::load_all('.')
library(igraph)
library(tictoc)

cat('Go...\n')


# d <- 10
# g <- igraph::make_full_graph(d)
# g <- igraph::delete.edges(g, sample(igraph::ecount(g), 5))

g <- igraph::make_ring(7)
g <- igraph::add.edges(g, c(1, 4, 1, 5))

tic()
# maxCliques <- igraph::max_cliques(g)
gList <- split_graph(g)
toc()

x <- lapply(gList, function(g){
    sort(getPids(g))
})


plot(g)
for(gg in gList) plot(gg, vertex.label = getPids(gg))


# tic()
# cliques <- igraph::cliques(g)
# toc()
