
devtools::load_all('.')
library(igraph)
library(tictoc)

try(dev.off(), silent=TRUE)


set.seed(1)

d0 <- 100
g0 <- generate_random_connected_graph(d0, p = 4/(d0+1))
g0 <- setPids(g0)

g0List <- split_graph(g0)
g0pList <- lapply(g0List, getPids)
g <- g0List[[1]]
d <- igraph::vcount(g)


print(d0)
print(length(g0List))
print(d)

# plot(g0)
# for(gg in g0List){
#   plot(gg, vertex.label = getPids(gg))
# }


minSeps <- igraph::min_separators(g)
# minSeps <- igraph::min_st_separators(g) # too slow!
ratios <- sapply(minSeps, function(ms){
  gsSep <- split_graph_at_sep(g, ms)
  partSizes <- sapply(gsSep, length)
  max(partSizes) / d
})

dists <- igraph::distances(g)
