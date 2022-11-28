
devtools::load_all('.')
library(igraph)
library(tictoc)

# set.seed(1)

newSeed <- floor(2^20 * runif(1))
# newSeed <- 856790
set.seed(newSeed)
cat('Seed:', newSeed, '\n')

cat('Go...\n')


d <- 100

# p <- runif(1, 4/(d+1), 1)
p <- runif(1) * 1/6 + 1/6

cat('p =', p, '\n')

# g <- generate_random_connected_graph(d, p = 4/(d+1))
g0 <- generate_random_connected_graph(d, p = p)
tic()
gList <- split_graph(g0)
toc()
g <- gList[[which.max(sapply(gList, igraph::vcount))]]
g <- setPids(g)
graph <- g

d <- igraph::vcount(g)
m <- igraph::ecount(g)
p1 <- m / (d*(d-1)/2)
cat(d, 'vertices,', m, 'edges, density:', p1, '\n')

if(d < 3 || p1 == 1){
    stop('Too small g')
}

A <- igraph::as_adjacency_matrix(g, sparse=FALSE)
B <- (A == 0)
diag(B) <- FALSE


# tic()
# gl1 <- make_graph_list(g)
# toc()
# print(length(gl1$graphs))
# sepList1 <- lapply(gl1$partitions, function(tmp) tmp$C)


tic()
sepList2 <- make_sep_list(g)
toc()
print(length(sepList2))

print(identical(sepList1, sepList2))

