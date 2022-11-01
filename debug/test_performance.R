
devtools::load_all('.')
library(igraph)
library(tictoc)


# g <- igraph::graph_from_edgelist(rbind(
#     c(1,2),
#     c(2,3),
#     c(3,4),
#     c(4,1),
#     c(3,5),
#     c(4,5)
# ), directed = FALSE)

# set.seed(1)
# g <- generate_random_cactus(60, 5, 15)

d <- 100

g <- generate_random_connected_graph(d, p = 4/(d+1))

d <- igraph::vcount(g)

A <- igraph::as_adjacency_matrix(g, sparse=FALSE)
B <- (A == 0)
diag(B) <- FALSE

gList <- split_graph(g)
print(length(gList))

g1 <- gList[[which.max(sapply(gList, length))]]

# plot(g)

G <- generate_random_Gamma(d)

TOL <- 1e-9
N <- 1e5

# old:
tic()
G_c <- complete_Gamma_general(G, g, N=N, tol=TOL)
toc()
P_c <- Gamma2Theta(G_c)
cat('error:', max(abs(P_c[B])), '\n')

# new:
tic()
G3 <- complete_Gamma_general_mc(G, g, N=N, tol=TOL)
toc()
P3 <- Gamma2Theta(G3)
cat('error:', max(abs(P3[B])), '\n')

# new_mc:
tic()
G3 <- complete_Gamma_general_mc(G, g, N=N, tol=TOL, mc.cores = 16)
toc()
P3 <- Gamma2Theta(G3)
cat('error:', max(abs(P3[B])), '\n')


g3 <- Gamma2graph(G3)

# d <- 20

# m <- 1.5*d

# g <- generate_random_chordal_graph(d, m=m, sMax = 1)
# # g <- igraph::sample_gnm(d, m=m)
# g <- setPids(g)

# subGraphs <- split_graph(g)

# print('')
# print(igraph::count_components(g))
# print(length(subGraphs))

# for(ggg in subGraphs) {
#     plot(ggg, vertex.label = getPids(ggg))
# }

# plot(g)
