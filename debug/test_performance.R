
devtools::load_all('.')
library(igraph)
library(tictoc)


# newSeed <- floor(2^20 * runif(1))
newSeed <- 231320
cat('Seed:', newSeed, '\n')
set.seed(newSeed)


d <- 10

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

# # old:
# tic()
# G_c <- complete_Gamma_general(G, g, N=N, tol=TOL)
# toc()
# P_c <- Gamma2Theta(G_c)
# cat('error:', max(abs(P_c[B])), '\n')

# new:
tic()
# G3 <- complete_Gamma_general_mc(G, g, N=N, tol=TOL)
G3 <- complete_Gamma_general_sc(G, g, N=N, tol=TOL)
toc()
P3 <- Gamma2Theta(G3)
cat('error:', max(abs(P3[B])), '\n')

sepList <- make_sep_list(g, FALSE)
sep <- sepList[[1]]
sepIds <- sep
sepDetails <- makeSepDetails(g, sep)

# # new_mc:
# tic()
# G3 <- complete_Gamma_general_mc(G, g, N=N, tol=TOL, mc.cores = 16)
# toc()
# P3 <- Gamma2Theta(G3)
# cat('error:', max(abs(P3[B])), '\n')


g3 <- Gamma2graph(G3)

