
devtools::load_all('.')
library(igraph)
library(tictoc)

plotWithPid <- function(g, ...){
    plot(g, vertex.label = getPids(g), ...)
}

newSeed <- floor(2^20 * runif(1))
# newSeed <- 815653
cat('Seed:', newSeed, '\n')
set.seed(newSeed)


d <- 100

g <- generate_random_connected_graph(d, p = 3/(d+1))

d <- igraph::vcount(g)

A <- igraph::as_adjacency_matrix(g, sparse=FALSE)
A <- (A == 1)
B <- !A
diag(B) <- FALSE

gList <- split_graph(g)
cat('gList (', length(gList), '):\n', sep='')
print(sapply(gList, igraph::vcount))

g1 <- gList[[which.max(sapply(gList, length))]]

# plot(g)

G <- generate_random_Gamma(d)

TOL <- 1e-9
N <- 1e5

# 0:
P <- Gamma2Theta(G)
cat('errorG:', max(abs(G - G)[A]), '\n')
cat('errorP:', max(abs(P[B])), '\n')

# old:
cat('\n=== old ===\n')
tic()
G_c <- complete_Gamma_general(G, g, N=N, tol=TOL)
toc()
P_c <- Gamma2Theta(G_c)
cat('errorG:', max(abs(G - G_c)[A]), '\n')
cat('errorP:', max(abs(P_c[B])), '\n')

# new:
cat('\n=== new ===\n')
tic()
G3 <- complete_Gamma_general_mc(G, g, N=N, tol=TOL)
# G3 <- complete_Gamma_general_sc(G, g, N=N, tol=TOL)
toc()
P3 <- Gamma2Theta(G3)
cat('errorG:', max(abs(G - G3)[A]), '\n')
cat('errorP:', max(abs(P3[B])), '\n')

sepList <- make_sep_list(g, FALSE)
sepDetailsList <- make_sep_list(g)
AList <- lapply(sepDetailsList, function(dets){
    gg <- dets$graph
    A <- igraph::as_adjacency_matrix(gg, sparse=FALSE)
})
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

