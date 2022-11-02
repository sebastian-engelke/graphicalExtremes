
# devtools::load_all('.')
library(igraph)
library(tictoc)

plotWithPid <- function(g, ...){
    plot(g, vertex.label = getPids(g), ...)
}

newSeed <- floor(2^20 * runif(1))
newSeed <- 494411
cat('Seed:', newSeed, '\n')
set.seed(newSeed)


d <- 200

# g <- generate_random_connected_graph(d, p = 3/(d+1))
g <- generate_random_chordal_graph(d)
G0 <- generate_random_Gamma(d)
G0 <- ensure_symmetry(G0, Inf)

B <- getNonEdgeIndices(g)
G <- G0
G[B] <- NA

# G2 <- edmcr::edmc(G, d = d)

G1 <- complete_Gamma(G)

Theta1 <- Gamma2Theta(G1)

print(is_sym_cnd(G1))

g1 <- Gamma2graph(G1)
print(graphs_equal(g1, g))
print(max(abs(getEdgeEntries(G1 - G, g))))
print(max(abs(getNonEdgeEntries(Theta1, g))))


