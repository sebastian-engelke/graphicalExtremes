
devtools::load_all('.')
library(igraph)
library(tictoc)

plotWithPid <- function(g, ...){
    plot(g, vertex.label = getPids(g), ...)
}

newSeed <- floor(2^20 * runif(1))
newSeed <- 959736
cat('Seed:', newSeed, '\n')
set.seed(newSeed)


d <- 20

g <- generate_random_connected_graph(d, p = 3/(d+1))
G0 <- generate_random_Gamma(d)

B <- getNonEdgeIndices(g)
G <- G0
G[B] <- NA

# G2 <- edmcr::edmc(G, d = d)

A <- 1*!is.na(G)

compFuncs <- list(
    dpf = function(G) (edmcr::dpf(G, d-1, toler=1e-20)),
    # sdp = function(G) (edmcr::sdp(G, A))$D,
    npf = function(G) (edmcr::npf(G, A, d-1))
)


comps <- lapply(compFuncs, function(func){
    print('X')
    tic()
    tmp <- func(G)
    tmp$t <- toc()
    tmp
})

names(comps) <- names(compFuncs)

for(n in names(comps)){
    cat('\n', n, ':\n', sep='')
    G2 <- comps[[n]]$D
    print(sum(eigen(G2)$values > 0))
    print(min(eigen(Gamma2Sigma(G2, k=1))$values))
    t <- comps[[n]]$t
    print(t$toc - t$tic)
}

# print(tmp$optval)

