library(profvis)

n <- 1e3
d <- 1e2

Sigma <- diag(d)
Gamma <- Sigma2Gamma(Sigma, full = T)
profvis(rmpareto(n, "HR", d, Gamma))
profvis(rmpareto(n, "logistic", d, 0.3))
profvis(rmpareto(n, "neglogistic", d, 1.3))
profvis(rmpareto(n, "dirichlet", d, runif(d)))


n <- 1e3
d <- 1e2
gg <- igraph::make_tree(n = d, 2, mode = "undirected")
G_tree <- complete_Gamma(gg, runif(d-1))
gg_check <- Gamma2Graph(G_tree)
igraph::V(gg_check)$color <- "white"
igraph::tkplot(gg_check)
profvis(rmpareto_tree(n, tree = gg, par = G_tree))
