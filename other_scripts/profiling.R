library(profvis)

n <- 1e3
d <- 1e2

Sigma <- diag(d)
Gamma <- Sigma2Gamma(Sigma, full = T)
profvis(rmpareto(n, "HR", d, Gamma))
profvis(rmpareto(n, "logistic", d, 0.3))
profvis(rmpareto(n, "neglogistic", d, 1.3))
profvis(rmpareto(n, "dirichlet", d, runif(d)))


n <- 50
d <- 20
gg <- igraph::make_tree(n = d, 2, mode = "undirected")
G_tree <- fullGamma(gg, runif(d - 1))
profvis(rmpareto_tree(n, tree = gg, par = G_tree))
