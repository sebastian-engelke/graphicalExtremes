library(profvis)

n <- 100
d <- 50

Sigma <- diag(d)
Gamma <- Sigma2Gamma(Sigma, full = T)
profvis(rmpareto(n, "HR", d, Gamma))
profvis(rmpareto(n, "logistic", d, 0.3))
profvis(rmpareto(n, "neglogistic", d, 1.3))
profvis(rmpareto(n, "dirichlet", d, runif(d)))
