
library(tidyverse)
library(igraph)
library(graphicalExtremes)


X <- danube$data
g <- graph_from_edgelist(danube$flow_edges)
g_und <- graph_from_edgelist(danube$flow_edges, directed=F)

n <- nrow(X)
p <- 1 - floor(n^.7)/n # gives p = 0.84

eglearn_fit <- eglearn(X, p, rholist=.09, reg_method="ns", complete_Gamma=T)

data.std <- data2mpareto(X, p)
G_emp <- emp_vario(data.std)

g_fitted <- eglearn_fit$graph[[1]]
A <- igraph::as_adjacency_matrix(g_fitted, sparse=FALSE) # correct adjacency matrix

G_comp <- complete_Gamma(G_emp, g_fitted, N=100000, check_tol=0)

Theta_comp <- Gamma2Theta(G_comp)
A2 <- (abs(Theta_comp) > 1e-7) - diag(nrow(Theta_comp)) # new adj. matrix

print(max(abs(A2 - A))) # check that the adj. matrices are identical

# plot(eglearn_fit$graph[[1]], layout = danube$coords_to_plot, edge.arrow.size=.3)
# plot(Gamma2graph(G_comp), layout = danube$coords_to_plot, edge.arrow.size=.3)


