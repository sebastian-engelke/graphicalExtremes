# Data definitions ####

# g is a square symmetric binary matrix
# interp. an undirected graph (UG)

g <- rbind(c(0, 1, 1, 0),
           c(1, 0, 0, 1),
           c(1, 0, 0, 0),
           c(0, 1, 0, 0))


# d is a positive integer
# interp. the number of dimensions

d <- 10


# no.simu is a positive integer
# interp. the number of simulations

n <- 1e3

# *** no.simu -> n_sim


# variogram is a numeric matrix with non-negative entries
# interp. the variogram parametrizing the Huesler-Reiss distribution

variogram <- rbind(c(0, 1, 1, 1),
                   c(1, 0, 2, 2),
                   c(1, 2, 0, 2),
                   c(1, 2, 2, 0))



arguments <- c("n", "model", "d", "par", "tree", "idx", "trend", "chol_mat",
               "res", "counter", "theta", "alpha", "G.vec", "A_mat", "A",
               "alpha.start", "alpha.end", "nb.edges", "Gamma", "to_plot",
               "graph", "data", "p", "S", "full", "k", "gamma", "u", "pot",
               "x", "K", "cens", "init", "maxit", "method", "convergence",
               "nnlik", "hessian", "time", "q", "thr", "sel.edges",
               "AIC", "added.edges")
length(arguments)
sort(arguments)

# Function definitions ####
# 1. Simulation functions
# rmpareto:       n model d par -> list(res, counter)
# rmpareto_tree:  n model tree par -> list(res, counter)
# rmstable        n model d par -> list(res, counter)
# rmstable_tree:  n model tree par -> list(res, counter)

# simu_px_HR:             n idx d trend chol_mat -> matrix
# simu_px_logistic:       n idx d theta -> matrix
# simu_px_neglosistic:    n idx d theta -> matrix
# simu_px_dirichlet:      n idx d alpha -> matrix
# simu_px_tree_HR:        n G.vec A_mat -> matrix
# simu_px_tree_logistic:  n idx nb.edges theta A -> matrix
# simu_px_tree_dirichlet: n alpha.start alpha.end A_mat -> matrix


# 2. Transformation functions
# Gamma2Graph:  Gamma to_plot -> graph (gamma_to_graph)
# data2mpareto: data p -> matrix       (data_to_mpareto)
# Sigma2Gamma:  S full -> matrix       (sigma_to_gamma)
# Gamma2Sigma:  Gamma k full -> matrix (gamma_to_sigma)
# fullGamma:    graph Gamma -> matrix   (rename e.g. block_gamma_completion)
# par2Gamma:    par -> matrix           (par_to_gamma) (par2gamma)
# Gamma2par:    Gamma -> matrix         (gamma_to_par) (gamma2par)
# chi2Gamma:    gamma -> numeric    (take template from Theta2Gamma)





# 3. Fitting/Estimation functions
# chi.est:        data u pot -> numeric
# Gamma2Chi_HR:   Gamma -> numeric    (take template from chi3D)
# vario.est:      data k p -> matrix
# V_HR:           x par -> vector
# logdV_HR:       x K par -> vector
# logLH_HR:       data Gamma cens -> vector
# mst_HR:         data cens -> graph
# fpareto_HR:     data cens init maxit graph method -> list(convergence, par,
#                   Gamma, nllik, hessian, time)  (wait)
# estGraph_HR:    graph data q thr cens sel.edges -> list(graph, Gamma, AIC,
#                   added.edges)  (wait)

# unif:           x -> vector
# selectEdges:    graph -> matrix


# 4. Sebastian's file
# est.chi3D:  data triplets u Gtrue pot main -> numeric
# plotChi:    Chi.emp Chi.theo is.con main PDF filename -> plot


# 5. Discarded functions
# Theta2Gamma:  theta -> numeric
# Gamma2Theta:  gamma -> numeric
# chi.mpd.est:  data pot -> numeric


