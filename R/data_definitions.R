# Data definitions

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
