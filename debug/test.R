
d <- 5

G <- generate_random_Gamma(d)

G[1,d] <- G[d,1] <- NA

S <- Gamma2Sigma(G, k=2, full=TRUE)

# S[1,d] <- S[d,1] <- NA

# Sigma2Sigma(S, k2=2)
