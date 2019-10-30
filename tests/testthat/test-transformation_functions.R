context("test-transformation_functions")


# Define variables
set.seed(1234)
non_decomposable <- igraph::graph_from_adjacency_matrix(
  rbind(c(0, 1, 0, 0, 0, 0),
        c(1, 0, 1, 0, 1, 0),
        c(0, 1, 0, 1, 0, 1),
        c(0, 0, 1, 0, 1, 1),
        c(0, 1, 0, 1, 0, 0),
        c(0, 0, 1, 1, 0, 0)),
  mode = "undirected"
)
igraph::is_chordal(non_decomposable)$chordal

non_block <- igraph::graph_from_adjacency_matrix(
  rbind(c(0, 1, 0, 0, 0, 0),
        c(1, 0, 1, 1, 1, 0),
        c(0, 1, 0, 1, 1, 1),
        c(0, 1, 1, 0, 1, 1),
        c(0, 1, 1, 1, 0, 0),
        c(0, 0, 1, 1, 0, 0)),
  mode = "undirected"
)
igraph::is_chordal(non_block)$chordal

block <- igraph::graph_from_adjacency_matrix(
  rbind(c(0, 1, 0, 0, 0, 0, 0),
        c(1, 0, 1, 1, 1, 0, 0),
        c(0, 1, 0, 1, 1, 1, 1),
        c(0, 1, 1, 0, 0, 0, 0),
        c(0, 1, 1, 1, 0, 0, 0),
        c(0, 0, 1, 0, 0, 0, 1),
        c(0, 0, 1, 0, 0, 1, 0)),
  mode = "undirected"
)
igraph::is_chordal(block)$chordal

block2 <- igraph::graph_from_adjacency_matrix(
  rbind(c(0, 1),
        c(1, 0)),
  mode = "undirected")
igraph::is_chordal(block2)$chordal

Gamma1 <- rbind(c(0, 2, 0, 0, 0, 0),
                c(2, 0, 2, 0, 2, 0),
                c(0, 2, 0, 2, 0, 2),
                c(0, 0, 2, 0, 2, 2),
                c(0, 2, 0, 2, 0, 0),
                c(0, 0, 2, 2, 0, 0))

Gamma2 <- rbind(c(0, 2, 0, 0, 0, 0),
                 c(2, 0, 2, 2, 2, 0),
                 c(0, 2, 0, 2, 2, 2),
                 c(0, 2, 2, 0, 2, 2),
                 c(0, 2, 2, 2, 0, 0),
                 c(0, 0, 2, 2, 0, 0))

Gamma3 <- rbind(c(0, 2, 0, 0, 0, 0, 0),
                c(2, 0, 2, 2, 2, 0, 0),
                c(0, 2, 0, 2, 2, 2, 2),
                c(0, 2, 2, 0, 2, 0, 0),
                c(0, 2, 2, 2, 0, 0, 0),
                c(0, 0, 2, 0, 0, 0, 2),
                c(0, 0, 2, 0, 0, 2, 0))

Gamma3_wrong <- Gamma3
Gamma3_wrong[2, 3] <- 1
Gamma3_vec <- Gamma3[upper.tri(Gamma3)]
Gamma3_vec <- Gamma3_vec[which(Gamma3_vec != 0)]
Gamma3_vec_wrong <- c(Gamma3_vec, 2)

empty_graph <- igraph::graph_from_adjacency_matrix(
  matrix(0, nrow = NROW(Gamma3), ncol = NCOL(Gamma3)),
  mode = "undirected"
)

Gamma3_completed <- rbind(c(0, 2, 4, 4, 4, 6, 6),
                          c(2, 0, 2, 2, 2, 4, 4),
                          c(4, 2, 0, 2, 2, 2, 2),
                          c(4, 2, 2, 0, 2, 4, 4),
                          c(4, 2, 2, 2, 0, 4, 4),
                          c(6, 4, 2, 4, 4, 0, 2),
                          c(6, 4, 2, 4, 4, 2, 0))

data <- mvtnorm::rmvnorm(n = 1e2, mean = runif(5))

G <-  cbind(c(0, 1.5, 1.5, 2),
            c(1.5, 0, 2, 1.5),
            c(1.5, 2, 0, 1.5),
            c(2, 1.5, 1.5, 0))

par <- G[upper.tri(G)]


# Run tests
test_that("fullGamma works", {

  expect_error(fullGamma(non_decomposable, Gamma1))
  expect_error(fullGamma(non_block, Gamma2))
  expect_warning(fullGamma(igraph::as.directed(block), Gamma3))
  expect_error(fullGamma(empty_graph, Gamma3))
  expect_error(fullGamma(block, Gamma3_wrong))
  expect_error(fullGamma(block, Gamma3_vec_wrong))
  expect_error(fullGamma(block, Gamma2))

  res1 <- fullGamma(block, Gamma3)
  res2 <- fullGamma(block, Gamma3_vec)
  expect_equal(res1, Gamma3_completed)
  expect_equal(res2, Gamma3_completed)

  res1 <- fullGamma(block2, Gamma3[1:2, 1:2])
  res2 <- fullGamma(block2, c(2))
  expect_equal(res1, Gamma3[1:2, 1:2])
  expect_equal(res2, Gamma3[1:2, 1:2])
})

test_that("Gamma2Graph works", {
  res1 <- fullGamma(block, Gamma3)
  expect_s3_class(Gamma2Graph(res1), "igraph")
  expect_s3_class(Gamma2Graph(res1, to_plot = T), "igraph")
  expect_s3_class(Gamma2Graph(res1, to_plot = F), "igraph")
})

test_that("data2rmpareto works", {
  p <- .999
  m <- 1 / (1 - apply(data, 2, unif))
  q <- 1 / (1 - p)
  idx <- which(apply(m, 1, max) > q)
  res <- m[idx, ] / q
  expect_equal(data2mpareto(data, p), res)

  p <- .95
  m <- 1 / (1 - apply(data, 2, unif))
  q <- 1 / (1 - p)
  idx <- which(apply(m, 1, max) > q)
  res <- m[idx, ] / q
  expect_equal(data2mpareto(data, p), res)

  p <- 0
  m <- 1 / (1 - apply(data, 2, unif))
  q <- 1 / (1 - p)
  idx <- which(apply(m, 1, max) > q)
  res <- m[idx, ] / q
  expect_equal(data2mpareto(data, p), res)

})

test_that("Sigma2Gamma works", {
  for (k in 1:NCOL(G)){
    S <- Gamma2Sigma(G, k = k, full = F)
    expect_equal(Sigma2Gamma(S = S, k = k, full = F), G)
  }

  for (k in 1:NCOL(G)){
    S <- Gamma2Sigma(G, k = k, full = T)
    expect_equal(Sigma2Gamma(S, k = k, full = T), G)
    expect_equal(Sigma2Gamma(S, full = T), G)
  }
})

test_that("Gamma2Sigma works", {
  d <- sample(2:10, 1)
  S <- matrix(runif(d ^ 2), nrow = d, ncol = d) + diag(d)

  for (k in 1:d){
    G <- Sigma2Gamma(S = S, k = k, full = F)
    expect_equal(Gamma2Sigma(Gamma = G, k = k, full = F), S)
  }

  S <- cbind(rep(0, d + 1), rbind(rep(0, d), S))
  G <- Sigma2Gamma(S = S, full = T)
  expect_equal(Gamma2Sigma(Gamma = G, full = T), S)

  for (k in 1:d){
    shuffle <- 1:(d + 1)
    shuffle[shuffle <= k] <- shuffle[shuffle <= k] - 1
    shuffle[1] <- k
    shuffle <- order(shuffle)
    S_shuffled <- S[shuffle, shuffle]
    G <- Sigma2Gamma(S = S_shuffled, k = 20, full = T)
    expect_equal(Gamma2Sigma(Gamma = G, k = k, full = T), S_shuffled)
  }
})

test_that("par2Gamma works", {
  G1 <- matrix(0, nrow = 4, ncol = 4)
  G1[upper.tri(G1)] <- par
  G1 <- G1 + t(G1)
  expect_error(par2Gamma(par = 1:4))
  expect_equal(par2Gamma(par = par), G1)
})

test_that("Gamma2par works", {
  expect_equal(Gamma2par(Gamma = G), G[upper.tri(G)])
  expect_equal(Gamma2par(Gamma = c(1.5, 2, 1.5)), c(1.5, 2, 1.5))
})

test_that("Chi2Gamma works", {
  chi <- runif(1)
  expect_error(Chi2Gamma(2))
  expect_equal(Chi2Gamma(0), Inf)
  expect_equal(Chi2Gamma(1), 0)
  expect_equal(Chi2Gamma(chi),  (2 * qnorm(1 - 0.5 * chi)) ^ 2)

  theta <- 2 - chi
  expect_equal(Chi2Gamma(chi), (2*qnorm(theta/2))^2)
})
