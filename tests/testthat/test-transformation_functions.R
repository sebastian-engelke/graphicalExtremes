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

G3 <- matrix(nrow = 7, ncol = 7)
for (i in 1:7){
  for (j in 1:7){
    G3[i, j] <- abs(i - j)/2
  }
}


# Run tests
test_that("complete_Gamma works", {

  expect_error(complete_Gamma(Gamma1, non_decomposable))
  expect_error(complete_Gamma(Gamma2, non_block))
  expect_warning(complete_Gamma(Gamma3, igraph::as.directed(block)))
  expect_error(complete_Gamma(Gamma3, empty_graph))
  expect_error(complete_Gamma(Gamma3_wrong, block))
  expect_error(complete_Gamma(Gamma3_vec_wrong, block))
  expect_error(complete_Gamma(Gamma2, block))

  res1 <- complete_Gamma(Gamma3, block)
  res2 <- complete_Gamma(Gamma3_vec, block)
  expect_equal(Gamma3_completed, res1)
  expect_equal(Gamma3_completed, res2)

  res1 <- complete_Gamma(Gamma3[1:2, 1:2], block2)
  res2 <- complete_Gamma(c(2), block2)
  expect_equal(Gamma3[1:2, 1:2], res1)
  expect_equal(Gamma3[1:2, 1:2], res2)
})

test_that("Gamma2graph works", {
  res1 <- complete_Gamma(Gamma3, block)
  expect_s3_class(Gamma2graph(res1), "igraph")
  expect_s3_class(Gamma2graph(res1, to_plot = T), "igraph")
  expect_s3_class(Gamma2graph(res1, to_plot = F), "igraph")
})

test_that("data2rmpareto works", {
  p <- .999
  m <- 1 / (1 - apply(data, 2, unif))
  q <- stats::quantile(m, p)
  idx <- which(apply(m, 1, max) > q)
  res <- m[idx, ] / q
  expect_equal(data2mpareto(data, p), res)

  p <- .95
  m <- 1 / (1 - apply(data, 2, unif))
  q <- stats::quantile(m, p)
  idx <- which(apply(m, 1, max) > q)
  res <- m[idx, ] / q
  expect_equal(data2mpareto(data, p), res)

  p <- 0
  m <- 1 / (1 - apply(data, 2, unif))
  q <- stats::quantile(m, p)
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

test_that("chi2Gamma works", {
  chi <- runif(1)
  expect_error(chi2Gamma(2))
  expect_equal(chi2Gamma(0), Inf)
  expect_equal(chi2Gamma(1), 0)
  expect_equal(chi2Gamma(chi),  (2 * stats::qnorm(1 - 0.5 * chi)) ^ 2)

  chi <- matrix(runif(4 * 4), nrow = 4, ncol = 4)
  diag(chi) <- 1
  res <- chi2Gamma(chi)
  expect_equal(res, (2 * stats::qnorm(1 - 0.5 * chi)) ^ 2)

  theta <- 2 - chi
  expect_equal(chi2Gamma(chi), (2*stats::qnorm(theta/2))^2)
})

test_that("Gamma2chi works", {
  chi <- runif(1)
  expect_error(chi2Gamma(2))
  expect_equal(chi2Gamma(0), Inf)
  expect_equal(chi2Gamma(1), 0)
  expect_equal(chi2Gamma(chi),  (2 * stats::qnorm(1 - 0.5 * chi)) ^ 2)

  theta <- 2 - chi
  expect_equal(chi2Gamma(chi), (2*stats::qnorm(theta/2))^2)
})

test_that("Gamma2chi_3D works", {
  expect_error(Gamma2chi_3D(G3))
  expect_gte(Gamma2chi_3D(1e16 * G3[1:3, 1:3]), 0)
  expect_lte(Gamma2chi_3D(1e-16 * G3[1:3, 1:3]), 1)
})

