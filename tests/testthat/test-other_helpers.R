context("test-other_helpers")

# Define variables 1
set.seed(1234)
x <- c(1, 5, -1, 2)
G1 <- matrix(nrow = 1, ncol = 2)
G2 <- matrix(nrow = 1, ncol = 1)
G3 <- matrix(nrow = 3, ncol = 3)


# Run tests
test_that("unif works", {
  expect_equal(unif(x), c(2 / 5, 4 / 5, 1 / 5, 3 / 5))
})

test_that("dim_Gamma works", {
  expect_error(dim_Gamma(G1))
  expect_equal(dim_Gamma(G2), 1)
  expect_equal(dim_Gamma(G3), 3)
})

test_that("select_edges works", {
  small_graph <- igraph::make_empty_graph(n = 4, directed = FALSE)
  small_graph <- igraph::add_edges(small_graph, c(1, 2, 2, 3, 2, 4))
  expect_equal(select_edges(small_graph), rbind(c(1, 3), c(1, 4), c(3, 4)))
})

test_that("fast_diag works", {
  n <- 2e3
  d <- 1e2
  y <- matrix(runif(n * d), nrow = n)
  M <- matrix(runif(d * d), nrow = d)
  expect_equal(fast_diag(y, M), diag(y %*% M %*% t(y)))
})

test_that("graphs_equal works", {
  g1 <- generate_random_connected_graph(d = 5)
  g2 <- igraph::add_edges(g1, c(1, 2, 2, 3, 2, 4))

  expect_equal(graphs_equal(g1, g2), FALSE)
  expect_equal(graphs_equal(g1, g1), TRUE)
})


