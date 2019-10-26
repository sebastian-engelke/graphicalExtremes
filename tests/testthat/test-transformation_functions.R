context("test-transformation_function")


# Define variables
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


})

test_that("Gamma2Graph works", {

})
