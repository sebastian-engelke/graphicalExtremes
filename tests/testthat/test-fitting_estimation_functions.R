context("test-fitting_estimation_functions")


# Define variables
set.seed(1234)
n <- 1e3
G1 <-  cbind(c(0, 1.5, 1.5, 2),
             c(1.5, 0, 2, 1.5),
             c(1.5, 2, 0, 1.5),
             c(2, 1.5, 1.5, 0))
data1 <- rmpareto(n = n, model = "HR", d = NCOL(G1), par = G1)
G2 <- G1[-(3:4), -(3:4)]
data2 <- rmpareto(n = n, model = "HR", d = NCOL(G2), par = G2)

d <- 7
G3 <- matrix(nrow = 7, ncol = 7)
for (i in 1:7){
  for (j in 1:7){
    G3[i, j] <- abs(i - j)/2
  }
}

data3 <- rmpareto(n = n, model = "HR", d = NCOL(G3), par = G3)

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

non_connected <- block
non_connected[1, 2] <- non_connected[2, 1] <- 0

block2 <- igraph::graph_from_adjacency_matrix(
  rbind(c(0, 1),
        c(1, 0)),
  mode = "undirected")
igraph::is_chordal(block2)$chordal

Gamma3_completed <- rbind(c(0, 2, 4, 4, 4, 6, 6),
                          c(2, 0, 2, 2, 2, 4, 4),
                          c(4, 2, 0, 2, 2, 2, 2),
                          c(4, 2, 2, 0, 2, 4, 4),
                          c(4, 2, 2, 2, 0, 4, 4),
                          c(6, 4, 2, 4, 4, 0, 2),
                          c(6, 4, 2, 4, 4, 2, 0))


# Run tests
test_that("emp_chi works", {
  expect_error(emp_chi(data = data1[, 1], p = .95, pot = F))
  expect_error(emp_chi(data = as.matrix(data1[, 1]), p = .95, pot = T))
  expect_length(emp_chi(data = data2, p = .95, pot = F), 1)
  expect_length(emp_chi(data = data2, p = .95, pot = T), 1)
})

test_that("emp_chi_mat works", {
  dat <- data1
  res <- emp_chi_mat(data = dat, p = .95, pot = T)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data1
  res <- emp_chi_mat(data = dat, p = .95, pot = F)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data2
  res <- emp_chi_mat(data = dat, p = .95, pot = T)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data2
  res <- emp_chi_mat(data = dat, p = .95, pot = F)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data3
  res <- emp_chi_mat(data = dat, p = .95, pot = T)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data3
  res <- emp_chi_mat(data = dat, p = .95, pot = F)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))
})

test_that("emp_vario works", {
  data <- rmpareto(1e1, "HR", d = 4, G1)
  dd <- NCOL(data)

  res <- emp_vario(data)
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- emp_vario(data, k = sample(1:dd, 1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- emp_vario(data, p = runif(1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- emp_vario(data, p = .99)
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- emp_vario(data, k = sample(1:dd, 1), p = runif(1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

})

test_that("V_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  expect_error(V_HR(x = rep(0, d), par = par))
  expect_error(V_HR(x = rep(1, d + 1), par = par))
  res <- V_HR(x = rep(1, d), par = par)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("logdV_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  expect_error(logdV_HR(x = rep(0, d), par = par))
  expect_error(logdV_HR(x = rep(1, d + 1), par = par))
  res <- logdV_HR(x = rep(1, d), par = par)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("logdVK_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  expect_error(logdVK_HR(x = rep(0, d), K = sample(1:d, ceiling(runif(1) * d)),
                         par = par))
  expect_error(logdVK_HR(x = rep(1, d + 1),
                         K = sample(1:d, ceiling(runif(1) * d)),
                         par = par))
  res <- logdVK_HR(x = rep(1, d), K = sample(1:d, ceiling(runif(1) * d)),
                   par = par)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("logLH_HR works", {

  res <- logLH_HR(data1, G1)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1, G1, cens = FALSE)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1, G1, cens = TRUE)
  expect_type(res, "double")
  expect_length(res, 1)

})

test_that("fmpareto_graph_HR works", {
  expect_warning(fmpareto_graph_HR(graph = igraph::as.directed(block),
                             data = data3, p = .95, cens = FALSE))
  # expect_warning(fmpareto_graph_HR(graph = igraph::as.directed(block),
                             # data = data3, p = .95, cens = TRUE))
  expect_error(fmpareto_graph_HR(graph = non_connected,
                           data = data3, p = .95, cens = FALSE))
  expect_error(fmpareto_graph_HR(graph = non_decomposable,
                           data = data3, p = .95, cens = FALSE))


  data <- rmpareto(1e3, "HR", d = 2, Gamma3_completed[1:2, 1:2])
  res <- fmpareto_graph_HR(block2, data, p = .95)
  expect_type(res, "list")
  expect_length(res, 2)
})

test_that("fmpareto_HR works", {

})

test_that("mst_HR works", {

})
