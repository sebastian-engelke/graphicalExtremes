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
test_that("chi.est works", {
  expect_warning(chi.est(data = data1$res, u = .95, pot = F))
  expect_warning(chi.est(data = data1$res, u = .95, pot = T))
  expect_length(chi.est(data = data2$res, u = .95, pot = F), 1)
  expect_length(chi.est(data = data2$res, u = .95, pot = T), 1)
})

test_that("chi_mat works", {
  dat <- data1$res
  res <- chi_mat(data = dat, u = .95, pot = T)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data1$res
  res <- chi_mat(data = dat, u = .95, pot = F)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data2$res
  res <- chi_mat(data = dat, u = .95, pot = T)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data2$res
  res <- chi_mat(data = dat, u = .95, pot = F)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data3$res
  res <- chi_mat(data = dat, u = .95, pot = T)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))

  dat <- data3$res
  res <- chi_mat(data = dat, u = .95, pot = F)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))
})

test_that("Gamma2Chi_HR works", {
  expect_error(Gamma2Chi_HR(G3))
  expect_gte(Gamma2Chi_HR(1e16 * G3[1:3, 1:3]), 0)
  expect_lte(Gamma2Chi_HR(1e-16 * G3[1:3, 1:3]), 1)
})

test_that("vario.est works", {
  data <- rmpareto(1e1, "HR", d = 4, G1)
  dd <- NCOL(data$res)

  res <- vario.est(data$res)
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- vario.est(data$res, k = sample(1:dd, 1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- vario.est(data$res, p = runif(1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- vario.est(data$res, p = .99)
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0), T)

  res <- vario.est(data$res, k = sample(1:dd, 1), p = runif(1))
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

  res <- logLH_HR(data1$res, G1)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1$res, G1, cens = FALSE)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1$res, G1, cens = TRUE)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("estGraph_HR works", {
  data <- rmpareto(1e2, "HR", 2, Gamma3_completed)
  estGraph_HR(block2, data$res, p = .95)
})

test_that("fpareto_HR works", {

})

test_that("mst_HR works", {

})
