context("test-simulation_helpers")

# Define variables
n <- 7
idx <- 2
d <- 4
trend <- c(0.75, 0.00, 1.00, 0.75)
G <-  cbind(c(0, 1.5, 1.5, 2),
            c(1.5, 0, 2, 1.5),
            c(1.5, 2, 0, 1.5),
            c(2, 1.5, 1.5, 0))
G.vec <- c(2, 3)
A <- list(
  rbind(c(0, 0),
        c(0, 1),
        c(1, 1))
)
# !!! continue with building A, tests with simulation_helpers, and
# change explanation of A in the function descriptions

cov_mat <- Gamma2Sigma(G, k = 3, full = FALSE)
chol_mat <- matrix(0, d, d)
chol_mat[-3, -3] <- chol(cov_mat)


# Run tests
test_that("simu_px_HR works", {
  expect_error(simu_px_HR(n, c(1, 3), d, trend, chol_mat))

  res <- simu_px_HR(n, idx, d, trend, chol_mat)
  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))
})

test_that("simu_px_logistic works", {
  expect_error(simu_px_logistic(n, idx = c(1, 3), d, theta = 0.2))
  expect_error(simu_px_logistic(n, idx = c(sample(1:d, n, replace = F), 1),
                                d, theta = 0.2))

  res <- simu_px_logistic(n, idx, d, theta = 0.2)
  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))

  res <- simu_px_logistic(7, sample(x = 1:4, size = 7, replace = T), d,
                          theta = 0.2)

  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))
})

test_that("simu_px_neglogistic works", {
  expect_error(simu_px_neglogistic(n, idx = c(1, 3), d, theta = 0.2))
  expect_error(simu_px_neglogistic(n, idx = c(sample(1:d, n, replace = T), 1),
                                d, theta = 0.2))


  res <- simu_px_neglogistic(n, idx, d, theta = 0.2)
  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))

  res <- simu_px_neglogistic(7, sample(x = 1:4, size = 7, replace = T), d,
                          theta = 0.2)

  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))

})

test_that("simu_px_dirchlet works", {
  expect_error(simu_px_dirichlet(n, idx = c(1, 3), d,
                                 alpha = c(0.2, 1, 1.2, 0.8)))
  expect_error(simu_px_dirichlet(n, idx = c(sample(1:d, n, replace = T), 1),
                                   d, alpha = c(0.2, 1, 1.2, 0.8)))


  res <- simu_px_dirichlet(n, idx, d, alpha = c(0.2, 1, 1.2, 0.8))
  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))

  res <- simu_px_dirichlet(7, sample(x = 1:4, size = 7, replace = T), d,
                           alpha = c(0.2, 1, 1.2, 0.1))

  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))
})


test_that("simu_px_tree_HR works", {
  res <- simu_px_tree_HR(n, G.vec = )
  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))
})
