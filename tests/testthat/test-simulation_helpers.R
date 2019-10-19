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
  expect_error(simu_px_logistic(n, idx = c(sample(1:n, n), 1), d, theta = 0.2))

  res <- simu_px_logistic(n, idx, d, theta = 0.2)
  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))

  res <- simu_px_logistic(7, 1:7, d, theta = 0.2)

  expect_type(res, "double")
  expect_equal(dim(res), c(n, d))
})
