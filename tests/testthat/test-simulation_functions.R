context("test-simulation_functions")


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
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
