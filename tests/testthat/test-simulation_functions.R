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
G_wrong <- 0.2
G_wrong2 <- matrix(rep(0, d * d), nrow = d)
G_wrong3 <- G[, -1]
theta_1 <- 0.3
theta_2 <- 1.5
theta_wrong1 <- 10
theta_wrong2 <- -3
alpha <- runif(d)
alpha_wrong1 <- 43
alpha_wrong2 <- c(1, 1, 1, -0.2)
cov_mat <- Gamma2Sigma(G, k = 3, full = FALSE)
chol_mat <- matrix(0, d, d)
chol_mat[-3, -3] <- chol(cov_mat)

# Run tests
test_that("rmpareto works", {
  expect_error(rmpareto(n, "HR", -1, par = G))
  expect_error(rmpareto(n, "HR", 1.2, par = G))
  expect_error(rmpareto(-1, "HR", d, par = G))
  expect_error(rmpareto(1.2, "HR", d, par = G))
  expect_error(rmpareto(n, "HRrrrrrrrr", d, par = G))
  expect_error(rmpareto(n, "HR", d, par = G_wrong))
  expect_error(rmpareto(n, "logistic", d, par = theta_wrong1))
  expect_error(rmpareto(n, "neglogistic", d, par = theta_wrong2))
  expect_error(rmpareto(n, "dirichlet", d, par = alpha_wrong1))
  expect_error(rmpareto(n, "dirichlet", d, par = alpha_wrong2))
  expect_error(rmpareto(n, "HR", d, par = G_wrong2))
  expect_error(rmpareto(n, "HR", d, par = G_wrong3))

  res <- rmpareto(n, "HR", d, par = G)
  expect_length(res, 2)
  expect_type(res$res, "double")
  expect_type(res$counter, "double")
  expect_equal(dim(res$res), c(n, d))

  res <- rmpareto(n, "logistic", d, par = theta_1)
  expect_length(res, 2)
  expect_type(res$res, "double")
  expect_type(res$counter, "double")
  expect_equal(dim(res$res), c(n, d))

  res <- rmpareto(n, "neglogistic", d, par = theta_2)
  expect_length(res, 2)
  expect_type(res$res, "double")
  expect_type(res$counter, "double")
  expect_equal(dim(res$res), c(n, d))

  res <- rmpareto(n, "dirichlet", d, par = alpha)
  expect_length(res, 2)
  expect_type(res$res, "double")
  expect_type(res$counter, "double")
  expect_equal(dim(res$res), c(n, d))
})

test_that("rmpareto_tree works", {

})
