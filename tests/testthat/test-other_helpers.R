context("test-other_helpers")

# Define variables 1
set.seed(1234)
x <- c(1, 5, -1, 2)
G1 <- matrix(nrow = 1, ncol = 2)
G2 <- matrix(nrow = 1, ncol = 1)
G3 <- matrix(nrow = 3, ncol = 3)


# Run tests
test_that("unif works", {
  expect_equal(unif(x), c(2/5, 4/5, 1/5, 3/5))
})

test_that("dim_Gamma works", {
  expect_error(dim_Gamma(G1))
  expect_equal(dim_Gamma(G2), 1)
  expect_equal(dim_Gamma(G3), 3)
})
