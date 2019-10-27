context("test-fitting_estimation_functions")


# Define variables
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
