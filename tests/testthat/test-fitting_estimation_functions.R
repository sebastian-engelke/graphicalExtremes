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

test_that("logdVK_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  logdVK_HR(x = rep(1, d), par = par)
  expect_error(V_HR(x = rep(0, d), par = par))
  expect_error(V_HR(x = rep(1, d + 1), par = par))
  res <- V_HR(x = rep(1, d), par = par)
  expect_type(res, "double")
  expect_length(res, 1)
})
