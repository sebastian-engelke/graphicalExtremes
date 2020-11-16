library(Rcpp)
sourceCpp("src/logdv_helper.cpp")

d <- 30
n <- 100e2
x <- matrix(runif(n * d), ncol = d)
G <- matrix(nrow = d, ncol = d)
for (i in 1:d){
  for (j in 1:d){
    G[i, j] <- 2 * abs(i - j)
  }
}
S <- Gamma2Sigma(G, k=i)
cholS <- chol(S)
Sm1 <- chol2inv(cholS)
logdetS <- 2*sum(log(diag(cholS)))
y <- (t(t(log(x/x[,i])) + G[,i]/2))[,-i, drop = FALSE]


# logdv <- - apply(log(x),1,sum) - log(x[,i]) -
  # ((d-1)/2)*log(2*pi) -1/2*logdetS - 1/2 * diag(y%*%Sm1%*%t(y))

bench::mark(
  diag(y%*%Sm1%*%t(y)),
  diag_cpp(y, Sm1),
  check = FALSE
)

# !!! try with Eigen or loop?
test_that("logdv_helper works", {
  expect_equal(2 * 2, 4)
})
