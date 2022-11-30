

p1 <- 0.3
p2 <- 0.4

makeP <- function(d){
  diag(d) - matrix(1/d, d, d)
}

# T <- rbind(
#   c(p1, -p1, 0),
#   c(-p1, p1+p2, -p2),
#   c(0, -p2, p2)
# )

T1 <- matrix(0, 3, 3)
T1[1:2, 1:2] <- rbind(
  c(p1, -p1),
  c(-p1, p1)
)

T2 <- matrix(0, 3, 3)
T2[1:2, 1:2] <- rbind(
  c(p2, -p2),
  c(-p2, p2)
)


P1 <- matrix(0, 3, 3)
P1[1:2, 1:2] <- makeP(2)

T <- T1 + T2


