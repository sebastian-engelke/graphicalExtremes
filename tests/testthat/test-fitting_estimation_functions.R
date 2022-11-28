context("test-fitting_estimation_functions")

# Define variables
set.seed(1234)
n <- 1e3
G1 <- cbind(
  c(0, 1.5, 1.5, 2),
  c(1.5, 0, 2, 1.5),
  c(1.5, 2, 0, 1.5),
  c(2, 1.5, 1.5, 0)
)
data1 <- rmpareto(n = n, model = "HR", d = NCOL(G1), par = G1)

G2 <- G1[-(3:4), -(3:4)]
data2 <- rmpareto(n = n, model = "HR", d = NCOL(G2), par = G2)

d <- 7
G3 <- matrix(nrow = 7, ncol = 7)
for (i in 1:7) {
  for (j in 1:7) {
    G3[i, j] <- abs(i - j) / 2
  }
}

data3 <- rmpareto(n = n, model = "HR", d = NCOL(G3), par = G3)

non_decomposable <- igraph::graph_from_adjacency_matrix(
  rbind(
    c(0, 1, 0, 0, 0, 0),
    c(1, 0, 1, 0, 1, 0),
    c(0, 1, 0, 1, 0, 1),
    c(0, 0, 1, 0, 1, 1),
    c(0, 1, 0, 1, 0, 0),
    c(0, 0, 1, 1, 0, 0)
  ),
  mode = "undirected"
)
igraph::is_chordal(non_decomposable)$chordal

non_block <- igraph::graph_from_adjacency_matrix(
  rbind(
    c(0, 1, 0, 0, 0, 0),
    c(1, 0, 1, 1, 1, 0),
    c(0, 1, 0, 1, 1, 1),
    c(0, 1, 1, 0, 1, 1),
    c(0, 1, 1, 1, 0, 0),
    c(0, 0, 1, 1, 0, 0)
  ),
  mode = "undirected"
)
igraph::is_chordal(non_block)$chordal

block <- igraph::graph_from_adjacency_matrix(
  rbind(
    c(0, 1, 0, 0, 0, 0, 0),
    c(1, 0, 1, 1, 1, 0, 0),
    c(0, 1, 0, 1, 1, 1, 1),
    c(0, 1, 1, 0, 0, 0, 0),
    c(0, 1, 1, 1, 0, 0, 0),
    c(0, 0, 1, 0, 0, 0, 1),
    c(0, 0, 1, 0, 0, 1, 0)
  ),
  mode = "undirected"
)
igraph::is_chordal(block)$chordal

non_connected <- block
non_connected[1, 2] <- non_connected[2, 1] <- 0

Gamma3_completed <- rbind(
  c(0, 2, 4, 4, 4, 6, 6),
  c(2, 0, 2, 2, 2, 4, 4),
  c(4, 2, 0, 2, 2, 2, 2),
  c(4, 2, 2, 0, 2, 4, 4),
  c(4, 2, 2, 2, 0, 4, 4),
  c(6, 4, 2, 4, 4, 0, 2),
  c(6, 4, 2, 4, 4, 2, 0)
)


# Run tests
test_that("emp_chi_multdim works", {
  expect_error(emp_chi_multdim(data = data1[, 1], p = .95))
  expect_error(emp_chi_multdim(data = as.matrix(data1[, 1]), p = .95))
  expect_length(emp_chi_multdim(data = data2, p = .95), 1)
})


test_that("emp_chi works", {
  expect_error(emp_chi(data = data1[, 1], p = .95))
  expect_error(emp_chi(data = as.matrix(data1[, 1]), p = .95))

  dat <- data1
  res <- emp_chi(data = dat, p = .9301)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))
  expect_equal(all(!is.na(res)), TRUE)

  dat <- data2
  res <- emp_chi(data = dat, p = .95)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))
  expect_equal(all(!is.na(res)), TRUE)


  dat <- data3
  res <- emp_chi(data = dat, p = .95)
  expect_equal(NROW(res), NCOL(dat))
  expect_equal(NCOL(res), NCOL(dat))
  expect_equal(all(!is.na(res)), TRUE)

  dat <- data2mpareto(data1, p = 0.95)
  expect_equal(emp_chi(data = dat), emp_chi(data1, p = 0.95))
})

test_that("emp_vario works", {
  data <- rmpareto(1e1, "HR", d = 4, G1)
  dd <- NCOL(data)

  res <- emp_vario(data)
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0 | is.na(res)), T)

  res <- emp_vario(data, k = sample(1:dd, 1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0 | is.na(res)), T)

  res <- emp_vario(data, p = runif(1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0 | is.na(res)), T)

  res <- emp_vario(data, p = .99)
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0 | is.na(res)), T)

  res <- emp_vario(data, k = sample(1:dd, 1), p = runif(1))
  expect_equal(NROW(res), dd)
  expect_equal(NCOL(res), dd)
  expect_type(res, "double")
  expect_equal(all(res >= 0 | is.na(res)), T)
})

test_that("V_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  expect_error(V_HR(x = rep(0, d), Gamma = par))
  expect_error(V_HR(x = rep(1, d + 1), Gamma = par))
  res <- V_HR(x = rep(1, d), Gamma = par)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("logdV_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  expect_error(logdV_HR(x = rep(0, d), Gamma = par))
  expect_error(logdV_HR(x = rep(1, d + 1), Gamma = par))
  res <- logdV_HR(x = rep(1, d), Gamma = par)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("logdVK_HR works", {
  d <- NROW(G1)
  par <- G1[upper.tri(G1)]
  expect_error(logdVK_HR(
    x = rep(0, d), K = sample(1:d, ceiling(runif(1) * d)),
    Gamma = par
  ))
  expect_error(logdVK_HR(
    x = rep(1, d + 1),
    K = sample(1:d, ceiling(runif(1) * d)),
    Gamma = par
  ))
  res <- logdVK_HR(
    x = rep(1, d), K = sample(1:d, ceiling(runif(1) * (d - 1))),
    Gamma = par
  )
  expect_type(res, "double")
  expect_length(res, 1)

  d <- NROW(G2)
  par <- G2[upper.tri(G2)]
  res <- logdVK_HR(
    x = rep(1, d), K = sample(1:d, ceiling(runif(1) * (d - 1))),
    Gamma = par
  )

  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("logLH_HR works", {
  res <- logLH_HR(data1, G1)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1[1, ], G1)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1, G1, cens = FALSE)
  expect_type(res, "double")
  expect_length(res, 1)

  res <- logLH_HR(data1, G1, cens = TRUE)
  expect_type(res, "double")
  expect_length(res, 1)
})

test_that("fmpareto_graph_HR works", {
  expect_error(fmpareto_graph_HR(
    graph = non_connected,
    data = data3, p = .95, cens = FALSE
  ))
  expect_error(fmpareto_graph_HR(
    graph = non_decomposable,
    data = data3, p = .95, cens = FALSE
  ))
  expect_error(fmpareto_graph_HR(
    graph = non_block,
    data = data3, p = .95, cens = FALSE
  ))
  expect_error(fmpareto_graph_HR(
    graph = block,
    data = data3[, -1], p = .95, cens = FALSE
  ))
})

test_that("fmpareto_HR_MLE works", {
  d <- 3
  v_idx <- c(2, 4, 1)
  Gamma_small_block <- Gamma3_completed[v_idx, v_idx]
  small_block <- igraph::induced_subgraph(block, v_idx)
  data <- rmstable(1e3, "HR", d = d, Gamma_small_block)

  # providing p
  init_param <- Gamma_small_block[upper.tri(Gamma_small_block)]
  res <- fmpareto_HR_MLE(
    data = data,
    p = 0.95,
    cens = FALSE,
    useTheta = FALSE,
    init = init_param
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)

  res <- fmpareto_HR_MLE(
    data = data, p = 0.95, cens = TRUE,
    useTheta = FALSE,
    init = init_param
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)

  init_param <- c(2, 2)
  res <- fmpareto_HR_MLE(
    data = data, p = 0.95, cens = FALSE,
    useTheta = FALSE,
    init = init_param, graph = small_block
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)


  res <- fmpareto_HR_MLE(
    data = data, p = 0.95, cens = TRUE,
    useTheta = FALSE,
    init = init_param, graph = small_block
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)

  set.seed(19)
  data <- rmpareto(1e3, "HR", d = d, Gamma_small_block)

  # not providing p
  init_param <- Gamma_small_block[upper.tri(Gamma_small_block)]
  res <- fmpareto_HR_MLE(
    data = data, p = NULL, cens = FALSE,
    useTheta = FALSE,
    init = init_param
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)

  res <- fmpareto_HR_MLE(
    data = data, p = NULL, cens = TRUE,
    useTheta = FALSE,
    init = init_param
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)

  init_param <- c(2, 2)
  res <- fmpareto_HR_MLE(
    data = data, p = NULL, cens = FALSE,
    useTheta = FALSE,
    init = init_param, graph = small_block
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)


  res <- fmpareto_HR_MLE(
    data = data, p = NULL, cens = TRUE,
    useTheta = FALSE,
    init = init_param, graph = small_block
  )
  expect_type(res, "list")
  expect_equal(t(res$Gamma), res$Gamma)
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(t(res$hessian), res$hessian)
})


test_that("emst works", {
  d <- 3
  v_idx <- c(2, 4, 1)
  Gamma_small_block <- Gamma3_completed[v_idx, v_idx]
  data <- rmpareto(1e3, "HR", d = d, Gamma_small_block)

  res <- emst(data = data, p = NULL, cens = FALSE)
  expect_length(res, 2)
  expect_equal(class(res$graph), "igraph")
  expect_equal(res$Gamma, t(res$Gamma))
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(is.numeric(res$Gamma) & is.matrix(res$Gamma), TRUE)

  res <- emst(data = data, p = 0.95, cens = FALSE)
  expect_length(res, 2)
  expect_equal(class(res$graph), "igraph")
  expect_equal(res$Gamma, t(res$Gamma))
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(is.numeric(res$Gamma) & is.matrix(res$Gamma), TRUE)

  res <- emst(data = data, p = NULL, cens = TRUE)
  expect_length(res, 2)
  expect_equal(class(res$graph), "igraph")
  expect_equal(res$Gamma, t(res$Gamma))
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(is.numeric(res$Gamma) & is.matrix(res$Gamma), TRUE)

  res <- emst(data = data, p = 0.95, cens = TRUE)
  expect_length(res, 2)
  expect_equal(class(res$graph), "igraph")
  expect_equal(res$Gamma, t(res$Gamma))
  expect_equal(all(!is.na(res$Gamma)), TRUE)
  expect_equal(is.numeric(res$Gamma) & is.matrix(res$Gamma), TRUE)
})

test_that("loglik_HR works", {
  d <- 3
  v_idx <- c(2, 4, 1)
  Gamma_small_block <- Gamma3_completed[v_idx, v_idx]
  data <- rmpareto(1e3, "HR", d = d, Gamma_small_block)

  expect_named(
    loglik_HR(data, graph = block, Gamma = Gamma_small_block, cens = FALSE),
    c("loglik", "aic", "bic")
  )
})


test_that("eglearn works", {
  d <- 4
  my_model <- generate_random_model(d = d)
  data <- rmpareto(1e3, "HR", d = d, my_model$Gamma)
  res1 <- eglearn(data, reg_method = "ns",
                 complete_Gamma = FALSE)

  res2 <- eglearn(data, reg_method = "ns",
                  complete_Gamma = FALSE)

  res3 <- eglearn(data, reg_method = "glasso",
                  complete_Gamma = FALSE)

  expect_length(res1, 5)
  expect_equal(class(res1$graph[[1]]), "igraph")
  expect_equal(res3$Gamma[[1]], NA)
  expect_error(eglearn(data, rholist = -1))
  expect_error(eglearn(data, rholist = c(2, -1)))
  expect_message(eglearn(data, rholist = c(1, 200),
                         complete_Gamma = TRUE))

})

