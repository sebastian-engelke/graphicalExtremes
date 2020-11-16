rmse = function(G, Ghat){
  idx = which(Ghat[upper.tri(Ghat)] !=Inf)
  sqrt(mean(((G[upper.tri(G)][idx] -  Ghat[upper.tri(Ghat)][idx])/ G[upper.tri(G)][idx])^2))
}

rand.tree <- function(d){
  # Generate a random tree with dimension d
  adj <- matrix(FALSE, d,d)
  i0 <- sample(1:d, size = 1)
  j0 <- sample(setdiff(1:d,i0), size = 1)
  adj[i0,j0] <- adj[j0,i0] <- TRUE
  
  while(length(which(rowSums(adj)==0)) > 0){
    idx.i <- which(rowSums(adj)<=2 & rowSums(adj)>=1)
    if(length(idx.i)==1) i <- idx.i
    else i <- sample(idx.i, 1)
    idx.j <- which(rowSums(adj)==0)
    if(length(idx.j)==1) j <- idx.j
    else j <- sample(idx.j, 1)
    adj[i,j] <- adj[j,i] <- TRUE
  }
  
  diag(adj) <- TRUE
  Tree <- graph_from_adjacency_matrix(adj, diag =FALSE, mode="undirected")
  #plot(Tree)
  return(list(adj=adj, Tree=Tree))
}




sim_study <- function(d = 10, n = 100, p = NULL, 
                      model = c("dirichlet", "HR",
                                "logistic"), method = c("maxstable", "mpareto"), 
                      edge.max = 0.3, ML = FALSE, cens = TRUE, chi.edge = NULL,
                      noise = c("none", "iid", "tree"), rng = NULL){
  ## perform a simulation study for tree learning with different configurations.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - n: sample size.
  ##    - p: if NULL, then no thresholding is performed and the data is assumed 
  ##       to be MGPD; if not NULL then this is the probability threshold that is
  ##       used to define exceedances.
  ##    - model: "HR", "dirichlet", or "logistic".
  ##       The edge weights for HR are defined using 
  ##       edge.max and for dirichlet it is defined in the code. Note that 
  ##       dirichlet is an asymmetric model.
  ##    - method: "mpareto" or "maxstable".
  ##    - edge.max: parameter to sample the HR tree parameters.
  ##    - ML: logical. If TRUE, ML tree estimation is perforemd, with or without 
  ##       censoring (this is slow).
  ##    - cens: logical. If TRUE, censoring is performed in ML trees.
  ##    - chi.edge: numeric.
  ##    - noise: logical.
  ##    - rng: integer, 7-element vector or NULL
  ## Returns:
  ##  a tibble with error, srr, and rmse
  
  # check arguments
  model <- match.arg(model)
  method <- match.arg(method)
  noise <- match.arg(noise)
  
  
  if(noise=="tree"){
    set.seed(34123)
    tree_noise <- rand.tree(d)$Tree
    if(model=="HR")
      par_noise <- complete_Gamma(graph=tree_noise, Gamma = runif(n=length(E(tree_noise)), 
                                                            min = 0.2, 
                                                            max = edge.max))
    else if(model=="dirichlet")
      par_noise = matrix(runif(n = 2*(d-1), min = 1, max = 10), nrow=d-1, ncol=2) 
  }
    
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  # perform simulation
  err.mat <- numeric(4)
  rmse.mat <- numeric(4)
  srr.mat <- numeric(4)
  time.mat <- numeric(4)
  
  graph.full <- make_full_graph(d)
  G <- NULL
  alpha.mat <- NULL
  
  tree <- rand.tree(d = d)$Tree
  
  if(model=="HR"){
    if(!is.null(chi.edge))
      par <- complete_Gamma(graph=tree, Gamma = rep(chi2Gamma(chi.edge),
                                                    times=d-1))
    else
      par <- complete_Gamma(graph=tree, Gamma = runif(n=length(E(tree)), 
                                                      min = 0.2, 
                                                      max = edge.max)) 
  }
  else if(model=="dirichlet"){
    par = matrix(runif(n = 2*(d-1), min = 1, max = 10), nrow=d-1, ncol=2) 
  }
  else if(model=="logistic"){
    if(!is.null(chi.edge))
      par <-rep(log(2-chi.edge)/log(2), times=d-1) 
    else
      par <- runif(d-1, min = .4, max = .6)
  }
  
  if(method=="mpareto")
    X <- rmpareto_tree(n=n, model=model, tree=tree, par=par)
  else if(method=="maxstable")
    X <- rmstable_tree(n=n, model=model, tree=tree, par=par)
  
  if(noise=="iid")
    Y <- matrix((-1/log(runif(d*n))), nrow=n, ncol=d)
  else if(noise=="tree")
    Y <- rmstable_tree(n=n, model=model, tree=tree_noise, par=par_noise)
  else if(noise=="none")
    Y=0
  
  X <- X + Y^{0.5}
  
  ptm <- proc.time()
  chi.est <- emp_chi_mat(data=X, p=p)
  tree.chi.est <- igraph::mst(graph=graph.full, 
                              weights = 2-chi.est[ends(graph.full,
                                                       E(graph.full))],
                              algorithm = "prim") 
  time.mat[1] <- proc.time() - ptm
  
  
  ptm <- proc.time()
  G1.est <- graphicalExtremes:::emp_vario(data=X, k=1, p=p)
  tree.G1.est <- igraph::mst(graph=graph.full, 
                             weights = 2-Gamma2chi(G1.est[ends(graph.full,
                                                               E(graph.full))]),
                             algorithm = "prim")
  time.mat[2] <- proc.time() - ptm
  
  
  ptm <- proc.time()
  G.est <- graphicalExtremes:::emp_vario(data=X, p=p) 
  tree.G.est <- igraph::mst(graph=graph.full, 
                            weights = 2- Gamma2chi(G.est[ends(graph.full,
                                                              E(graph.full))]),
                            algorithm = "prim")
  time.mat[3] <- proc.time() - ptm
  
  if(ML){
    ptm <- proc.time()
    xx <- graphicalExtremes::mst_HR(data = X, cens = cens, p = p)
    time.mat[4] <- proc.time() - ptm
    G.ML.est <- xx$Gamma
    tree.ML.est <- xx$tree
  }
    
  
  ## compute errors
  err.mat[1] <- 1-ecount(intersection(tree.chi.est, tree))/(d-1)
  err.mat[2] <- 1-ecount(intersection(tree.G1.est, tree))/(d-1)
  err.mat[3] <- 1-ecount(intersection(tree.G.est, tree))/(d-1)
  if(ML){
    err.mat[4] <- 1-ecount(intersection(tree.ML.est, tree))/(d-1)
  }
  
  srr.mat <- (err.mat != 0) * 1
  
  if(model=="HR"){
    rmse.mat[1] = rmse(par,chi2Gamma(chi.est))
    rmse.mat[2] = rmse(par,G1.est)
    rmse.mat[3] = rmse(par,G.est)
    if (ML){
      rmse.mat[4] = rmse(par,G.ML.est)
    }
  }
  
  tbl <- tibble(type = paste0("err", 1:length(err.mat)),
                value = err.mat) %>% 
    bind_rows(tibble(type = paste0("rmse", 1:length(rmse.mat)),
                     value = rmse.mat)) %>% 
    bind_rows(tibble(type = paste0("srr", 1:length(srr.mat)),
                     value = srr.mat)) %>% 
    bind_rows(tibble(type = paste0("time", 1:length(time.mat)),
                     value = time.mat))
  
  return(tbl)
}


wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Args:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_study, args = fun_args[i, ]) %>% 
    mutate(rowid = rowid)
}

set_rng <- function(tbl, seed){
  ## adds to tbl a column with seeds to generate independent streams of random
  ## numbers.
  ##
  ## Args:
  ##     - tbl: a tibble where the columns contain the parameter settings and the
  ##       rows contain the simulation runs.
  ##
  ## Returns:
  ##     The function returns tbl appending a column with seeds used to generate
  ##     independent streams of random numbers.
  ##
  ## Note:
  ##     This function creates ensures that the simulations are fully repeatable.
  ##     This is possible because it assigns to each simulation run a unique 
  ##     random seed (generated with L'Ecuyer RNG method, which is suitable 
  ##     for parallel processing, too).
  
  m <- n_groups(tbl)
  group_idxs <- group_indices(tbl)
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m, seed = seed, simplify = FALSE)
  rng <- rng[group_idxs]
  
  # add RNG streams to tbl
  tbl$rng <- rng
  
  # return tibble
  return(tbl)
}

rep_tibble <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    select(-rowname)
  
} 

assign_random_seed <- function(tbl, grouping_vars, seed){
  ## tibble character_vector integer -> tibble
  ## assign random seed according to the variables in grouping_vars
  if (is.null(grouping_vars)){
    tbl <- tbl %>%
      rowwise()
  } else {
    tbl <- tbl %>% 
      group_by(across(all_of(grouping_vars)))
  }
  
  tbl %>% 
    set_rng(seed) %>% 
    ungroup()
}


