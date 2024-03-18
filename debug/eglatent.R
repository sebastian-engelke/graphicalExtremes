
generate_latent_model <- function(p, h) {
  W <- matrix(1, p + h, p + h)
  S <- matrix(0, p, p)
  for (i in 1:(p - 1)) {
    S[i, i + 1] <- 1
  }
  S[1, p] <- 1
  S <- S + t(S)
  diag(S) <- 0
  
  W[1:p, 1:p] <- S
  diag(W) <- 0
  L <- matrix(runif((p + h)^2, 2, 2), nrow = p + h)
  if(h!=0){
  L[(p + 1):(p + h), (p + 1):(p + h)] <- diag(h)
  
  
  z <- 30 / (sqrt(as.integer(p / h))) * matrix(runif(as.integer(p / h), 1, 2), nrow = as.integer(p / h))
  
  for (i in 1:h) {
    L[1:p, (p + i):(p + i)] <- 0
    L[(p + i):(p + i), 1:p] <- 0
    L[seq(i, p, h), (p + i):(p + i)] <- z
    L[(p + i):(p + i), seq(i, p, h)] <- z
  }
  }
  W <- W * L
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  Theta <- diag(rowSums(W)) - W
  #Lst <- (Theta[1:p, ((p + 1):(p + h))]) %*% solve(Theta[(p + 1):(p + h), (p + 1):(p + h)]) %*% t((Theta[1:p, ((p + 1):(p + h))]))
  #inc <- max(diag(svd(Lst)$u[, 1:h] %*% t(svd(Lst)$u[, 1:h])))
  G <- Theta2Gamma(Theta)
  #Lst <- (Theta[1:p, ((p + 1):(p + h))]) %*% solve(Theta[(p + 1):(p + h), (p + 1):(p + h)]) %*% t((Theta[1:p, ((p + 1):(p + h))]))
  return(list(Gamma = G, graph = Gamma2graph(G), Lst = NULL))
}




mychol <- function(M){
  d <- nrow(M)
  n <- rankMatrix(M)
  if(n==d) R <- chol(M)
  else{
    R1 <- chol(M[1:n, 1:n])
    R <- cbind(R1, solve(t(R1)) %*% M[1:n, (n+1):d])
  }
  return(R)
}

aic <- function(n, p) 2
bic <- function(n, p) log(n)

# modified BIC of Wang & Leng, JRSSB 2009
mbic <- function(n, p) log(n) * log(log(p))


BIC_score <- function(Gamma_est, graph_est, rk_est, n, Gamma_emp) {
  d <- nrow(Gamma_est)
  U <- svd(diag(d) - 1 / d * matrix(1, d, 1) %*% t(matrix(1, d, 1)))$u[, 1:as.numeric(d - 1)]
  Theta_est <- Gamma2Theta(Gamma_est)
  BIC <- (-log(det(t(U) %*% Theta_est %*% U)) - 1 / 2 * sum(diag(Gamma_emp %*% Theta_est))) + (2 * length(E(graph_est)) + 2 * rk_est * d - rk_est^2) * log(n) / (2 * n)
  return(BIC)
}







eglatent <- function(Gamma,
                     lam1_list = c(0.1, 0.15, 0.19, 0.205),
                     lam2_list = c(0.1, 0.15, 0.19, 0.205),
                     refit = TRUE) {
  d <- nrow(Gamma)
  U <- svd(diag(d) - 1 / d * matrix(1, d, 1) %*% t(matrix(1, d, 1)))$u[, 1:as.numeric(d - 1)]
  r <- 1
  Gamma_obs <- list()
  graph_obs <- list()
  rk_vec <- numeric()
  G_obs <- list()
  G_obs_refit <- list()
  lambda_list <- list()

  for (lambda1_iter in 1:length(lam1_list)) {
    for (lambda2_iter in 1:length(lam2_list)) {
      lambda_1 <- lam1_list[lambda1_iter]
      lambda_2 <- lam2_list[lambda2_iter]

      lambda_list[[r]] <- c(lambda_1, lambda_2)

      # run the sparse+low-rank estimator using CVXR
      P <- Variable(d, d, PSD = TRUE)
      L <- Variable(d, d, PSD = TRUE)
      S <- Variable(d, d, PSD = TRUE)
      R <- -log_det(t(U) %*% P %*% U) - 1 / 2 * sum(diag(P %*% Gamma)) + lambda_1 * (sum(sum(abs(S))) + lambda_2 * sum(diag(L))) # Bitrate
      objective <- Minimize(R)
      constraints <- list(P == S - L, U %*% t(U) %*% (P) %*% U %*% t(U) == P)
      prob <- Problem(objective, constraints)
      result <- solve(prob)
      G_obs[[r]]<-Theta2Gamma(result$getValue(P))
      G_obs[[r]] <- graphicalExtremes:::ensure_matrix_symmetry(Theta2Gamma(result$getValue(P), check=FALSE))

      rk_vec[r] <- length(which(eigen(result$getValue(L))$values >= 10^(-3)))
      if (rk_vec[r] == 0){
        subspace_est<-matrix(0,d,1)
      }else{
      subspace_est <- svd(result$getValue(L))$u[, 1:rk_vec[r], drop=FALSE]
      }
      output_1 <- result$getValue(S)
      output_1[which(abs(output_1) <= 10^(-3))] <- 0
      output_1[which(abs(output_1) > 10^(-3))] <- 1
      output_1 <- output_1 - diag(diag(output_1))
      off_support_est <- which(output_1 + diag(d) == 0, arr.ind = TRUE)
      graph_obs[[r]] <- graph_from_adjacency_matrix(output_1, mode = "undirected")
      
      if (refit) {
        
        A <- matrix(0,1,d^2)
        if (nrow(off_support_est)>0){
        A <- matrix(0,nrow(off_support_est),d^2)
        for (i in 1:nrow(off_support_est)){
          S <- matrix(0,d,d)
          S[off_support_est[i,1],off_support_est[i,2]]<-1
          A[i,]<- c(S)
        }
        }
        rk <- rk_vec[r]
        if (rk !=0){
        P <- Variable(d, d, PSD = TRUE)
        M <- Variable(rk, rk, PSD = TRUE)
        S <- Variable(d, d, PSD = TRUE)
        R <- -log_det(t(U) %*% P %*% U) - 1 / 2 * sum(diag(P %*% Gamma)) # Bitrate
        objective <- Minimize(R)
        constraints <- list(P == S - subspace_est %*% M %*% t(subspace_est), U %*% t(U) %*% (P) %*% U %*% t(U) == P,A%*%reshape_expr(S,c(d^2,1)) == 0)

        prob <- Problem(objective, constraints)
        result <- solve(prob)
        }
        else{
          P <- Variable(d, d, PSD = TRUE)
          S <- Variable(d, d, PSD = TRUE)
          R <- -log_det(t(U) %*% P %*% U) - 1 / 2 * sum(diag(P %*% Gamma)) # Bitrate
          objective <- Minimize(R)
          constraints <- list(P == S, U %*% t(U) %*% (P) %*% U %*% t(U) == P,A%*%reshape_expr(S,c(d^2,1)) == 0)
          
          prob <- Problem(objective, constraints)
          result <- solve(prob)
        }
        
        G_obs_refit[[r]] <- graphicalExtremes:::ensure_matrix_symmetry(Theta2Gamma(result$getValue(P), check=FALSE))
      }
      r <- r + 1
    }
  }

  if(!refit) low_rank <- 0 # added this; Armeen check!
  return(list(G_obs = G_obs, G_obs_refit = G_obs_refit, graph = graph_obs, rk = rk_vec))
}

glasso_mb2 <- function(data, samp_size, lambda){
  
  # Initialize variables
  dd <- ncol(data)
  # data_std <- scale(data)
  S_tmp <- t(data) %*% data
  data_std <- data %*% diag(diag(S_tmp)^(-1/2))
  adj.est <- array(NA, dim = c(dd, dd, length(lambda)))
  adj.ic.est <- array(NA, dim = c(dd, dd, 3))
  lambda_order <- order(lambda, decreasing = TRUE)
  lambda_dec <- sort(lambda, decreasing = TRUE)
  # there was some subtlety with the ordering since glmnet always gives back
  # in a certain order, that's why I have this here
  
  # Loop through variables
  for(i in (1:dd)){
    X <- data_std[,-i]
    Y <- data_std[,i]
    lasso_fit <- glmnet::glmnet(x = X, y = Y, family = "gaussian",
                                lambda = lambda_dec/nrow(X)*samp_size/(samp_size-1),
                                standardize=F, intercept=F)
    if(i==1){
      null.vote <- array(0, dim = c(dd, dd, length(lambda)))
      null.vote.ic <- array(0, dim = c(dd, dd, 3))
    }
    # make sure consistent with default value
    null.vote[i, -i, ] <- null.vote[i, -i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= 1e-10)
    null.vote[-i, i, ] <- null.vote[-i, i, ] +
      (abs(as.matrix(lasso_fit$beta)) <= 1e-10)
    
    dev <- (samp_size-1) * (1 - lasso_fit$dev.ratio) * lasso_fit$nulldev
    aic.idx <- which.min( dev + aic(samp_size, ncol(X) + 1) * lasso_fit$df )
    bic.idx <- which.min( dev + bic(samp_size, ncol(X) + 1) * lasso_fit$df )
    mbic.idx <- which.min( dev + mbic(samp_size, ncol(X) + 1) * lasso_fit$df )
    # aic.idx <- which.min( samp_size * log( (1 - lasso_fit$dev.ratio) * lasso_fit$nulldev )
    #                       + aic(samp_size, ncol(X) + 1) * lasso_fit$df )
    # bic.idx <- which.min( samp_size * log( (1 - lasso_fit$dev.ratio) * lasso_fit$nulldev )
    #                       + bic(samp_size, ncol(X) + 1) * lasso_fit$df )
    # mbic.idx <- which.min( samp_size * log( (1 - lasso_fit$dev.ratio) * lasso_fit$nulldev )
    #                        + mbic(samp_size, ncol(X) + 1) * lasso_fit$df)
    
    null.vote.ic[i, -i, ] <- null.vote.ic[i, -i,] +
      (abs(as.matrix(lasso_fit$beta[,c(aic.idx, bic.idx, mbic.idx)])) <= 1e-10)
    null.vote.ic[-i, i, ] <- null.vote.ic[-i, i,] +
      (abs(as.matrix(lasso_fit$beta[,c(aic.idx, bic.idx, mbic.idx)])) <= 1e-10)
  }
  adj.est[,,lambda_order] <- null.vote <= 1
  adj.ic.est <- null.vote.ic <= 1
  return(list(adj.est=adj.est, adj.ic.est = adj.ic.est))
}







sim_study_latent <- function(d = 5, 
                          n = 100,
                          p = NULL,
                          method = c("maxstable", "mpareto"), 
                          m = 2,
                          h = 2, 
                          lambda_2,
                          gen_model = c("BA", "latent"),
                          reg_method = c("eglearn", "MTP2", "eglatent"),
                          rhostring = "seq(0.01,0.15,length=20)",
                          rng = NULL){
  ## perform a simulation study to measure performance of the EMTP2 block descent algorithm.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - n: number of samples.
  ##    - p: threshold probability.
  ##    - gen_method: data generation method.
  ##    - m: the number of edges added in each setp of the Barabasi-Albert model (m=1 is a tree)
  ##    - reg_method: regression method used to estimate the extremal graphical structure.
  ##    - rholist: the list of penality parameters; must be given a string.
  ## Returns:
  ##  a tibble with F scores
  
  # check arguments
  
  rholist <-  eval(parse(text = rhostring))
  
  F1 <- numeric(length(rholist))
  F1_val <- NA
  loglik_max <- NA
  rk_max <- NA
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  
  if(gen_model=="BA"){
    BA_model <- generate_BA_model(d = d, m = m)
    G <- BA_model$G
    g <- BA_model$graph
  }
  if (gen_model == "latent") {
    latent_model <- generate_latent_model(p = d, h)
    G <- 10 * latent_model$Gamma[1:d, 1:d]
    g <- subgraph(latent_model$graph, 1:d)
  }
      
  # perform simulation
  if(method=="maxstable")  X <- rmstable(n=2*n, d=d, model="HR", par=G)
  if(method=="mpareto")  X <- rmpareto(n=2*n, d=d, model="HR", par=G)
  
  X_train <- X[1:n,]
  X_val <- X[(n+1):(2*n),]
  G_emp <- emp_vario(data = X_train, p = p)


  if(reg_method=="MTP2"){
    ptm <- proc.time()[1]
    G_emtp2 <- emtp2(G_emp, tol=1e-6,verbose = FALSE)$G_emtp2 
    time <- proc.time()[1] - ptm
    adj_emtp2 <- (abs(Gamma2Theta(G_emtp2)) >= 1e-4) 
    graph_emtp2 <- igraph::graph_from_adjacency_matrix(adj_emtp2, mode = "undirected", diag = FALSE)
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=graph_emtp2))
    loglik_val <- rep(loglik_HR(data = X_val, p=p, Gamma = G_emtp2, cens = FALSE)[1], times = length(rholist))
  }
  else if(reg_method=="eglearn"){
    ptm <- proc.time()[1]
    fit <- eglearn(data = X_train, p=p, rholist = 7*rholist, reg_method = "ns", complete_Gamma = FALSE)
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g=g, gest=fit$graph[[i]]))  
    connected_eglearn <- sapply(1:length(rholist), FUN = function(i) ifelse(is_connected(fit$graph[[i]]), 1,0))
    loglik_val <- sapply(1:length(rholist), FUN = function(i) ifelse(connected_eglearn[i]==1, loglik_HR(data = X_val, p=p, Gamma = complete_Gamma(Gamma = G_emp, graph = fit$graph[[i]]), cens = FALSE)[1], NA))
  }
  else if(reg_method=="eglatent"){
    ptm <- proc.time()[1]
    fit_latent <- eglatent(Gamma = G_emp, lam1_list = rholist, lam2_list = lambda_2, refit = TRUE)    
    time <- proc.time()[1] - ptm
    F1 <- sapply(1:length(rholist), FUN = function(i) F1_score(g = g, gest = fit_latent$graph[[i]]))
    rk <- fit_latent$rk
    loglik_val <- sapply(1:length(rholist), FUN = function(i) {loglik_HR(data = X_val, p=p, Gamma = fit_latent$G_obs_refit[[i]], cens = FALSE)[1]})
  }



  F1_max <- max(F1)
  idx_val_max <- which.max(loglik_val)
  F1_val <- F1[idx_val_max]
  loglik_max <- max(loglik_val, na.rm = TRUE)
  if(reg_method=="eglatent") rk_max <- rk[idx_val_max] 

  
  tbl <- tibble(type = paste0("time"), 
                value = time) %>% 
    bind_rows(tibble(type = paste0("F1_score", 1:length(rholist)),
                     value =  F1)) %>% 
     bind_rows(tibble(type = c("F1_max"),
                      value =  F1_max)) %>% 
    bind_rows(tibble(type = c("F1_val"),
                     value =  F1_val)) %>% 
    bind_rows(tibble(type = c("loglik_max"),
                     value =  loglik_max)) %>% 
    bind_rows(tibble(type = c("rk_max"),
                     value =  rk_max))
  
  return(tbl)
}



F1_score <- function(g, gest) {
  (2*ecount(intersection(gest, g)))/( 2*ecount(intersection(gest, g)) + 
                                        ecount(intersection(complementer(g), gest)) +
                                        ecount(intersection(g, complementer(gest))))
}


  
# traditional criteria
aic <- function(n, p) 2
bic <- function(n, p) log(n)

# modified BIC of Wang & Leng, JRSSB 2009
mbic <- function(n, p) log(n) * log(log(p))



wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Ahttps://uniart1.wixsite.com/uni-artrgs:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_fn, args = sim_fn_args[i, ]) %>% #ML: fixed the name of the last argument. It used to be fun_args
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
    dplyr::select(-rowname)
  
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


#rep_tibble_new solves an issue with package intersection

rep_tibble_new <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
}



generate_BA_model <- function(d,m){
  g <- sample_pa(n=d, m=m, zero.appeal=1,directed=FALSE)
  W <- as_adj(g, sparse=F) * matrix(runif(d^2,2, 5), nrow=d) #matrix(2 + rexp(d^2, rate = 1), nrow=d) #matrix(runif(d^2,2, 5), nrow=d) #  # 
  W[lower.tri(W)] <- t(W)[lower.tri(W)]
  O <- diag(rowSums(W)) - W
  G <- Theta2Gamma(O)
  return(list(G = G, graph = Gamma2graph(G)))
}



generate_block_model <- function(ncliques, clique_size, alphad = 1){
  kk <- clique_size
  GG <- matrix(NA, ncliques*(kk-1) + 1, ncliques*(kk-1) + 1)
  for(i in 1:ncliques){
    bigS <- rcorrmatrix(kk, alphad = alphad)
    G1 <- Sigma2Gamma(bigS, full = TRUE)
    if(i==1) GG[1:kk, 1:kk] <- G1
    else GG[(kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1), (kk + (i-2)*(kk-1)):(kk + (i-2)*(kk-1) + kk - 1)] <- G1
  }
  G <- complete_Gamma(GG)
  round(Gamma2Theta(G),2)
  sum(sum(round(Gamma2Theta(G),2) > 0)) - nrow(G)
  sum(sum(round(Gamma2Theta(G),2) < 0))
  
  return(list(G=G, graph = Gamma2graph(G,to_plot = FALSE)))
}

save_myplot <- function(plt, plt_nm,
                        width, height, 
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"), 
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}



lst_methods <- list("F1_val_eglearn" = "eglearn",
                    "F1_val_MTP2" = "emtp2",
                    "F1_val_eglatent" = "eglatent",
                    "F1_max_eglearn" = "eglearn_max",
                    "F1_max_MTP2" = "emtp2_max",
                    "F1_max_eglatent" = "eglatent_max"
)

my_palette <- list(
  "red" = "#D55E00",
  "blue" = "#0072B2", 
  "green" = "#009E73",
  "yellow" = "#E69F00",
  "pink" = "#CC79A7",
  "light_blue" = "#56B4E9",
  "grey" = "#999999",
  "background" = "#332288"
)



my_palette_methods <- list(
  c("reg_method" = lst_methods$F1_val_eglearn, "color" = my_palette$red, "fill" = "white"),
  c("reg_method" = lst_methods$F1_val_eglatent, "color" = my_palette$blue, "fill" = "white"),
  c("reg_method" = lst_methods$F1_val_MTP2, "color" = my_palette$green, "fill" = "white"),
    c("reg_method" = lst_methods$F1_max_eglearn, "color" = my_palette$red, "fill" = my_palette$red),
  c("reg_method" = lst_methods$F1_max_eglatent, "color" = my_palette$blue, "fill" = my_palette$blue),
  c("reg_method" = lst_methods$F1_max_MTP2, "color" = my_palette$green, "fill" = my_palette$green)
) %>% 
  purrr::transpose() %>% 
  as_tibble() %>% 
  unnest(cols = c(reg_method, color, fill)) 

my_col <-  my_palette_methods %>% 
  dplyr::select(reg_method, color) %>% 
  deframe()

my_fill <-  my_palette_methods %>% 
  dplyr::select(reg_method, fill) %>% 
  deframe()


refactor_methods <- function(methods, lst_methods){
  ## character_vector list with mapping -> factor
  ## refactor column with methods
  
  lst_methods <- lst_methods
  
  
  unique_methods <- unique(methods)
  
  new_levels <- names(lst_methods)
  new_labels <- lst_methods %>% unlist() %>% unname()
  
  factor(methods,
         levels = new_levels,
         labels = new_labels)
}


theme_fct <- function(font_size1=11,  font_size2=7.5){
  theme_set(theme_bw() +
              theme(
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.background = element_blank(),
                strip.background = element_rect(fill = "white"),
                plot.caption=element_text(size=font_size2, hjust=0, 
                                          margin=margin(t=15)),
                text = element_text(size = font_size1),
                axis.ticks = element_blank(),
                axis.text = element_text(size = font_size1),
                panel.grid.major = element_line(size = 0.25)
              ) 
            #+
            # theme_cowplot(font_size = 11)
  )
}

theme_fct()


create_palette_levels <- function(reg_method_n, palette_tbl){
  
  str_spl <- strsplit(reg_method_n, "__")
  
  my_tbl <- tibble(
    reg_method = purrr::map_chr(str_spl, function(el){el[1]}),
    level =  purrr::map_chr(str_spl, function(el){el[2]})
  ) %>%
    left_join(palette_tbl, by = "reg_method") %>%
    mutate(fill = if_else(level == "2", color, fill)) %>%
    mutate(reg_method_lev = paste(reg_method, level, sep = "__"))
  
  my_col <-  my_tbl %>%
    dplyr::select(reg_method_lev, color) %>%
    deframe()
  
  my_fill <-  my_tbl %>%
    dplyr::select(reg_method_lev, fill) %>%
    deframe()
  
  
  list(cols = my_col, fills = my_fill)
  
}



save_myplot <- function(plt, plt_nm,
                        width, height, 
                        width_pdf = 50, height_pdf = 50,
                        crop = TRUE, cairo = TRUE) {
  
  dir_name <- dirname(plt_nm)
  if (!file.exists(dir_name)){
    dir.create(dir_name)
  }
  
  if (cairo) {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"), 
           device = cairo_pdf, family = "Arial")
  } else {
    ggsave(plt_nm, 
           egg::set_panel_size(p = plt, 
                               width = unit(width, "in"), 
                               height = unit(height, "in")),
           width = width_pdf, height = height_pdf,
           limitsize = FALSE, units = c("in"))
  }
  
  if (crop){
    knitr::plot_crop(plt_nm)
  } 
}