
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
  devtools::load_all('.')
}
library(igraph)
library(tictoc)

newSeed <- floor(2^20 * runif(1))
newSeed <- 875916
cat('Seed:', newSeed, '\n')
set.seed(newSeed)


d <- 10
n <- 100
graphType <- c('general', 'decomposable', 'tree')[1]

m <- generate_random_model(d, graphType)

G0 <- m$Gamma
graph <- m$graph

data <- rmpareto(n, 'HR', d, G0)

Gammas <- list()

solvers <- list(
  function(data){
    fmpareto_HR_MLE(
      data,
      graph = graph,
      useTheta = FALSE
    )
  },
  function(data){
    fmpareto_HR_MLE(
      data,
      graph = graph,
      useTheta = TRUE
    )
  },
  function(data){
    fmpareto_HR_MLE(
      data,
      graph = graph,
      cens = TRUE
    )
  },
  function(data){
    fmpareto_graph_HR(
      data,
      graph,
      method = 'ML',
      handleCliques = 'full'
    )
  },
  function(data){
    fmpareto_graph_HR(
      data,
      graph,
      method = 'ML',
      handleCliques = 'average'
    )
  },
  function(data){
    fmpareto_graph_HR(
      data,
      graph,
      method = 'ML',
      handleCliques = 'sequential'
    )
  },
  function(data){
    fmpareto_graph_HR(
      data,
      graph,
      method = 'vario',
      handleCliques = 'full'
    )
  },
  function(data){
    fmpareto_graph_HR(
      data,
      graph,
      method = 'vario',
      handleCliques = 'average'
    )
  }
)

nSolvers <- length(solvers)
results <- replicate(nSolvers, NULL)

for(i in seq_along(solvers)){
  solver <- solvers[[i]]
  fBody <- paste0(format(body(solver)), collapse = '\n')
  cat(i, '/', nSolvers, ':\n', fBody, '\n...\n', sep='')
  ret <- list()
  t0 <- as.numeric(Sys.time())
  tmp <- tryCatch(
    list(ret0 = solver(data)),
    error = function(e){
      list(err = e)
    }
  )
  ret$ret0 <- tmp$ret0
  ret$err <- tmp$err
  t1 <- as.numeric(Sys.time())
  if(is.matrix(ret$ret0)){
    ret$Gamma <- ret$ret0
  } else if(is.list(ret$ret0)){
    ret$Gamma <- ret$ret0$Gamma
  }
  ret$t <- (t1 - t0)
  results[[i]] <- ret
}

for(i in seq_along(results)){
  Gamma <- results[[i]]$Gamma
  if(is.matrix(Gamma)){
    Theta <- Gamma2Theta(Gamma)
    results[[i]]$Theta <- Theta
    results[[i]]$maxGammaDiff <- max(abs(getEdgeEntries(Gamma - G0, graph)))
    results[[i]]$meanGammaDiff <- mean(abs(getEdgeEntries(Gamma - G0, graph)))
    results[[i]]$maxThetaDiff <- max(abs(getNonEdgeEntries(Theta, graph)))
  }
}

print(sapply(results, '[[', 'maxGammaDiff'))
print(sapply(results, '[[', 'meanGammaDiff'))
print(sapply(results, '[[', 'maxThetaDiff'))


