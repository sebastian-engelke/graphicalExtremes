
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}
library(igraph)
library(tictoc)

newSeed <- floor(2^20 * runif(1))
newSeed <- 494411
cat('Seed:', newSeed, '\n')
set.seed(newSeed)


d <- 5
n <- 100
graphType <- c('general', 'decomposable', 'tree')[1]

m <- generate_random_model(d, graphType)

G0 <- m$Gamma
graph <- m$graph

data <- rmpareto(n, 'HR', d, G0)

ret <- fmpareto_graph_HR(
    data,
    graph,
    method = 'ML',
    handleCliques = 'average'
)


