


d <- 6
n <- 20
MAX_IT <- 100
nSimulations <- 20
modelSeed <- 1
seed <- 14

set.seed(modelSeed)
m <- generate_random_model(d, 'decomposable', cMax=4)

graph <- m$graph

set.seed(seed)
data <- graphicalExtremes::rmpareto(n, 'HR', d, m$Gamma)


fmpareto_HR_MLE(
    data,
    graph = graph,
    maxit = MAX_IT,
    useTheta = TRUE
)



