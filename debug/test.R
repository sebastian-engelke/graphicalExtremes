
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

# g <- generate_random_connected_graph(d, p = 3/(d+1))
# g <- generate_random_connected_graph(d)
# g <- igraph::make_ring(d)
g <- generate_random_tree(d)
G0 <- ensure_symmetry(generate_random_graphical_Gamma(g))
# G0 <- ensure_symmetry(generate_random_Gamma(d), Inf)
par <- G0


data <- rmpareto(n, 'HR', d, par)

G_emp <- emp_vario(data)

init <- upper.tri.val(G_emp)

tic()
cat('MLE Gamma (fix)...\n')
par2f <- fmpareto_HR_MLE(
    data,
    init = init,
    fixParams = 1,
    graph = g,
    useTheta = FALSE,
    cens = FALSE,
    p = NULL
)
toc()

tic()
cat('MLE Gamma...\n')
par2 <- fmpareto_HR_MLE(
    data,
    init = init,
    graph = g,
    useTheta = FALSE,
    cens = FALSE,
    p = NULL
)
toc()


# tic()
# cat('MLE Theta...\n')
# par3 <- fmpareto_HR_MLE(data, graph = g, useTheta = TRUE, cens = TRUE, p = 0.9)
# toc()

