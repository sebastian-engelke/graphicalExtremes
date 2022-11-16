
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}
# library(igraph)
# library(tictoc)

load('debug/graph.Rdata')
load('debug/variogram.Rdata')

# loads:
vario_emp
vario_emp <- ensure_symmetry(vario_emp)
graph


G1 <- complete_Gamma_general_sc(vario_emp, graph, N=100000)

T1 <- Gamma2Theta(G1)

max(abs(getNonEdgeEntries(T1, graph)))

# g1 <- Gamma2graph(G1, tol = 1e-3)


