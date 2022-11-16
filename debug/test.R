
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

# Check Theta-zeros in non-edge entries:
max(abs(getNonEdgeEntries(T1, graph)))

# Check unchanged Gamma in edge entries:
max(abs(getEdgeEntries(G1 - vario_emp, graph)))

# Check that G1 is valid:
is_sym_cnd(G1)



