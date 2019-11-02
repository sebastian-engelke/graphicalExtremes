###################################################################
#### Code for simulation study and application of
#### Engelke & Hitz, Graphical Models for Extremes (2018, preprint)
###################################################################
setwd("other_scripts") # we will delete
library("graphicalExtremes")
library("igraph")
# library("matrixcalc")

##############################################################
#### Function definitions
##############################################################
### Estimates empirically the chi coefficient in 3 dimensions
#data: nxd data matrix
#triplets: mx3 matrix with locations out of 1:d to be evaluated
#p: probability threshold for emp_chi
#Gtrue: if supplied then the estimated chi are plotted against the once implied
#from this HR matrix
#pot: if TRUE, then pot-type estimation of chi is used
est.chi3D <- function(data, triplets, p, Gtrue=NULL, pot=FALSE, main=""){
  d <- ncol(data)
  chi <- apply(triplets, 1, function(x) emp_chi(data[,x], p=p, pot=pot))
  if(!is.null(Gtrue)){
    chi.theo = apply(triplets, 1, function(x) Gamma2chi_3D(Gamma=Gtrue[x,x]))

    par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
        mar = c(5,5,4,2) +.1)
    plot(chi, chi.theo, xlim = c(0.1,.9), ylim = c(0.1,.9), main=main,
         xlab="Fitted Model", ylab="Empirical")
    abline(0,1, xlim=c(1,2))
  }
  return(chi)
}



### Plots two chi-estimates against each other
plotChi <- function(Chi.emp,
                    Chi.theo,
                    is.con,
                    main="",
                    PDF = FALSE,
                    filename = "")
{

  if(PDF) pdf(filename, width = 7)
  par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
      mar = c(5,5,4,2) +.1)
  upper.idx <- upper.tri(Chi.emp)
  con.idx <- as.logical(upper.idx * is.con)
  plot(Chi.theo[upper.idx], Chi.emp[upper.idx], xlim = c(0.23,1), ylim = c(0.23,1),main=main,xlab="Fitted Model", ylab="Empirical")
  abline(b = 1, a = 0)
  points(Chi.theo[con.idx],Chi.emp[con.idx],col="blue")
  #legend("topleft", col=c("blue","black"), pch=19, legend=c("flow-connected pairs","flow-unconnected pairs"), bty = "n", cex = 1.5, pt.cex = 1.5)
  #print(mean((Chi.theo - Chi.emp)^2))
  if(PDF) dev.off()
}



##############################################################
#### Simulation study: first experiment
#### comparison of full MLE and MLE over cliques
##############################################################

edges = rbind(c(1,2),c(1,3),c(2,4), c(2,5))
d = max(edges)
graph0 = make_empty_graph(n = d, directed = FALSE)
for(i in 1:nrow(edges))  graph0 = add_edges(graph = graph0, edges = edges[i,])
graph0 <- set_graph_parameters(graph0)
plot.igraph(graph0)

G.vec = c(1,2,1,2)
G0 = complete_Gamma(graph = graph0, Gamma = G.vec)
p = length(G.vec)

nexp <- 2 # !!! 200
n.vec <- c(100)
method.vec <- c("full", "graph")
est.par <- array(NA, dim = c(nexp, p, length(method.vec),length(n.vec)))

######### Do not run #########
for(i in 1:nexp){
  cat("\r Simulation", i, " of ", nexp)
  for(l in 1:length(n.vec)){
    n <- n.vec[l]
    set.seed(i+10*l)
    data = rmpareto(n=n, model="HR", d=d,  par=G0)
    G.est <- emp_vario(data=data)
    est.par[i,,1,l] = (fmpareto_HR(data=data, init=G.est[edges], cens=TRUE,
                                  graph=graph0)$Gamma)[edges]
    est.par[i,,2,l] = (estGraph_HR(graph=graph0, data=data,
                                   cens=TRUE)$Gamma)[edges]
  }
}
#############################


#############################
#### Figure 4 in the paper
#############################

load("data/est.par_n100_nexp200.RData")
boxplot(matrix(rbind(t(t(est.par[,,1,1]) - G.vec),t(t(est.par[,,2,1]) - G.vec)),nrow=nexp,ncol=length(method.vec)*p),
        main = "",
        at = 1:(length(method.vec)*p) + rep((0:(p-1))/2, times=rep(length(method.vec),p)),
        las = 2,
        col = c("orange", "cyan"),
        xaxt='n', ann=FALSE,
        cex.axis=2,
        ylim=c(-1,2))
axis(side = 1, at = 1/2 + length(method.vec)*(0:(p-1)) + length(method.vec)/2 + (0:(p-1))/2,
     cex.axis = 2,
     labels = c(expression(Gamma[12]), expression(Gamma[13]), expression(Gamma[24]),expression(Gamma[25])))

load("data/est.par_n200_nexp200.RData")
boxplot(matrix(rbind(t(t(est.par[,,1,1]) - G.vec),t(t(est.par[,,2,1]) - G.vec)),nrow=nexp,ncol=length(method.vec)*p),
        main = "",
        at = 1:(length(method.vec)*p) + rep((0:(p-1))/2, times=rep(length(method.vec),p)),
        las = 2,
        col = c("orange", "cyan"),
        xaxt='n', ann=FALSE,
        cex.axis=2,
        ylim=c(-1,2))
axis(side = 1, at = 1/2 + length(method.vec)*(0:(p-1)) + length(method.vec)/2 + (0:(p-1))/2,
     cex.axis = 2,
     labels = c(expression(Gamma[12]), expression(Gamma[13]), expression(Gamma[24]),expression(Gamma[25])))


##############################################################
#### Simulation study: second experiment
#### model selection on larger graph
##############################################################


edges.tree = rbind(c(1,2),c(1,3),c(1,4),c(1,5),c(2,6),c(2,7),c(3,8),c(3,9),c(4,10),c(4,11),c(5,12),c(5,13),c(1,14),c(14,15),c(14,16))
d=max(edges.tree)
star.tree = make_empty_graph(n = d, directed = FALSE)
for(i in 1:nrow(edges.tree))  star.tree = add_edges(graph = star.tree, edges = edges.tree[i,])
added.edges.true = rbind(c(2,3),c(15,16),c(6,7))

graph.true = star.tree
for(i in 1:nrow(added.edges.true))  graph.true = add_edges(graph = graph.true, edges = added.edges.true[i,])

set.seed(9119770)
G.vec = runif(ecount(graph.true), min = 0.5, max = 1)
Gamma = complete_Gamma(graph = graph.true, Gamma = G.vec) ## choose seed such that all cliques have valid Gamma matrix

Gamma2graph(Gamma = Gamma)

nexp <- 1 # !!! 100
n <- c(100)
est.AIC <- array(NA, dim = c(nexp, 2*d))
nb.chosen.edges = matrix(0,d,d)

#### Do not run ####
for(i in 1:nexp){
  cat("\r Simulation", i, " of ", nexp)
  set.seed(i)
  data = rmpareto(n=n, model="HR", d=d, par=Gamma)
  tree.tmp = mst_HR(data = data, cens = TRUE)
  sel.edges = select_edges(graph=tree.tmp)
  fit.tmp = estGraph_HR(graph=tree.tmp, data=data, cens=TRUE, edges_to_add =sel.edges)
  est.AIC[i,1:length(fit.tmp$AIC)] = fit.tmp$AIC
  chosen.edges.tmp = (rbind(ends(tree.tmp,E(tree.tmp)),fit.tmp$added.edges))[1:vcount(graph.true),]
  nb.chosen.edges[chosen.edges.tmp] = nb.chosen.edges[chosen.edges.tmp] + 1
  nb.chosen.edges[chosen.edges.tmp[,c(2,1)]] = nb.chosen.edges[chosen.edges.tmp[,c(2,1)]] + 1
}
nb.chosen.edges = round(nb.chosen.edges / max(nb.chosen.edges),2) * 100
###################

load(file="data/est.AIC.RData")
load(file="data/nb.chosen.edges.RData")

#############################
#### Data for Figure 5
#### in the paper
#############################

nb.chosen.edges

AIC.freq = hist(apply(est.AIC,MARGIN=1,FUN=which.min),breaks=1:8)$counts / 100
AIC.freq


##############################################################
#### Application: Danube river data
##############################################################


##### Data from Asadi et al. (2015, Annals of Applied Statistics)
load("data/DataEvents.RData")
load("data/FlowCon.Rdata")  ### flow connection matrix != adjacency matrix
load("data/GfitM4.Rdata")    ### results from spatial model in Asadi et al. (2015)
load("data/coordinates_river.RData")

d <- ncol(DataEvents)

edges = rbind(c(12,11),c(11,10),c(10,9),c(9,8),c(8,7),c(7,6),c(6,5),c(5,4),c(4,3),c(3,2),c(2,1),c(22,21),c(21,20),c(20,7),c(19,18),
              c(18,17),c(17,16),c(16,15),c(15,14),c(14,2),c(29,28),c(28,31),c(31,30),c(30,13),c(13,1),c(27,26),c(26,25),c(25,4),c(24,23),c(23,4))

FlowGraphDirected = make_empty_graph(n = d, directed = TRUE)
for(i in 1:nrow(edges)) FlowGraphDirected = add_edges(graph = FlowGraphDirected, edges = edges[i,])
FlowGraph = set_graph_parameters(as.undirected(FlowGraphDirected))


##################################
#### Figure 6 (left) in the paper
##################################

plot(FlowGraph,  layout = coordinates_river, vertex.color="lightblue", vertex.size=10, edge.width=2)


X <- data2mpareto(data=DataEvents, p=0.9)
mstFit = mst_HR(data = X, cens = TRUE)

##################################
#### Figure 9 (left) in the paper
##################################

plot(mstFit,  layout = coordinates_river)



##################################
#### Figure 9 (right) in the paper
#### Gaussian MST
##################################

graph.full <- make_full_graph(d)
logDataEvents = log(DataEvents)
Gauss_mst = set_graph_parameters(igraph::mst(graph=graph.full, weights = log(1- (cor(logDataEvents)[ends(graph.full,E(graph.full))])^2), algorithm = "prim"))

plot(Gauss_mst, layout = coordinates_river)


##################################
#### Figure 10 in the paper
##################################

thr.vector = c(.4,.5,.6,.7,.75,.775,.8,.825,.85,.875,.9,.925,.95)

#### Do not run ####
MSTAll = list()
for(t in 1:length(thr.vector)){
  print(t)
  Xtmp <- data2mpareto(data=DataEvents, p=thr.vector[t])
  MSTAll[[t]] = mst_HR(data = Xtmp, cens = TRUE)
}
###################

load("data/MSTAll.RData")

A90 = as.matrix(as_adjacency_matrix(mstFit))
sameEdges = numeric()
for(t in 1:length(thr.vector)){
  AA = as.matrix(as_adjacency_matrix(MSTAll[[t]]))
  sameEdges[t] = sum(AA[upper.tri(AA)] & A90[upper.tri(A90)])
}

AGauss = as.matrix(as_adjacency_matrix(Gauss_mst))
sameEdgesGauss = sum(AGauss[upper.tri(AGauss)] & A90[upper.tri(A90)])

par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, mar = c(5,5,4,2) +.1)
plot(thr.vector, sameEdges, type="b", ylab="Identical edges", main="",
        xlab="Threshold quantile", pch=19, lwd=2, col=c('blue'), ylim = c(20,30), xlim=c(.4,.95)
)
abline(sameEdgesGauss, 0, lty=2, col="orange", lwd=2)




########################################
#### Figure 6 (right) and 7 in the paper
########################################

### Do not run ####
sel.edges = select_edges(FlowGraph)
Mfit_flow = estGraph_HR(graph=FlowGraph, data=X, cens=TRUE, edges_to_add=sel.edges)

sel.edges = select_edges(mstFit)
Mfit_mst = estGraph_HR(graph=mstFit, data=X, cens=TRUE, edges_to_add=sel.edges)
###################

load("data/Mfit_flow.Rdata")
load("data/Mfit_mst.Rdata")

L = length(Mfit_flow$AIC)
nX = nrow(X)
p.vec = ecount(FlowGraph)+0:(L-1)

AIC.min.idx = which.min(Mfit_flow$AIC)
Mfit_flow$added.edges[1:AIC.min.idx,]

plot(Mfit_mst$graph[[1]],  layout = coordinates_river)

plot(Mfit_flow$graph[[1]],  layout = coordinates_river)

plot(Mfit_flow$graph[[AIC.min.idx]], layout = coordinates_river)

results = matrix(NA, nrow=4, ncol=3)
colnames(results) = c("twice neg logLH", "nb par", "AIC")
rownames(results) = c("tree", "best graph", "largest graph", "Asadi et al.")

results[,1] = -2*c(logLH_HR(data=X,Gamma=Mfit_flow$Gamma[[1]],cens=TRUE),
                   logLH_HR(data=X,Gamma=Mfit_flow$Gamma[[AIC.min.idx]],cens=TRUE),
                   logLH_HR(data=X,Gamma=Mfit_flow$Gamma[[L]],cens=TRUE),
                   logLH_HR(data=X,Gamma=GfitM4, cens=TRUE))

results[,2] = c(p.vec[1], p.vec[AIC.min.idx], p.vec[L], 6)
results[,3] = 2*results[,2] + results[,1]
round(results,2)

par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
    mar = c(5,5,4,2) +.1)
matplot(p.vec, cbind(Mfit_flow$AIC, c(Mfit_mst$AIC,NA)), type="b", ylab="AIC", main="",
        xlab="number of edges", pch=19, lwd=2, col=c('black', 'blue')
)
abline(results[4,3],0, lty=2, col="orange", lwd=2)


chi.emp = emp_chi_mat(data=DataEvents, p=.9, pot=TRUE)
plotChi(Chi.emp = chi.emp, Chi.theo = Gamma2chi(Mfit_flow$Gamma[[AIC.min.idx]]),
        main="Hüsler-Reiss graphical model",
        PDF = FALSE,
        is.con = FlowCon)
plotChi(Chi.emp = chi.emp, Chi.theo = Gamma2chi(GfitM4),
        main="Hüsler-Reiss model in Asadi et al. (2015)",
        PDF = FALSE,
        is.con = FlowCon)


##############################
#### Figure 11 in the paper
##############################

set.seed(222)
triplets <- as.matrix(expand.grid(1:d,1:d, 1:d))[sample(1:d^3, size = 400, replace = FALSE),]
triplets = triplets[which(apply(triplets,1, function(x) length(unique(x)))==3),,drop=FALSE]

chi3D.spatial = est.chi3D(data=DataEvents, triplets=triplets, p=.9, Gtrue=GfitM4,
                          pot=TRUE, main="Hüsler-Reiss model in Asadi et al. (2015)")
chi3D.graph = est.chi3D(data=DataEvents, triplets=triplets, p=.9,
                        Gtrue=Mfit_flow$Gamma[[AIC.min.idx]], pot=TRUE, main="Hüsler-Reiss graphical model")





