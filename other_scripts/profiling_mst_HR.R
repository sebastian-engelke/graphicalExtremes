library(profvis)
devtools::load_all()
library(igraph)
library(matrixcalc)
library(mvtnorm)
library(tidyverse)
source("other_scripts/Tree_functions.R")


profvis(
  sim_study(d = 20, n = 1e2, p = .9, model = "HR", method = "maxstable",
          ML = TRUE, cens = TRUE, edge.max = 1, noise = "tree"))
