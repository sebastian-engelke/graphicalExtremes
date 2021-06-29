knitr::opts_chunk$set(
collapse = TRUE,
comment = ">"
)

library(graphicalExtremes)
library(dplyr)
library(ggplot2)
theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank(),
                  strip.background = element_rect(fill = "white"),
                  plot.caption=element_text(size=7.5, hjust=0,
                                            margin=margin(t=15)),
                  text = element_text(size = 11),
                  axis.ticks = element_blank(),
                  panel.grid.major = element_line(size = 0.25)))

# Load the dataset
data <- flights$data
airports <- flights$airports
connections <- flights$connections

mat <- data %>%
  select(-DATE) %>%
  as.matrix()

# Plot flight connection
ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = connections,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest,
                 size = N_FLIGHTS),
             alpha = .2, curvature = 0)


flight_graph <- connections %>%
  select(ID.origin, ID.dest) %>%
  as.matrix() %>%
  igraph::graph_from_edgelist(directed = FALSE)


# Function definition
plot_connections <- function(graph, airports) {
  ## igraph tibble tibble -> ggplot
  ## plots the given `graph` on the US map

  # name the graph nodes
  igraph::V(graph)$name <- airports$IATA_CODE[airports$ID == igraph::V(graph)]

  # write flight connections
  flights_connections_est <- igraph::get.edgelist(graph) %>%
  as_tibble(.name_repair = ~ c("ORIGIN_AIRPORT", "DESTINATION_AIRPORT")) %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT,
         LATITUDE.origin, LONGITUDE.origin,
         LATITUDE.dest, LONGITUDE.dest)

  # plot connections
  ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections_est,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .4, curvature = 0)


}

# Plot an igraph object
plot_connections(flight_graph, airports)

p <- .7
flights_emst_fit <- emst(data = mat, p = p, method = "vario")


plot_connections(flights_emst_fit$graph, airports)


flights_loglik_tree <- loglik_HR(data=mat, p=p,
                                 Gamma = flights_emst_fit$Gamma,
                                 graph = flights_emst_fit$graph)

paste0("Tree BIC = ", flights_loglik_tree[3] %>% round(2))


emp_chi_mat <- emp_chi(mat, p = p)

ggplot() +
  geom_point(aes(x = c(Gamma2chi(flights_emst_fit$Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")

model_fit <- fmpareto_graph_HR(data = mat,
                               graph = flight_graph, p = p, method = "vario")

flights_loglik_graph <- loglik_HR(data = mat,
                                  p = p, graph = flight_graph,
                                  Gamma = model_fit$Gamma)
paste0("BIC = ", flights_loglik_graph[3] %>% round(2))


ggplot() +
  geom_point(aes(x = c(Gamma2chi(model_fit$Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")


Gamma <- emp_vario(mat, p = p)
rholist = seq(1e-4, 0.10, length.out = 10)
flights_eglasso_fit <- eglasso(Gamma, rholist = rholist, complete_Gamma = TRUE)


plot_connections(flights_eglasso_fit$graph[[10]], airports)


flights_loglik <- sapply(seq_along(rholist), FUN = function(j) {
 loglik_HR(data=mat, p=p,
           Gamma = flights_eglasso_fit$Gamma[[j]],
           graph = flights_eglasso_fit$graph[[j]] )
})

ggplot(mapping = aes(x = rholist, y = flights_loglik[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_fit$graph,
                                        igraph::gsize),
                        name="Number of edges")
  )


best_Gamma <- flights_eglasso_fit$Gamma[[which.min(flights_loglik[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")


Gamma_vario_k_1 <- emp_vario(data = mat, k = 1, p = p)
flights_eglasso_k_1 <- eglasso(Gamma_vario_k_1, rholist,
                               complete_Gamma = TRUE)

Gamma_ml <- graphicalExtremes:::ml_weight_matrix(data = mat, p = p)
flights_eglasso_ml <- eglasso(Gamma_ml$est_gamma, rholist,
                               complete_Gamma = TRUE)


flights_loglik_k_1 <- sapply(seq_along(rholist), FUN = function(j) {
 loglik_HR(data=mat, p=p,
           Gamma = flights_eglasso_k_1$Gamma[[j]],
           graph = flights_eglasso_k_1$graph[[j]] )
})

not_na <- !is.na(flights_eglasso_ml$Gamma)

flights_loglik_ml <- sapply(seq_along(rholist[not_na]), FUN = function(j) {
 loglik_HR(data=mat, p=p,
           Gamma = flights_eglasso_ml$Gamma[[j]],
           graph = flights_eglasso_ml$graph[[j]] )
})


ggplot(mapping = aes(x = rholist, y = flights_loglik[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_fit$graph,
                                        igraph::gsize),
                        name="Number of edges")
  ) +
  ggtitle("Empirical variogram")

ggplot(mapping = aes(x = rholist, y = flights_loglik_k_1[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_k_1$graph,
                                        igraph::gsize),
                        name="Number of edges")
  ) +
  ggtitle("Empirical variogram with k = 1")



ggplot(mapping = aes(x = rholist[not_na], y = flights_loglik_ml[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_ml$graph,
                                        igraph::gsize),
                        name="Number of edges")
  ) +
  ggtitle("ML variogram")




best_Gamma <- flights_eglasso_fit$Gamma[[which.min(flights_loglik[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical") +
  ggtitle("Empirical variogram")

best_Gamma <- flights_eglasso_k_1$Gamma[[which.min(flights_loglik_k_1[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical") +
  ggtitle("Empirical variogram with k = 1")

best_Gamma <- flights_eglasso_ml$Gamma[[which.min(flights_loglik_ml[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical") +
  ggtitle("ML variogram")
library(graphicalExtremes)
library(dplyr)
library(ggplot2)
theme_set(theme_bw() +
            theme(plot.background = element_blank(),
                  legend.background = element_blank(),
                  strip.background = element_rect(fill = "white"),
                  plot.caption=element_text(size=7.5, hjust=0,
                                            margin=margin(t=15)),
                  text = element_text(size = 11),
                  axis.ticks = element_blank(),
                  panel.grid.major = element_line(size = 0.25)))

# Load the dataset
data <- flights$data
airports <- flights$airports
connections <- flights$connections

mat <- data %>%
  select(-DATE) %>%
  as.matrix()

# Plot flight connection
ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = connections,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest,
                 size = N_FLIGHTS),
             alpha = .2, curvature = 0)


flight_graph <- connections %>%
  select(ID.origin, ID.dest) %>%
  as.matrix() %>%
  igraph::graph_from_edgelist(directed = FALSE)


# Function definition
plot_connections <- function(graph, airports) {
  ## igraph tibble tibble -> ggplot
  ## plots the given `graph` on the US map

  # name the graph nodes
  igraph::V(graph)$name <- airports$IATA_CODE[airports$ID == igraph::V(graph)]

  # write flight connections
  flights_connections_est <- igraph::get.edgelist(graph) %>%
  as_tibble(.name_repair = ~ c("ORIGIN_AIRPORT", "DESTINATION_AIRPORT")) %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT,
         LATITUDE.origin, LONGITUDE.origin,
         LATITUDE.dest, LONGITUDE.dest)

  # plot connections
  ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections_est,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .4, curvature = 0)


}

# Plot an igraph object
plot_connections(flight_graph, airports)

p <- .7
flights_emst_fit <- emst(data = mat, p = p, method = "vario")


plot_connections(flights_emst_fit$graph, airports)


flights_loglik_tree <- loglik_HR(data=mat, p=p,
                                 Gamma = flights_emst_fit$Gamma,
                                 graph = flights_emst_fit$graph)

paste0("Tree BIC = ", flights_loglik_tree[3] %>% round(2))


emp_chi_mat <- emp_chi(mat, p = p)

ggplot() +
  geom_point(aes(x = c(Gamma2chi(flights_emst_fit$Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")

model_fit <- fmpareto_graph_HR(data = mat,
                               graph = flight_graph, p = p, method = "vario")

flights_loglik_graph <- loglik_HR(data = mat,
                                  p = p, graph = flight_graph,
                                  Gamma = model_fit$Gamma)
paste0("BIC = ", flights_loglik_graph[3] %>% round(2))


ggplot() +
  geom_point(aes(x = c(Gamma2chi(model_fit$Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")


Gamma <- emp_vario(mat, p = p)
rholist = seq(1e-4, 0.10, length.out = 10)
flights_eglasso_fit <- eglasso(Gamma, rholist = rholist, complete_Gamma = TRUE)


plot_connections(flights_eglasso_fit$graph[[10]], airports)


flights_loglik <- sapply(seq_along(rholist), FUN = function(j) {
 loglik_HR(data=mat, p=p,
           Gamma = flights_eglasso_fit$Gamma[[j]],
           graph = flights_eglasso_fit$graph[[j]] )
})

ggplot(mapping = aes(x = rholist, y = flights_loglik[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_fit$graph,
                                        igraph::gsize),
                        name="Number of edges")
  )


best_Gamma <- flights_eglasso_fit$Gamma[[which.min(flights_loglik[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")


Gamma_vario_k_1 <- emp_vario(data = mat, k = 1, p = p)
flights_eglasso_k_1 <- eglasso(Gamma_vario_k_1, rholist,
                               complete_Gamma = TRUE)

Gamma_ml <- graphicalExtremes:::ml_weight_matrix(data = mat, p = p)
flights_eglasso_ml <- eglasso(Gamma_ml$est_gamma, rholist,
                               complete_Gamma = TRUE)


flights_loglik_k_1 <- sapply(seq_along(rholist), FUN = function(j) {
 loglik_HR(data=mat, p=p,
           Gamma = flights_eglasso_k_1$Gamma[[j]],
           graph = flights_eglasso_k_1$graph[[j]] )
})

not_na <- !is.na(flights_eglasso_ml$Gamma)

flights_loglik_ml <- sapply(seq_along(rholist[not_na]), FUN = function(j) {
 loglik_HR(data=mat, p=p,
           Gamma = flights_eglasso_ml$Gamma[[j]],
           graph = flights_eglasso_ml$graph[[j]] )
})


ggplot(mapping = aes(x = rholist, y = flights_loglik[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_fit$graph,
                                        igraph::gsize),
                        name="Number of edges")
  ) +
  ggtitle("Empirical variogram")

ggplot(mapping = aes(x = rholist, y = flights_loglik_k_1[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_k_1$graph,
                                        igraph::gsize),
                        name="Number of edges")
  ) +
  ggtitle("Empirical variogram with k = 1")



ggplot(mapping = aes(x = rholist[not_na], y = flights_loglik_ml[3, ])) +
  geom_line() +
  geom_point(shape = 21, size = 3, stroke = 1, fill = "white") +
  geom_hline(aes(yintercept = flights_loglik_tree[3]), lty = "dashed") +
  xlab("rho") +
  ylab("BIC") +
  scale_x_continuous(
    breaks = rholist,
    labels = round(rholist, 3),
    sec.axis = sec_axis(trans=~., breaks = rholist,
                        labels = sapply(flights_eglasso_ml$graph,
                                        igraph::gsize),
                        name="Number of edges")
  ) +
  ggtitle("ML variogram")




best_Gamma <- flights_eglasso_fit$Gamma[[which.min(flights_loglik[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical") +
  ggtitle("Empirical variogram")

best_Gamma <- flights_eglasso_k_1$Gamma[[which.min(flights_loglik_k_1[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical") +
  ggtitle("Empirical variogram with k = 1")

best_Gamma <- flights_eglasso_ml$Gamma[[which.min(flights_loglik_ml[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi_mat))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical") +
  ggtitle("ML variogram")
