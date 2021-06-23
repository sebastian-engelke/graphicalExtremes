library(graphicalExtremes)
library(tidyverse)
library(glasso)

# Function definition
haversine_dist <- function(lat1, lon1, lat2, lon2){
  ## numeric (4x) -> numeric
  ## compute Haversine distance between to points on Earth

  R = 6371 # Earth's mean radius in KM
  phi1 = lat1 * pi / 180 # φ, λ in radians
  phi2 = lat2 * pi / 180
  dphi = (lat2-lat1) * pi / 180
  dlambda = (lon2-lon1) * pi / 180

  a = sin(dphi/2) * sin(dphi/2) +
    cos(phi1) * cos(phi2) *sin(dlambda/2) * sin(dlambda/2)

  c = 2 * atan2(sqrt(a), sqrt(1-a))

  d = R * c # in KM

  return(d)
}

compute_dist <- function(lat1, lon1, lat_list, lon_list){
  ## numeric (2x) numeric_vector (2x) -> numeric_vector
  ## compute Haversine distance between (lat1, lon1) and (lat_list, lon_list)

  if (length(lat_list) != length(lon_list)){
    stop("lat2_list and lon2_list must have same length.", call. = FALSE)
  }

  sapply(seq_along(lat_list), function(i){
    haversine_dist(lat1, lon1, lat_list[i], lon_list[i])
  })


}

is_within_helper <- function(lat1, lon1, lat_list, lon_list, d){
  ## numeric (2x) numeric_vector (2x) numeric -> boolean
  ## checks whether (lat1, lon1) is within d KM of some (lat_list, lon_list)
  any(compute_dist(lat1, lon1, lat_list, lon_list) < d)

}

is_within <- function(lat1_list, lon1_list, lat2_list, lon2_list, d){
  ## numeric_vector (4x) numeric -> boolean_vector
  ## produce true if the i-th entry of (lat1_list, lon1_list) is within d KM
  ## of some (lat_list, lon_list)

  if (length(lat1_list) != length(lon1_list)){
    stop("lat1_list and lon1_list must have same length.", call. = FALSE)
  }

  sapply(seq_along(lat1_list), function(i){
    is_within_helper(lat1_list[i], lon1_list[i], lat2_list, lon2_list, d)
  })
}

# haversine_dist(airports$LATITUDE[1], airports$LONGITUDE[1],
#                airports$LATITUDE[2], airports$LONGITUDE[2])
#
# is_within_helper(lat1, lon1, lat_list, lon_list, 100)
# is_within(airports2$LATITUDE, airports2$LONGITUDE,
#           large_airports$LATITUDE, large_airports$LONGITUDE,
#           100)

eglasso <- function(Gamma, rholist=c(0.1, 0.15, 0.19, 0.205), ret.graph=FALSE,
                    eps=.5, reg_method = c("mb", "glasso", "ct", "tiger")){
  r <- length(rholist)
  d <- ncol(Gamma)
  null.vote <- array(0, dim=c(d, d, length(rholist))) # votes for EXCLUDING the edge
  for(k in 1:d){
    Sk <- cov2cor(Gamma2Sigma(Gamma=Gamma, k=k))
    ###### Added regularization
    # Tk <- solve(Sk)
    # Tk.reg <- Tk + diag(rep(eps, nrow(Tk)))
    # Sk.reg <- solve(Tk.reg)
    # Ck <- Sk.reg #cov2cor(Sk.reg)

    ###### Same regularization, but does not require Sk to be invertible
    tmp <- solve(diag(ncol(Sk)) + eps*Sk)
    Ck <- tmp %*% Sk

    ###### Using "glasso" package
    approx <- (reg_method == "mb")
    if(reg_method != "glasso" && reg_method != "mb") warning(paste("Method", reg_method, "not implemended in glasso. Regular glasso was used instead."))
    invisible(capture.output(gl.tmp <- glassopath(Ck, rholist = rholist, approx=approx))) #penalize.diagonal = FALSE
    null.vote[-k,-k, ] <-  null.vote[-k,-k,] + (abs(gl.tmp$wi)==0)  ## change back to == 0 .. <=1e-4

    # ###### Using "huge" package
    # gl.tmp <- huge(Ck, lambda = rholist, method=reg_method, scr=T, verbose=F)$path
    # null.vote[-k,-k, ] <- null.vote[-k,-k,] + 1 - array(unlist(lapply(gl.tmp, function(M) as.matrix(M))), c(d-1, d-1, r))  ## change back to == 0 .. <=1e-4

  }
  adj.est <- (null.vote/(ncol(null.vote)-2)) < .49
  if(ret.graph){
    graph <- list()
    for(j in 1:r)
      graph[[j]] <- graph_from_adjacency_matrix(adj.est[,,j], mode="undirected", diag=FALSE)
    return(graph)
  }
  else
    return(adj.est)
}

# Import data
df <- read_csv("other_scripts/data/flights_data/flights.csv")
airports <- read_csv("other_scripts/data/flights_data/airports.csv")
airlines <- read_csv("other_scripts/data/flights_data/airlines.csv")

# Constansts
n_flights_large <- 1.5e5
n_flights_small <- 5e3
distance <- 400
large_airports_vec <- c("DFW", "LAX", "ORD")

# Select airports
airports2 <- df %>%
  group_by(ORIGIN_AIRPORT) %>%
  summarise(N_FLIGHTS = n()) %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  drop_na()

large_airports <- airports2 %>%
  filter(N_FLIGHTS > n_flights_large)

large_airports <- airports2 %>%
  filter(ORIGIN_AIRPORT %in% large_airports_vec)

small_airports <- airports2 %>%
  filter(is_within(LATITUDE, LONGITUDE,
                   large_airports$LATITUDE, large_airports$LONGITUDE, distance),
         N_FLIGHTS < n_flights_small)

selected_airports <- bind_rows(
  large_airports,
  small_airports
)

# Select flights
selected_flights <- df %>%
  filter(ORIGIN_AIRPORT %in% selected_airports$ORIGIN_AIRPORT,
         DESTINATION_AIRPORT %in% selected_airports$ORIGIN_AIRPORT) %>%
  filter(ORIGIN_AIRPORT %in% large_airports_vec | DESTINATION_AIRPORT %in%
           large_airports_vec)


# Edges
flights_connections <- selected_flights %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT) %>%
  distinct() %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT,
         LATITUDE.origin, LONGITUDE.origin,
         LATITUDE.dest, LONGITUDE.dest)


# Take only connections from/to main airports
selected_airports <- tibble(
  ORIGIN_AIRPORT = c(flights_connections$ORIGIN_AIRPORT,
                     flights_connections$DESTINATION_AIRPORT) %>% unique()
) %>%
  left_join(selected_airports, by = c("ORIGIN_AIRPORT")) %>%
  drop_na()

# Plot map
ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = selected_airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .1, curvature = 0) +
  theme_bw()


# Set up matrix
dat <- df %>%
  # mutate(DELAY = if_else(ARRIVAL_DELAY < 0, 0, ARRIVAL_DELAY)) %>%
  mutate(DELAY = ARRIVAL_DELAY) %>%
  group_by(DAY, MONTH, ORIGIN_AIRPORT) %>%
  summarise(SUM_DELAY = sum(DELAY)) %>%
  replace_na(list(SUM_DELAY = 0)) %>%
  filter(ORIGIN_AIRPORT %in% selected_airports$ORIGIN_AIRPORT) %>%
  ungroup()

mat <- dat %>%
  pivot_wider(id_cols = c("DAY", "MONTH"),
              names_from = ORIGIN_AIRPORT,
              values_from = SUM_DELAY) %>%
  select(-MONTH, -DAY)


# Minimum spanning tree
my_fit <- emst(data = mat, p = .8, method = "vario")
igraph::V(my_fit$tree)$name <- names(mat)

flights_connections_est <- igraph::get.edgelist(my_fit$tree) %>%
  as_tibble(.name_repair = ~ c("ORIGIN_AIRPORT", "DESTINATION_AIRPORT")) %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT,
         LATITUDE.origin, LONGITUDE.origin,
         LATITUDE.dest, LONGITUDE.dest)


# Plot estimated map
ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = selected_airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections_est,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .1, curvature = 0) +
  theme_bw()


# Glasso
Gamma <- emp_vario(mat %>% as.matrix(), p = .8)
adj_mat <- eglasso(Gamma, )
est_graph <- igraph::graph_from_adjacency_matrix(adj_mat[, , 1])
igraph::V(est_graph)$name <- names(mat)

flights_connections_est <- igraph::get.edgelist(est_graph) %>%
  as_tibble(.name_repair = ~ c("ORIGIN_AIRPORT", "DESTINATION_AIRPORT")) %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT,
         LATITUDE.origin, LONGITUDE.origin,
         LATITUDE.dest, LONGITUDE.dest)

# Plot estimated map
ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = selected_airports,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections_est,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .1, curvature = 0) +
  theme_bw()
