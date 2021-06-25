library(graphicalExtremes)
library(tidyverse)

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

  if (length(lat1_list) == 0) {
    return(FALSE)
  }

  sapply(seq_along(lat1_list), function(i){
    is_within_helper(lat1_list[i], lon1_list[i], lat2_list, lon2_list, d)
  })
}

select_airports <- function(df,
                            airline, state = NULL, n_large_airports = NULL,
                            large_airports_vec,
                            n_flights_small, distance){
  ## tibble character (2x) numeric character_vector numeric numeric -> tibble
  ## select list of airports
  airports2 <- df %>%
    filter(AIRLINE %in% airline) %>%
    group_by(ORIGIN_AIRPORT) %>%
    summarise(N_FLIGHTS = n()) %>%
    left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
    drop_na()

    if (is.null(state)) {
      if (is.null(n_large_airports)) {
        large_airports <- airports2 %>%
          filter(ORIGIN_AIRPORT %in% large_airports_vec)
      } else {
        large_airports <- airports2 %>%
          slice_max(N_FLIGHTS, n = n_large_airports)
      }


      small_airports <- airports2 %>%
        filter(is_within(LATITUDE, LONGITUDE,
                         large_airports$LATITUDE, large_airports$LONGITUDE,
                         distance)) %>%
      filter(N_FLIGHTS < n_flights_small)

      selected_airports <- bind_rows(
        large_airports,
        small_airports
      ) %>%
        distinct()

    } else {
      selected_airports <- airports2 %>%
        filter(STATE %in% state)
    }

  return(selected_airports)
}

# Import data
df <- read_csv("other_scripts/data/flights_data/flights.csv")
airports <- read_csv("other_scripts/data/flights_data/airports.csv")
airlines <- read_csv("other_scripts/data/flights_data/airlines.csv")

# Constansts
n_large_airports <- NULL
large_airports_vec <- c("JFK", "DFW") #c("LAX", "SAN", "SFO") #c("DFW", "LAX", "ORD")
n_flights_small <- 3e4
distance <- 300
airline <- "WN"
state <- c("TX", "CO", "AZ", "GA", "IL", "NM", "NV", "OK")
state <- c("CA", "NV", "AZ", "UT", "TX")


# Select airports
selected_airports <- select_airports(df, airline,
                                     state, n_large_airports, large_airports_vec,
                                     n_flights_small, distance)

# Select flights
selected_flights <- df %>%
  filter(AIRLINE %in% airline) %>%
  filter(ORIGIN_AIRPORT %in% selected_airports$ORIGIN_AIRPORT,
         DESTINATION_AIRPORT %in% selected_airports$ORIGIN_AIRPORT)


# Edges
flights_connections <- selected_flights %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT) %>%
  group_by(ORIGIN_AIRPORT, DESTINATION_AIRPORT) %>%
  summarise(N_FLIGHTS = n()) %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT, N_FLIGHTS,
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
                 y = LATITUDE.origin, yend = LATITUDE.dest,
                 size = N_FLIGHTS),
             alpha = .2, curvature = 0) +
  theme_bw()

# Set up matrix
dat <- selected_flights %>%
  mutate(DELAY = ARRIVAL_DELAY) %>%
  group_by(DAY, MONTH, ORIGIN_AIRPORT) %>%
  summarise(SUM_DELAY = sum(DELAY, na.rm = TRUE)) %>%
  filter(ORIGIN_AIRPORT %in% selected_airports$ORIGIN_AIRPORT) %>%
  ungroup() %>%
  rename(day = DAY)

mat <- dat %>%
  pivot_wider(id_cols = c("day", "MONTH"),
              names_from = "ORIGIN_AIRPORT",
              values_from = "SUM_DELAY") %>%
  select(-MONTH, -day)

# pairs(mat)
# plot(mat$AMA, mat$LGA)
# plot(mat$DAL, type = "line")


# Minimum spanning tree
p <- .7
my_fit_tree <- emst(data = mat, p = p, method = "vario")
igraph::V(my_fit_tree$graph)$name <- names(mat)

flights_connections_est <- igraph::get.edgelist(my_fit_tree$graph) %>%
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
             alpha = .4, curvature = 0) +
  theme_bw()


# Glasso
Gamma <- emp_vario(mat %>% as.matrix(), p = p)
rholist = seq(1e-4, 0.09, length.out = 10)
my_fit <- eglasso(Gamma, rholist = rholist, complete_Gamma = TRUE)
est_graph <- my_fit$graph[[6]]
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
             alpha = .4, curvature = 0) +
  theme_bw()


# BIC
flights_loglik <- sapply(1:length(rholist), FUN = function(j)
  loglik_HR(data=mat, p=p, Gamma = my_fit$Gamma[[j]],
            graph =  my_fit$graph[[j]]))


flights_loglik_tree <- loglik_HR(data=mat, p=p,
                                 Gamma = my_fit_tree$Gamma,
                                 graph = my_fit_tree$graph)

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
                        labels = sapply(my_fit$graph, igraph::gsize),
                        name="Number of edges")
  )

# introduce dataset + how to plot flight connections
# emst
ggplot() +
  geom_point(aes(x = c(Gamma2chi(my_fit_tree$Gamma)),
                 y = c(emp_chi(mat %>% as.matrix(), p = p)))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")

# eglasso
best_Gamma <- my_fit$Gamma[[which.min(flights_loglik[3,])]]

ggplot() +
  geom_point(aes(x = c(Gamma2chi(best_Gamma)),
                 y = c(emp_chi(mat %>% as.matrix(), p = p)))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")

# fix flights_connections to remove duplicates
flight_graph <- flights_connections %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT) %>%
  as.matrix() %>%
  igraph::graph_from_edgelist(directed = FALSE) %>%
  igraph::simplify()


plot(flight_graph)
# fmpareto_graph_HR: check completed Gamma agrees with given graph
model_fit <- fmpareto_graph_HR(data = mat %>% as.matrix(),
                               graph = flight_graph, p = p, method = "vario")


loglik_HR(data = mat %>% as.matrix(), p = p, graph = flight_graph,
          Gamma = model_fit$Gamma)

igraph::gsize(flight_graph)

ggplot() +
  geom_point(aes(x = c(Gamma2chi(model_fit$Gamma)),
                 y = c(emp_chi(mat %>% as.matrix(), p = p)))) +
  geom_abline(slope = 1, intercept = 0) +
  xlab("Fitted") +
  ylab("Empirical")


# Scatter
plot(mat[, c(4, 9)])

# Emp chi
emp_chi(mat %>% as.matrix(), p = .5) %>% hist(xlim=c(0, 1))
emp_chi(mat %>% as.matrix(), p = .7) %>% hist(xlim=c(0, 1))
emp_chi(mat %>% as.matrix(), p = .8) %>% hist(xlim=c(0, 1))

evd::chiplot(mat[, c(1, 3)])
