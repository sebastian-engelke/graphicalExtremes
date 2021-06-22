library(tidyverse)

# Import data
df <- read_csv("other_scripts/data/flights_data/flights.csv")
airports <- read_csv("other_scripts/data/flights_data/airports.csv")
airlines <- read_csv("other_scripts/data/flights_data/airlines.csv")


# Prepare data
df <- df %>%
  group_by(ORIGIN_AIRPORT) %>%
  mutate(N_FLIGHTS = n())

airline <- "US"
n_flights <- 150000

df_large_airports <- df %>%
  filter(AIRLINE == airline,
         N_FLIGHTS > n_flights) %>%
  ungroup()


# N_flights
n_flights <- df %>%
  select(ORIGIN_AIRPORT, N_FLIGHTS) %>%
  distinct()


# Edges
flights_connections <- df_large_airports %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT) %>%
  distinct() %>%
  left_join(airports, by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(airports, by = c("DESTINATION_AIRPORT" = "IATA_CODE"), suffix = c(".origin", ".dest")) %>%
  select(ORIGIN_AIRPORT, DESTINATION_AIRPORT,
         LATITUDE.origin, LONGITUDE.origin,
         LATITUDE.dest, LONGITUDE.dest)

# Nodes
airports_nodes <- bind_rows(
  flights_connections %>%
    select(ORIGIN_AIRPORT, LATITUDE.origin, LONGITUDE.origin) %>%
    rename(AIRPORT = ORIGIN_AIRPORT,
           LATITUDE = LATITUDE.origin, LONGITUDE = LONGITUDE.origin,) %>%
    distinct(),
  flights_connections %>%
    select(DESTINATION_AIRPORT, LATITUDE.dest, LONGITUDE.dest) %>%
    rename(AIRPORT = DESTINATION_AIRPORT,
           LATITUDE = LATITUDE.dest, LONGITUDE = LONGITUDE.dest) %>%
    distinct(),

) %>%
  left_join(n_flights, by = c("AIRPORT" = "ORIGIN_AIRPORT"))


# Plot map
ggplot() +
  geom_polygon(data = map_data("usa"),
               aes(x = long, y = lat, group = group), color = "grey65",
               fill = "#f9f9f9", size = 0.2) +
  geom_point(data = airports_nodes,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .5, curvature = 0) +
  theme_bw()


# Compute Gamma
dat <- df %>%
  filter(AIRLINE == airline,
         ORIGIN_AIRPORT %in% airports_nodes$AIRPORT,
         DESTINATION_AIRPORT %in% airports_nodes$AIRPORT) %>%
  mutate(DELAY = if_else(ARRIVAL_DELAY < 0, 0, ARRIVAL_DELAY)) %>%
  group_by(DAY, MONTH, ORIGIN_AIRPORT) %>%
  summarise(SUM_DELAY = sum(DELAY)) %>%
  replace_na(list(SUM_DELAY = 0)) %>%
  ungroup()

mat <- dat %>%
  pivot_wider(id_cols = c("DAY", "MONTH"),
              names_from = ORIGIN_AIRPORT,
              values_from = SUM_DELAY)


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
  geom_point(data = airports_nodes,
             aes(x = LONGITUDE, y = LATITUDE, size = N_FLIGHTS),
             alpha = 1) +
  geom_curve(data = flights_connections_est,
             aes(x = LONGITUDE.origin, xend = LONGITUDE.dest,
                 y = LATITUDE.origin, yend = LATITUDE.dest),
             alpha = .5, curvature = 0) +
  theme_bw()
