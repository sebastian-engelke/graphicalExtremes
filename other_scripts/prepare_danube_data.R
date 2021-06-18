library(dplyr)

info <- danube$StsInfo %>%
  as_tibble() %>%
  mutate(Lat_Center = danube$StsCenter[, 2],
         Long_Center = danube$StsCenter[, 1],
         Alt = danube$StsAlt,
         Area = danube$StsArea,
         Chos = danube$StsChos,
         Density = danube$StsDensity,
         Slope = danube$StsSlope)

data <- danube$DataEvents
flow_edges <- danube$flow_connections
coords_to_plot <- danube$coords_to_plot

danube2 <- list()
danube2$data <- data
danube2$info <- info
danube2$flow_edges <- flow_edges
danube2$coords_to_plot <- coords_to_plot

danube <- danube2

str(danube)

usethis::use_data(danube, overwrite = TRUE)
