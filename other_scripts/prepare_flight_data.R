library(dplyr)

flights <- list()
flights$data <- mat
flights$airports <- selected_airports
flights$connections <- flights_connections %>%
  select(-ids.x, -ids.y)

usethis::use_data(flights, overwrite = TRUE)
