library(dplyr)


flights <- list()
flights$data <- mat
flights$airports <- selected_airports %>%
  rename(IATA_CODE = ORIGIN_AIRPORT) %>%
  arrange(IATA_CODE) %>%
  mutate(ID = seq_len(nrow(selected_airports))) %>%
  select(ID, everything())
flights$connections <- flights_connections %>%
  select(-ids.x, -ids.y) %>%
  ungroup() %>%
  left_join(flights$airports %>% select(ID, IATA_CODE),
            by = c("ORIGIN_AIRPORT" = "IATA_CODE")) %>%
  left_join(flights$airports %>% select(ID, IATA_CODE),
            by = c("DESTINATION_AIRPORT" = "IATA_CODE"),
            suffix = c(".origin", ".dest")) %>%
  select(ID.origin, ID.dest, everything())

usethis::use_data(flights, overwrite = TRUE)
