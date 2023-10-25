
# file names
RAW_DIR <- 'data-raw'
OUT_FILE <- 'data/flights.rda'

FILE_DELAYS <- file.path(RAW_DIR, 'delays.rds')
FILE_FLICHT_COUNTS <- file.path(RAW_DIR, 'nFlightsYearly.rds')
FILE_CSV_AIRPORTS <- file.path(RAW_DIR, 'airportsUsed.csv')


# # specify years to keep (or `NULL` to keep all)
# KEEP_YEARS <- seq(2010, 2013) # used for CRAN
KEEP_YEARS <- NULL # used on GitHub


# read data
airports <- read.csv(FILE_CSV_AIRPORTS)
airports$Timezone <- as.numeric(airports$Timezone)
rownames(airports) <- airports$IATA

flightCounts <- readRDS(FILE_FLICHT_COUNTS)
delays <- readRDS(FILE_DELAYS)


# make sure everything is ordered by IATA
iatas1 <- rownames(airports)
airports <- airports[order(iatas1),]

iatas2a <- dimnames(flightCounts)[[1]]
iatas2b <- dimnames(flightCounts)[[2]]
flightCounts <- flightCounts[order(iatas2a), order(iatas2b),]

iatas3 <- dimnames(delays)[[2]]
delays <- delays[,order(iatas3),]


# make sure the list of IATAs is the same for all objects
iatasList <- list(
  rownames(airports),
  dimnames(flightCounts)[[1]],
  dimnames(flightCounts)[[2]],
  dimnames(delays)[[2]]
)

for(i in seq_along(iatasList)){
  for(j in seq_along(iatasList)){
    if(!identical(iatasList[[i]], iatasList[[j]])){
      stop('IATAS not identical:', i, j)
    }
  }
}


# select only some years (-> smaller dataset on CRAN)
if(!is.null(KEEP_YEARS)){
  flightCounts <- flightCounts[,,as.character(KEEP_YEARS)]
  dates <- as.Date(dimnames(delays)[[1]])
  years <- format(dates, '%Y')
  indDelays <- years %in% KEEP_YEARS
  delays <- delays[indDelays,,]
}


# create and save rda object
flights <- list(
  airports = airports,
  flightCounts = flightCounts,
  delays = delays
)

save(
  flights,
  file=OUT_FILE,
  compress = 'bzip2'
)


