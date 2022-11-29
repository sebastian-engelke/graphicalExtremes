
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}

library(ggplot2)


plotFlights(1:5, clipMap = 1)
