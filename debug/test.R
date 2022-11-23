
if(!nchar(Sys.getenv('VSCODE_DEBUG_SESSION'))){
    devtools::load_all('.')
}



plotFlights(plotConnections = FALSE, map='world')


