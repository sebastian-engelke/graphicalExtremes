
# Memoise functions
myHash <- function(x){
  if(igraph::is_igraph(x)){
    x <- igraph::as_adjacency_matrix(x)
  } else if(is.list(x)){
    x <- lapply(x, myHash)
  }
  return(rlang::hash(x))
}
myCache <- cachem::cache_disk(
  'memoise_cache',
  # logfile = stdout(),
  evict = 'fifo'
)
myMemoise <- function(f){
  if('memoised' %in% class(f)){
    return(f)
  }
  return(memoise::memoise(
    f,
    hash = myHash,
    cache = myCache
  ))
}

# Memoise expensive functions:
complete_Gamma <- myMemoise(complete_Gamma)
loglik_HR <- myMemoise(loglik_HR)
emp_chi <- myMemoise(emp_chi)
emp_vario <- myMemoise(emp_vario)
eglearn <- myMemoise(eglearn)
emp_vario_pairwise <- myMemoise(emp_vario_pairwise)
emp_chi_pairwise <- myMemoise(emp_chi_pairwise)

