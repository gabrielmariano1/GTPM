a <- c("b", "a", "b", "a", "b", "a", "b", "c", "C")
b <- c("a", "b", "j", "c", "c")

a <- c("G1 Entry", "G1 Exit",  "G1 Entry", "G1 Exit" )
b <- c("Reception Exit",  "Reception Entry")
print(LCS(a, b))

x <- lcs_(a, b, verbose )

lcs_ <- function(a, b, verbose = FALSE)
{
  LCS <- LCS(a, b)
  map_a <- list(list())
  v_a <- list(list())
  map_b <- list(list())
  v_b <- list(list())
  
  if(LCS$LLCS == 0)
  {
    returnList <- list(LCS = c(), LLCS = 0, v_a = list(), v_b = list())
    return(returnList)
  }
  
  
  #retorna todas as ocorrencias de LCS em A e todas as ocorrencias de LCS em B
  
  for(idx_lcs in 1:LCS$LLCS)
  {
    map_a[[idx_lcs]] <- numeric()
    
    for(idx_a in 1:length(a))
    {
      if(LCS$LCS[[idx_lcs]] == a[idx_a])
      {
          map_a[[idx_lcs]] <- c(map_a[[idx_lcs]], idx_a)
      }
    }
  }
  
  v <- list()
  v_a <- ps(c(), map_a, v)
  
  for(idx_lcs in 1:LCS$LLCS)
  {
    map_b[[idx_lcs]] <- numeric()
    
    for(idx_b in 1:length(b))
    {
      if(LCS$LCS[[idx_lcs]] == b[idx_b])
      {
        map_b[[idx_lcs]] <- c(map_b[[idx_lcs]], idx_b)
      }
    }
  }
  
  v <- list()
  v_b <- ps(c(), map_b, v)
  
  returnList <- list(LCS = LCS$LCS, LLCS = LCS$LLCS, v_a = v_a, v_b = v_b)
  
  
  return(returnList) 
}


ps <- function(p, map, v)
{
  if(length(map) > 0)
    for(idx_m in 1:length(map[[1]]))
    {
      if((length(p) == 0) || (map[[1]][idx_m] > max(p)))
      {
        p_ <- c(p, map[[1]][idx_m])
        map_ <- map[-1] #remove first row
        v<- ps(p_, map_, v)
      }
    }
  else
  {
    v[[length(v)+1]] <- p
  }

  return(v)
}
