#'randd()
#'
#'This function generates a random design from a bib design D(v0,b,k0,r) with 1 control added to each block
#'@param v0 b k0 r b>=v+1
#'@keywords
#'#export
#' examples
#' randd()

randd <- function(v0,b,k0,r){
  library(ibd)
  #adding 1 control to each block
  k = k0 + 1 
  v = v0 + 1
  t = r%%k
  lambda = r*(k0-1)/(v0-1)
  
  D = bibd(v0,b,r,k0,lambda,pbar=FALSE)$design
  
  #add control
  D1 = cbind(rep(0,b),D) 
  a = matrix(NA,nrow = NROW(D1),ncol = k)
  for(i in 1:NROW(D1)){
    a[i,] = sample(D1[i,])
  }
  return(a)
  
}

