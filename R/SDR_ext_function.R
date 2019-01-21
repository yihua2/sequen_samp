#'SDR_ext()
#'
#'This function is a auxiliary function for ext_yd. It tests if a system of distinct representative exists for each block of design x
#'@param x block_by_row(Default to True)
#'@keywords
#'#export
#' examples
#' SDR_ext()

SDR_ext <- function(x, block_by_row = T){
  if (block_by_row!=T){
    x = t(x)
  }
  if(length(unique(as.vector(x)))<dim(x)[1]){
    return(F)
  }
  else return(T)
}




