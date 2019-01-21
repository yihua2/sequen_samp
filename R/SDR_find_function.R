#'SDR_find()
#'
#'This function is a auxiliary function for ext_yd. It finds the system of distinct representatives for each block of design x
#'@param x block_by_row(Default to True)
#'@keywords
#'#export
#' examples
#' SDR_find()

SDR_find <- function(x, block_by_row = T){
  nr = dim(x)[1]
  nc = dim(x)[2]
  if (nc>1){
    fill = NULL
   # i=1
    remain = matrix(0,nrow = nr,ncol = nc-1)
    ind = 0
    while (ind==0){
      for (i in 1:nr){
        #fill[i] = sample(setdiff(x[i,],fill))[1]
        fill[i] = sample(x[i,])[1]
      }
      if (length(unique(fill))==length(fill)) {
        ind=1
      }
    }
    return(fill)
  }
  else return(x)

}



