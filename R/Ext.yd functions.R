#'ext_yd()
#'
#'This function generates the extended Youden Square from a bib design D(v0,b,k0,r) with 1 control added to each block
#'@param v0 b k0 r b>=v+1
#'@keywords
#'@export
#' examples
#' ext_yd()

ext_yd <- function(v0,b,k0,r){
  library(ibd)
  #adding 1 control to each block
  k = k0 + 1 
  v = v0 + 1
  t = r%%k
  lambda = r*(k0-1)/(v0-1)
  
  D = bibd(v0,b,r,k0,lambda,pbar=FALSE)$design
  
  #add control
  D1 = cbind(rep(0,b),D) 
  
  if(t!=0 || b%%k!=0){
    s = (k*v-t*v0-b%%k)/k # assume starting from a bib design
    
    #add dummy columns
    dummy = rep(0,k-b%%k)
    for (i in 1:v0){
      dummy = c(dummy,rep(i,k-t))
    }
    D2 = rbind(D1,matrix(dummy,nrow = s))
    
  } else{
    D2=D1
    s=0}
  
  
  #re-label trtments
  #trt 0 
  m = floor(b/k)
  t0 = b%%k
  D3 = D2
  
  #apply(D3==0,FUN = sum,MARGIN = 1) * c(1:15) # blocks that this trtment is in
  
  ind  = sample(which((apply(D3==0,FUN = sum, MARGIN = 1) * c(1:(b+s))>0 )* (apply(D3==0,FUN = sum,MARGIN = 1) * c(1:(b+s))<b+1)==1))
  
  if(t0>0){
    for(i in 1:m){
      D3[ind[(i*k+1):min((i+1)*k,b)],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,b)],] ==0, 0 + i*v ,D3[ind[(i*k+1):min((i+1)*k,b)],])
    }
    D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==0, 0 + i*v ,D3[(b+1):(b+s),])
    
  }else{
    for(i in 1:(m-1)){
      D3[ind[(i*k+1):min((i+1)*k,b)],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,b)],] ==0, 0 + i*v ,D3[ind[(i*k+1):min((i+1)*k,b)],])
    }
  }
  
  #test treatments
  for (j in 1:v0){
    ind  = sample(which((apply(D3==j,FUN = sum, MARGIN = 1) * c(1:(b+s))>0 )* (apply(D3==j,FUN = sum,MARGIN = 1) * c(1:(b+s))<b+1)==1))
    #m = floor(r/k) 
    if (s>0){
      r1 = r+k-t
    }else{
      r1 = r
    }
    m = floor(r1/k) 
    if (m>1 && r1%%k>0){
      for(i in 1:m){
        D3[ind[(i*k+1):min((i+1)*k,length(ind))],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,length(ind))],] ==j, j + i*v ,D3[ind[(i*k+1):min((i+1)*k,length(ind))],])
      }
      # D3[ind[(k+1):length(ind)],] = ifelse(D3[ind[(k+1):length(ind)],] ==i, i + v ,D3[ind[(k+1):length(ind)],])
      D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==j,j + i*v ,D3[(b+1):(b+s),])
    }else if (m>1 && r1%%k==0){
      for(i in 1:(m-1)){
        D3[ind[(i*k+1):min((i+1)*k,length(ind))],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,length(ind))],] ==j, j + i*v ,D3[ind[(i*k+1):min((i+1)*k,length(ind))],])
      }
      D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==j, j + i*v ,D3[(b+1):(b+s),])
      
    }
    
  }
  
  
  #################  Preparation Completed 
  
  
  new_fill = NULL 
  
  #fill row   
  j=1
  remain = D3
  while (j <=k){
    new_fill = cbind(new_fill,SDR_find(remain)) 
    remain = matrix(NA,nrow = b+s, ncol = k-dim(new_fill)[2] )
    for ( i in 1:(b+s)){
      remain[i,] = setdiff(D3[i,],new_fill[i,])
    }
    
    if(j<k && SDR_ext(remain)==T) {
      j = j+1
    }
    else if(j==k) break
    else  if(j<k &&  SDR_ext(remain)==F){
      j=1
      new_fill = NULL}
  }
  
  
  #transform back to original treatment label
  a=new_fill
  if(r/k>0|| b/k>0){
    a=(a%%v)[1:b,]
  }
  return(a)
  
}



#'SDR_ext()
#'
#'This function is a auxiliary function for ext_yd. It tests if a system of distinct representative exists for each block of design x
#'@param x block_by_row(Default to True)
#'@keywords
#'@export
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




#'SDR_find()
#'
#'This function is a auxiliary function for ext_yd. It finds the system of distinct representatives for each block of design x
#'@param x block_by_row(Default to True)
#'@keywords
#'@export
#' examples
#' SDR_find()

SDR_find <- function(x, block_by_row = T){
  nr = dim(x)[1]
  nc = dim(x)[2]
  if (nc>1){
    fill = NULL
    i=1
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



#'randd()
#'
#'This function generates a random design from a bib design D(v0,b,k0,r) with 1 control added to each block
#'@param v0 b k0 r b>=v+1
#'@keywords
#'@export
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

