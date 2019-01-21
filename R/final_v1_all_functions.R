# 01/21/2019

#Create file as final version of the SDR generation method combining stand_new and code_bioassay


########################  prev versions ########################################


# SDR_ext() : in code_bioassay.R
# SDR_find(): in SDR_find_new.R   alias SDR_find1()
# ext_yd() : in stand_new.R no control ext_yd2() ; in correction.R with 1 control ext_yd1()
# informat(): in code_bioassay.R 1 control ; in stan_new.R no control informat2()
# eff() : in stand_new.R no control eff2() ; in code_bioassay.R 1 control eff()
#         (actually same, one with informat2 and one with informat, MERGE )

########################   required packages ########################################

library(ibd)
library(data.table)
library(reshape2)
library(MASS)

##############  gen_eyd(v,b,k, ctrl=1)   ###############


# Generate an EYD design with ctrl controls required in each block, default 1.
# v: number of treatments and controls
# b: number of blocks
# k: block size
# ctrl: number of controls required in each block, 1 for 1 control, 0 for no control

gen_eyd <- function(v,b,k, ctrl){

  #adding 1 control to each block

  if(ctrl == 0){
    r = b*k/v
    try(if(((b*k) %%v>0)|((r*(k-1))%%(v-1)>0)|(b<v)) stop("No valid design available"))
    k0 = k
    v0 = v
    t = r%%k
    lambda = r*(k-1)/(v-1)

    D = bibd(v,b,r,k,lambda,pbar=FALSE)$design
    D1 = D
    s = (k*v-t*v)/k # assume starting from a bib design


    ### add dummy columns
    dummy = NULL
    for (i in 1:v0){
      dummy = c(dummy,rep(i,k-t))
    }
    D2 = rbind(D1,matrix(dummy,nrow = s))


    ### re-label trtments/control
    D3 = D2

    for (j in 1:v0){ #0 no control
      # shuffle the indices of original blocks in which treatment j exist
      ind  = sample(which((apply(D3==j,FUN = sum, MARGIN = 1) * c(1:(b+s))>0 )* (apply(D3==j,FUN = sum,MARGIN = 1) * c(1:(b+s))<b+1)==1))

      if (j>0) {
        m = floor(r/k)
        r1 = r+k-t
        tt=t
      }else {
        m = floor(b/k)
        r1 = b+k-t0
        tt=t0
      }

      if(m==1 ){
        if (tt>0){
          D3[ind[(m*k+1):min((m+1)*k,length(ind))],] = ifelse(D3[ind[(m*k+1):min((m+1)*k,length(ind))],] ==j, j + m*v ,D3[ind[(m*k+1):min((m+1)*k,length(ind))],])

        }
        D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==j, j + v ,D3[(b+1):(b+s),])
      }else if (m>1) {
        for(i in 1:(m-1)){
          D3[ind[(i*k+1):min((i+1)*k,length(ind))],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,length(ind))],] ==j, j + i*v ,D3[ind[(i*k+1):min((i+1)*k,length(ind))],])
        }
        if (tt>0){
          D3[ind[(m*k+1):min((m+1)*k,length(ind))],] = ifelse(D3[ind[(m*k+1):min((m+1)*k,length(ind))],] ==j, j + m*v ,D3[ind[(m*k+1):min((m+1)*k,length(ind))],])

        }
        D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==j,j + m*v ,D3[(b+1):(b+s),])
      }

    }
    #fill row
    new_fill = NULL
    j=1
    remain = D3

    while (j <=k){
      new_fill = cbind(new_fill,SDR_find(remain)) # in SDR_find_new
      remain = matrix(NA,nrow = b+s, ncol = k-dim(new_fill)[2] )
      for ( i in 1:(b+s)){
        remain[i,] = c(fsetdiff(data.table(D3[i,]),data.table(new_fill[i,]),all = TRUE),recursive = T,use.names = F)
      }

      if(j<k && SDR_ext(remain)==T) {  #SDR_ext in code_bioassay
        j = j+1
      }
      else if(j==k) break
      else  if(j<k &&  SDR_ext(remain)==F){
        j=1
        new_fill = NULL}
    }


    #transform back to original treatment label
    a=new_fill
    if(r/k>0){#|| b/k>0){
      a=(a%%v)[1:b,]
    }

  }else {
    # 1 control
    k0 = k - 1
    v0 = v - 1
    r = b*k0/v0
    try(if(((b*k0) %%v0>0)|((r*(k0-1))%%(v0-1)>0)|(b<v0)) stop("No valid design available"))

    t = r%%k
    t0 = b%%k

    lambda = r*(k0-1)/(v0-1)
    D = bibd(v0,b,r,k0,lambda,pbar=FALSE)$design

    #add control
    D1 = cbind(rep(0,b),D)
    s = (k*v-t*v0-b%%k)/k # assume starting from a bib design

    #add dummy columns
    dummy = rep(0,k-b%%k)
    for (i in 1:v0){
      dummy = c(dummy,rep(i,k-t))
    }
    D2 = rbind(D1,matrix(dummy,nrow = s))


    #re-label trtments/control
    D3 = D2
    for (j in 0:v0){
      # shuffle the indices of original blocks in which treatment j exist
      ind  = sample(which((apply(D3==j,FUN = sum, MARGIN = 1) * c(1:(b+s))>0 )* (apply(D3==j,FUN = sum,MARGIN = 1) * c(1:(b+s))<b+1)==1))
      if (j>0) {
        m = floor(r/k)
        r1 = r+k-t
        tt=t
      }else {
        m = floor(b/k)
        r1 = b+k-t0
        tt=t0
      }

      if(m==1 ){
        if (tt>0){
          D3[ind[(m*k+1):min((m+1)*k,length(ind))],] = ifelse(D3[ind[(m*k+1):min((m+1)*k,length(ind))],] ==j, j + m*v ,D3[ind[(m*k+1):min((m+1)*k,length(ind))],])
        }
        D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==j, j + v ,D3[(b+1):(b+s),])
      }else if (m>1) {
        for(i in 1:(m-1)){
          D3[ind[(i*k+1):min((i+1)*k,length(ind))],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,length(ind))],] ==j, j + i*v ,D3[ind[(i*k+1):min((i+1)*k,length(ind))],])
        }
        if (tt>0){
          D3[ind[(m*k+1):min((m+1)*k,length(ind))],] = ifelse(D3[ind[(m*k+1):min((m+1)*k,length(ind))],] ==j, j + m*v ,D3[ind[(m*k+1):min((m+1)*k,length(ind))],])
        }
        D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==j,j + m*v ,D3[(b+1):(b+s),])
      }

    }

    #fill row
    new_fill = NULL
    j=1
    remain = D3

    while (j <=k){
      new_fill = cbind(new_fill,SDR_find(remain)) # using new SDR_find function
      remain = matrix(NA,nrow = b+s, ncol = k-dim(new_fill)[2] )
      for ( i in 1:(b+s)){
        remain[i,] = c(fsetdiff(data.table(D3[i,]),data.table(new_fill[i,]),all = TRUE),recursive = T,use.names = F)
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
  }


  return(a)
}


# gen_eyd(6,20,3,0)
# gen_eyd(5,6,3,1)
# gen_eyd(7,10,4,1)



##############   SDR_find(x, block_by_row = T) ###############

# Find a SDR of a matrix
# x: matrix
# block_by_row: indication of orientation of the finding

SDR_find <- function(x, block_by_row = T){
  if(!block_by_row){
    x = t(x)
  }
  nr = dim(x)[1]
  nc = dim(x)[2]

  # Build a record for number of remaining available counts in the play
  sym = sort(unique(c(x)))
  counts = list()
  for (i in sym){
    counts[[as.character(i)]] = nc
  }


  fill = NULL
  ct_update = counts
  cand = x[1,]
  i=1
  while (i <= (nr)){
    prob = NULL
    for (j in 1:length(cand)){prob[j] = ct_update[[as.character(cand[j])]]}
    tmp = cand[sample(prob = prob/(sum(prob)), x = length(cand), size = 1)]
    if (i<nr){
      sdr_remain = length(setdiff(c(x[(i+1):nr,]), c(fill,tmp)))>=nr-i
    }else{
      sdr_remain = TRUE
    }


    if (sdr_remain && (! (tmp%in% fill))){
      fill[i] = tmp

      # update counts
      ct_update[[as.character(tmp)]] = ct_update[[as.character(tmp)]] - 1
      i = i+1
      if (i<=nr){cand = x[i,]} else cand = NULL

      next
    }else if (!sdr_remain){
      fill = NULL
      i = 1
      ct_update = counts
      cand = x[1,]
    }else{
      cand = setdiff(cand,tmp)

    }
    if (length(cand)==0){
      cand = x[1,]
      fill = NULL

      ct_update = counts
      #print(paste("start over at i=", i))
      i = 1
    }

  }
  return(fill)

}


##############   SDR_ext(x, block_by_row = T) ###############

# Check if the remaining design has SDR for each block
# x: matrix
# block_by_row: indication of orientation of the finding

###
SDR_ext <- function(x, block_by_row = T){
  if (block_by_row!=T){
    x = t(x)
  }
  if(length(unique(as.vector(x)))<dim(x)[1]){
    return(F)
  }
  else return(T)
}

##############   eff(design,v,b,k,ctrl, type = "D") ###############

# Calculate the relative efficiency of design with given optimality criterion.
# Definition of RE see paper.

# v: number of treatments and controls
# b: number of blocks
# k: block size
# ctrl: number of controls required in each block, 1 for 1 control, 0 for no control


eff<-function(design,v,b,k,ctrl, type = "D"){

  if(type=="D"){
    info = informat(design,v,b,k,ctrl)
    return(det(info$Md)/det(info$upper_M))
  }else if(type=="A"){
    info =informat(design,v,b,k,ctrl)
    return(sum(diag(info$Md))/sum(diag(info$upper_M)))
   }else if(type=="E"){
    info = informat(design,v,b,k,ctrl)
    return(max(diag(info$Md))/max(diag(info$upper_M)))
   }


}


##############   informat(design,v,b,k,ctrl)  ###############

# Calculate the information matrix of the design

# v: number of treatments and controls
# b: number of blocks
# k: block size
# ctrl: number of controls required in each block, 1 for 1 control, 0 for no control




informat <- function(design,v,b,k,ctrl){
  if (ctrl ==0){
    r = b*k/v
    k0 = k
    v0 = v
  }else{
    k0 = k - 1
    v0 = v - 1
    r = b*k0/v0
  }
  long = melt(t(design))

  long_b = long[,2:3]
  long_r = long[,c(1,3)]

  N_b = t(table(long_b)) # incidence matrix for block(plate)
  N_r = t(table(long_r)) # incidence matrix for row

  # Calculate  information matrix for beta
  if (ctrl ==0){
    Cd = diag(c(rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(rep(r,v0)),ncol = 1)%*%matrix(c(rep(r,v0)),ncol = v)/b/k
  }else{
    Cd = diag(c(b,rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(b,rep(r,v0)),ncol = 1)%*%matrix(c(b,rep(r,v0)),ncol = v)/b/k

  }
  # Calculate upper bound
  Td = matrix(0,nrow = b*k, ncol = v)
  U = matrix(0,nrow = b*k, ncol = b)
  P = matrix(0,nrow = b*k, ncol = k)
  for(i in 1:NROW(long_b)){
    Td[i,long_b[i,]$value+1] = 1
    U[i,ceiling(i/k)]=1
    P[i,ifelse(i%%k==0,k,i%%k)] =1
  }

  upper = t(Td)%*%(diag(rep(1,b*k))-U%*%ginv(t(U)%*%U)%*%t(U))%*%Td
  Md = Cd[2:v,2:v]
  upper_M = upper[2:v,2:v]

  return(list(Cd = Cd,upper = upper,Md = Md,upper_M = upper_M))
}

##############   shuffle(x)  ###############

# Generate a randome shuffle within each row of a matrix

# x: matrix

shuffle <- function(x){

  #add control
  #D1 = cbind(rep(0,b),D)
  a = matrix(NA,nrow = NROW(x),ncol = dim(x)[2])
  for(i in 1:NROW(x)){
    a[i,] = sample(x[i,])
  }
  return(a)

}


## test

# design = gen_eyd(7,10,4,1); N = 1000
# randeff = c()
# for(i in 1:N){
#   randeff[i] = eff(shuffle(design),7,10,4,1,type = "D")
# }
# summary(randeff)
# eff(design,7,10,4,1,type = "D")


# design = gen_eyd(6,20,3,0); N = 1000
# randeff = c()
# for(i in 1:N){
#   randeff[i] = eff(shuffle(design),6,20,3,0,type = "D")
# }
# summary(randeff)
# eff(design,6,20,3,0,type = "D")

