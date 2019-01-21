### Generate system of distinct representatives


### Start with bib design (v,b,k,r) b>=v+1, adding 1 control to each block
## Future: add b>=v+1 condition
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


  if(t!=0 && b%%k!=0){

    s = (k*v-t*v0-b%%k)/k # assume starting from a bib design

    #add dummy columns
    dummy = rep(0,k-b%%k)
    for (i in 1:v0){
      dummy = c(dummy,rep(i,k-t))
    }
    D2 = rbind(D1,matrix(dummy,nrow = s))

  } else if (b%%k==0 ){
  	s = (k*v-t*v0-b%%k-k)/k # assume starting from a bib design

    #add dummy columns
    dummy = c()
    for (i in 1:v0){
      dummy = c(dummy,rep(i,k-t))
    }
    D2 = rbind(D1,matrix(dummy,nrow = s))
  }else {
    D2=D1
    s=0}


    #re-label trtments
    #trt 0
    r0 = b+k-b%%k
    #m = floor(r0/k)
    m = floor(b/k)
    #t0 = (r0)%%k
    t0 = b%%k
    D3 = D2

    #apply(D3==0,FUN = sum,MARGIN = 1) * c(1:15) # blocks that this trtment is in

    ind  = sample(which((apply(D3==0,FUN = sum, MARGIN = 1) * c(1:(b+s))>0 )* (apply(D3==0,FUN = sum,MARGIN = 1) * c(1:(b+s))<b+1)==1))



    if(m>0 && t0>0){
      for(i in 1:m){
        D3[ind[(i*k+1):min((i+1)*k,b)],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,b)],] ==0, 0 + i*v ,D3[ind[(i*k+1):min((i+1)*k,b)],])
      }
      D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==0, 0 + i*v ,D3[(b+1):(b+s),])

    }else if (t0==0 && m>1){
      for(i in 1:(m-1)){
        D3[ind[(i*k+1):min((i+1)*k,b)],] = ifelse(D3[ind[(i*k+1):min((i+1)*k,b)],] ==0, 0 + i*v ,D3[ind[(i*k+1):min((i+1)*k,b)],])
      }
    }else if (m==1 && t0 == 0){
    	D3[(b+1):(b+s),] = ifelse(D3[(b+1):(b+s),] ==0, 0 + i*v ,D3[(b+1):(b+s),])
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


    #################  Preparation Completed  ###############


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

### to test if the remaining design has SDR for each block
SDR_ext <- function(x, block_by_row = T){
  if (block_by_row!=T){
    x = t(x)
  }
  if(length(unique(as.vector(x)))<dim(x)[1]){
    return(F)
  }
  else return(T)
}


### generate SDR

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


#generate a random design with control
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

### Transform to incidence matrix
informat <- function(design,v0,b,k0,r){
  #design1 = as.data.frame(design,row.names = c('b1','b2','b3','b4','b5','b6','b7','b8','b9','b10'))
  #design1 = as.data.frame(design,row.names = c('b1','b2','b3','b4','b5','b6'))
  k = k0+1
  v = v0+1

  library(reshape2)
  #table(melt(design1)[-1])
  #table(design)

  long = melt(t(design))

  long_b = long[,2:3]
  long_r = long[,c(1,3)]

  N_b = t(table(long_b)) # incidence matrix for block(plate)
  N_r = t(table(long_r)) # incidence matrix for row

  # Calculate  information matrix for beta
  Cd = diag(c(b,rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(b,rep(r,v0)),ncol = 1)%*%matrix(c(b,rep(r,v0)),ncol = v)/b/k
  # Calculate upper bound
  Td = matrix(0,nrow = b*k, ncol = v)
  U = matrix(0,nrow = b*k, ncol = b)
  P = matrix(0,nrow = b*k, ncol = k)
  for(i in 1:NROW(long_b)){
    Td[i,long_b[i,]$value+1] = 1
    U[i,ceiling(i/k)]=1
    P[i,ifelse(i%%k==0,k,i%%k)] =1
  }


  upper = t(Td)%*%(diag(rep(1,b*k))-U%*%MASS::ginv(t(U)%*%U)%*%t(U))%*%Td
  Md = Cd[2:v,2:v]
  upper_M = upper[2:v,2:v]

  return(list(Cd = Cd,upper = upper,Md = Md,upper_M = upper_M))
}

eff<-function(design,v0,b,k0,r,type = "D"){
  if(type=="D"){
    info = informat(design,v0,b,k0,r)
    #print("D efficiency")
    #return(det(solve(info$upper_M))/det(solve(info$Md)))
    return(det(info$Md)/det(info$upper_M))
  }else if(type=="A"){

      info = informat(design,v0,b,k0,r)
      #print("A efficiency")
      #return(sum(diag(solve(info$upper_M)))/sum(diag(solve(info$Md))))
      return(sum(diag(info$Md))/sum(diag(info$upper_M)))
  }else if(type=="E"){

    info = informat(design,v0,b,k0,r)
    #print("E efficiency")
    #return(max(diag(solve(info$upper_M)))/max(diag(solve(info$Md))))
    return(max(diag(info$Md))/max(diag(info$upper_M)))
  }


}

#Test

#k=4;v=6+1;b=10;s=5
time = system.time({design = ext_yd(v0=6,b=10,k0=3,r=5)})

N = 1000
randeff = c()
for(i in 1:N){
  randeff[i] = eff(randd(6,10,3,5),6,10,3,5,type = "D")
}
summary(randeff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.09522 0.42600 0.51200 0.50310 0.58260 0.81040

eff(design,6,10,3,5,type = "D")
#[1] "D efficiency"
#[1] 0.8889194


design1 = ext_yd(v0=4,b=6,k0=2,r=3)
N = 1000
randeff1 = c()
for(i in 1:N){
  randeff1[i] = eff(randd(4,6,2,3),4,6,2,3,type = "D")
}
summary(randeff1)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.05952 0.31888 0.46769 0.43795 0.56463 1.00000

eff(design1,4,6,2,3,type = "D")
#[1] "D efficiency"
#[1] 1

design2 = ext_yd(v0=5,b=5,k0=4,r=4)
N = 1000
randeff2 = c()
for(i in 1:N){
  randeff2[i] = eff(randd(5,5,4,4),5,5,4,4,type = "D")
}
summary(randeff2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.01449 0.19595 0.28084 0.28029 0.36762 0.64770

eff(design2,5,5,4,4,type = "D")
#[1] 0.8055187

#A

N = 1000
randeff = c()
for(i in 1:N){
  randeff[i] = eff(randd(6,10,3,5),6,10,3,5,type = "A")
}
summary(randeff)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.5016  0.8394  0.8807  0.8697  0.9108  0.9676

eff(design,6,10,3,5,type = "A")
#[1] "A efficiency"
#[1] 0.9820226


design1 = ext_yd(v0=4,b=6,k0=2,r=3)
N = 1000
randeff1 = c()
for(i in 1:N){
  randeff1[i] = eff(randd(4,6,2,3),4,6,2,3,type = "A")
}
summary(randeff1)

#      Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.4732  0.5426  0.6840  0.6930  0.8621  0.9347
eff(design1,4,6,2,3,type = "A")
#[1] "A efficiency"
#[1] 1

design2 = ext_yd(v0=5,b=5,k0=4,r=4)
N = 1000
randeff2 = c()
for(i in 1:N){
  randeff2[i] = eff(randd(5,5,4,4),5,5,4,4,type = "A")
}
summary(randeff2)



eff(design2,5,5,4,4,type = "A")

#E

N = 1000
randeff = c()
for(i in 1:N){
  randeff[i] = eff(randd(6,10,3,5),6,10,3,5,type = "E")
}
summary(randeff)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.3390  0.7200  0.7847  0.7706  0.8390  0.9406


eff(design,6,10,3,5,type = "E")
#[1] "D efficiency"
#[1] 0.9774513


design1 = ext_yd(v0=4,b=6,k0=2,r=3)
N = 1000
randeff1 = c()
for(i in 1:N){
  randeff1[i] = eff(randd(4,6,2,3),4,6,2,3,type = "E") #singularity?
}
summary(randeff1)

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.1385  0.4858  0.6982  0.6268  0.7905  1.0000

eff(design1,4,6,2,3,type = "E")
#[1] "D efficiency"
#[1] 1

design2 = ext_yd(v0=5,b=5,k0=4,r=4)
N = 1000
randeff2 = c()
for(i in 1:N){
  randeff2[i] = eff(randd(5,5,4,4),5,5,4,4,type = "E")
}
summary(randeff2)



eff(design2,5,5,4,4,type = "E")

