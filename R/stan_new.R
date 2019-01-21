v0 = 6
b = 20
k0 = 3
k = 3
r = 10

# no control
ext_yd2 <- function(v0,b,k0,r){
library(ibd)
library(data.table)
# no control
k = k0
v = v0
t = r%%k
#t0 = b%%k
lambda = r*(k0-1)/(v0-1)

D = bibd(v0,b,r,k0,lambda,pbar=FALSE)$design

# no add control
D1 = D
s = (k*v-t*v0)/k # assume starting from a bib design

#add dummy columns
dummy = NULL
for (i in 1:v0){
  dummy = c(dummy,rep(i,k-t))
}
D2 = rbind(D1,matrix(dummy,nrow = s))


#re-label trtments/control
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


#################  Preparation Completed

#fill row
new_fill = NULL
j=1
remain = D3

while (j <=k){
  new_fill = cbind(new_fill,SDR_find1(remain)) # in SDR_find_new
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

return(a)
}

### Transform to incidence matrix
informat2 <- function(design,v0,b,k0,r){

  k = k0#+1
  v = v0#+1
  library(reshape2)

  long = melt(t(design))

  long_b = long[,2:3]
  long_r = long[,c(1,3)]

  N_b = t(table(long_b)) # incidence matrix for block(plate)
  N_r = t(table(long_r)) # incidence matrix for row

  # Calculate  information matrix for beta
 # Cd = diag(c(b,rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(b,rep(r,v0)),ncol = 1)%*%matrix(c(b,rep(r,v0)),ncol = v)/b/k
  Cd = diag(c(rep(r,v0)))-N_b%*%t(N_b)/k - N_r%*%t(N_r)/b + matrix(c(rep(r,v0)),ncol = 1)%*%matrix(c(rep(r,v0)),ncol = v)/b/k

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

eff2<-function(design,v0,b,k0,r,type = "D"){
  if(type=="D"){
    info = informat2(design,v0,b,k0,r)
    #print("D efficiency")
    #return(det(solve(info$upper_M))/det(solve(info$Md)))
    return(det(info$Md)/det(info$upper_M))

  }else if(type=="A"){

    info = informat2(design,v0,b,k0,r)
    #print("A efficiency")
    #return(sum(diag(solve(info$upper_M)))/sum(diag(solve(info$Md))))
    return(sum(diag(info$Md))/sum(diag(info$upper_M)))

  }else if(type=="E"){

    info = informat2(design,v0,b,k0,r)
    #print("E efficiency")
    #return(max(diag(solve(info$upper_M)))/max(diag(solve(info$Md))))
    return(max(diag(info$Md))/max(diag(info$upper_M)))
  }

}



design = ext_yd2(v0=6,b=20,k0=3,r=10)
eff2(design,6,20,3,10,type = "D")

design_stan = matrix(c(1,2,3,
                       4,1,2,
                       3,5,1,
                       1,4,0,
                       0,1,5,
                       2,3,0,
                       5,2,4,
                       2,0,5,
                       0,4,3,
                       3,5,4,
                       4,0,5,
                       3,0,5,
                       0,2,4,
                       5,3,2,
                       4,2,3,
                       5,4,1,
                       3,1,0,
                       1,3,4,
                       2,5,1,
                       0,1,2)
                       , nrow = 20,byrow = T)
eff2(design_stan,6,20,3,10,type = "E")
