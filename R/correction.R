v0=4;b=6;k0=2;r=3
v0=6;b=10;k0=3;r=5
v0=6;b=20;k0=3;r=10
v0=5;b=5;k0=4;r=4

ext_yd1 <- function(v0,b,k0,r){
library(ibd)
library(data.table)
#adding 1 control to each block
k = k0 + 1
v = v0 + 1
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


#################  Preparation Completed

#fill row
new_fill = NULL
j=1
remain = D3

while (j <=k){
  new_fill = cbind(new_fill,SDR_find1(remain)) # using new SDR_find function
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
return(a)
}

ext_yd1(4,6,2,3)

ext_yd1(6,10,3,5)
ext_yd1(6,20,3,10)
ext_yd1(5,5,4,4)
