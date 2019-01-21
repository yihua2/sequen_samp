shuffle <- function(D1){
  # library(ibd)
  # #adding 1 control to each block
  # k = k0 + 1
  # v = v0 + 1
  # t = r%%k
  # lambda = r*(k0-1)/(v0-1)
  #
  #D = bibd(v0,b,r,k0,lambda,pbar=FALSE)$design

  #add control
  #D1 = cbind(rep(0,b),D)
  a = matrix(NA,nrow = NROW(D1),ncol = dim(D1)[2])
  for(i in 1:NROW(D1)){
    a[i,] = sample(D1[i,])
  }
  return(a)

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

# EX1 :k=4;v=6+1;b=10;s=5

design = ext_yd1(v0=6,b=10,k0=3,r=5)

N = 10000

randeff = c()
for(i in 1:N){
  randeff[i] = eff(shuffle(design),6,10,3,5,type = "D")
}
summary(randeff)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.03815 0.41904 0.50285 0.49924 0.58265 0.88892

eff(design,6,10,3,5,type = "D")
#[1] "D efficiency"
#[1] 0.8889194

randeff = c()
for(i in 1:N){
  randeff[i] = eff(shuffle(design),6,10,3,5,type = "A")
}
summary(randeff)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.6956  0.8822  0.9000  0.8998  0.9178  0.9800

eff(design,6,10,3,5,type = "A")
#[1] "A efficiency"
#[1] 0.98


randeff = c()
for(i in 1:N){
  randeff[i] = eff(shuffle(design),6,10,3,5,type = "E")
}
summary(randeff)

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.8200  0.9800  0.9800  0.9689  0.9800  0.9800

eff(design,6,10,3,5,type = "E")
#[1] 0.98


# EX3
design1 = ext_yd1(v0=4,b=6,k0=2,r=3)
randeff1 = c()
for(i in 1:N){
  randeff1[i] = eff(shuffle(design1),4,6,2,3,type = "D")
}
summary(randeff1)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.05952 0.31888 0.46769 0.43795 0.56463 1.00000

eff(design1,4,6,2,3,type = "D")
#[1] "D efficiency"
#[1] 1

randeff1 = c()
for(i in 1:N){
  randeff1[i] = eff(shuffle(design1),4,6,2,3,type = "A")
}
summary(randeff1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.5833  0.7917  0.8333  0.8332  0.8750  1.0000
eff(design1,4,6,2,3,type = "A")
#1

randeff1 = c()
for(i in 1:N){
  randeff1[i] = eff(shuffle(design1),4,6,2,3,type = "E")
}
summary(randeff1)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.5833  0.7917  0.8333  0.8332  0.8750  1.0000
eff(design1,4,6,2,3,type = "E")
#1


#EX2
design2 = ext_yd1(v0=5,b=5,k0=4,r=4)

randeff2 = c()
for(i in 1:N){
  randeff2[i] = eff(shuffle(design2),5,5,4,4,type = "D")
}
summary(randeff2)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.01449 0.19595 0.28084 0.28029 0.36762 0.64770

eff(design2,5,5,4,4,type = "D")
#[1] 0.8055187


randeff2 = c()
for(i in 1:N){
  randeff2[i] = eff(shuffle(design2),5,5,4,4,type = "A")
}
summary(randeff2)

eff(design2,5,5,4,4,type = "A")

randeff2 = c()
for(i in 1:N){
  randeff2[i] = eff(shuffle(design2),5,5,4,4,type = "E")
}
summary(randeff2)

eff(design2,5,5,4,4,type = "E")

