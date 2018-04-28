
#### 
#### version correlated theta + multi iterations(not updating g_part)
##  ##

####code _functions same as code_ver fin3.R ####

library(truncnorm)
library(MASS)
library(prodlim)
library(profvis)
library(foreach)
library(doParallel)
library(plyr)

tmp = diag(0.1,10,10)
for (i in 1:10){
  for (j in i:10){
    if (j<10){
      tmp[i,j+1] = tmp[i,j]*0.9
    }
    tmp[j,i] = tmp[i,j]
  }
}

theta1 = mvrnorm(n = 1, mu = rep(0,10), Sigma = tmp, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
mu = rnorm(mean = qnorm(p = 0.05), sd = sqrt(0.1), n = 1)
theta= matrix(c(mu,theta1),nrow = 11)
order(pnorm(X%*%theta),decreasing = T)[1:20] 
#26  6 27 31  7 21 11  1 32 22 12  2 28  8 33 23 13  3 30 29
sort(pnorm(X%*%theta),decreasing = T)[1:20] #11 12 31 32 13  6 14  7 51  1 52 33 26 34  2 27 21  8  9 22
#  [1] 0.7667918 0.7598655 0.7409427 0.7343259 0.7336213 0.7320760 0.7269110 0.7246300 0.7067386 0.7043787 0.6989662 0.6965781 0.6904190 0.6824590 0.6534551 0.6509254 0.6451330
# [18] 0.6425816 0.6028165 0.6021526

X = matrix(c(rep(1,120),#int
             rep(1,60),rep(0,60),#factor1: 1
             rep(1,20),rep(0,40),rep(1,20),rep(0,40),#factor2: 1
             rep(0,20),rep(1,20),rep(0,20),rep(0,20),rep(1,20), rep(0,20),#factor2: 2
             rep(c(rep(1,5),rep(0,15)),6),#factor3: 1
             rep(c(rep(0,5),rep(1,5),rep(0,10)),6),#factor3: 2
             rep(c(rep(0,10),rep(1,5),rep(0,5)),6),#factor3: 3
             rep(c(1,0,0,0,0),24),#factor4: 1
             rep(c(0,1,0,0,0),24),#factor4: 2
             rep(c(0,0,1,0,0),24),#factor4: 3
             rep(c(0,0,0,1,0),24)#factor4: 
), nrow=120, ncol=11)

#### DOE function ##

doe = function(design_stg1, theta_stg1, n1 ,pp,space,ind_no){
 
  #######  generate design
  
  d1d <- opt_design(design_stg1,theta_stg1 ,n1 ,pp=0,space)
  
  
  #######  execute the experiments using the generated design 
  
  
 
  ##   data preparation
  X1 = unique(X[,ind_no])
  yt = c()
  config = c()
  nt = c()
  # mu0 = integer(NROW(d1d$design_exact))
  
  for (i in 1:NROW(d1d$design_exact)){
    xt_ = matrix(d1d$design_exact[i,],ncol = 11) ###
    #xt_ = matrix(d1d$design_exact[i,ind_no],ncol = length(ind_no)) ###
    
    config[i] = row.match(as.numeric(xt_[ind_no]),X1)
    
    #yt[i] = rbinom(n=1, size= d1d$weight_exact[i],prob = pnorm(xt_%*%theta[ind_no]))
    yt[i] = rbinom(n=1, size= d1d$weight_exact[i],prob = pnorm(xt_%*%theta))
    nt[i] = d1d$weight_exact[i]
    
  }
  
  
  exp = data.frame(config = config, y = yt, n = nt)#,mu0 = mu0)
  
  
  
  if(length(exp$config)!=length(unique(exp$config))){
    df = aggregate(x = cbind(exp$y,exp$n), by = list(exp$config),FUN = 'sum')
    names(df) = c("config", "y","n")
    exp = df
  }
  
  
  #total number of levels k0
  k0 = NROW(X1)
  
  if (NROW(exp)<k0){
    exp1 = data.frame(config = setdiff(1:k0,exp$config),y=0,n=0)
    exp = rbind(exp,exp1)
    
  }
  exp = exp[order(exp$config),]
  rownames(exp) <- NULL
  # mu0 = exp$mu0
  mu0 = rep(0,k0)
  
  design0 = c()
  for (i in 1:NROW(d1d$design_exact)){
    design0 = rbind(design0,matrix(rep(d1d$design_exact[i,],d1d$weight_exact[i]), byrow = T, ncol =  dim(theta)[1]))
  }
  
 return(list(exp = exp,design = design0))  
}
#### the iterations ####

N=30
arm_opt_mult_rpm = integer(N)
n1 = integer(N)
for(R in 1:N){
  
  theta_prior = matrix(rnorm(n = 11,mean=rep(0,11), sd = rep(1,11)),nrow = 11)
  #w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
  theta_stg1 = theta_prior
  for (j in 1:10){
    
    exp1 = doe(design_stg1, theta_stg1, n1 = 100 ,pp = 0,space = X,ind_no = c(1:11))  #design_stg1 update only for previous step?
    #exp1 = doe(design_stg1 =exp1$design, theta_stg1, n1 = 100 ,pp = 0,space = X,ind_no = c(1:11))
    exp = exp1$exp
    n_acc = exp$n
    y_acc = exp$y
    
    
    k0 = 120
    mu0 = rep(0,k0)
    
    
    #calculate theta: (similar probability matching method, to make it comparable to TS)
    theta_hat = apply(sim.post(nruns = 2000,ndraws = 500,n_acc,y_acc,X,theta_stg1, last_N = 500,k0),MARGIN = 1,FUN = mean) ###
    #data.frame(theta_hat,theta,theta_prior)
    #sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]
    #order(pnorm(X%*%theta_hat),decreasing = T)[1:20]
    theta_stg1 = theta_hat
    design_stg1 = exp1$design
  }
  
  
  
  #######phase 3 end ###
  
  
  arm_opt_mult_rpm[R] = order(pnorm(X%*%theta_hat),decreasing = T)[1]
  #if( arm_opt_mult_rpm[R]!=13) browser()
  
} 

table(arm_opt_mult_rpm)
#arm_opt_mult_rpm  8+300
#1  6  7 11 12 13 21 26 27 28 31 32 
#1  3  1  1  5  1  2  8  2  3  1  2 
