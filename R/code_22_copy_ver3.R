
#######################
# version fin: multistage -contrast matrix identity, change theta to length 11
########################
library(truncnorm)
library(MASS)
library(prodlim)
library(profvis)
library(foreach)
library(doParallel)
library(plyr)

N= matrix(c(rep(0,60),rep(1,60), 
            rep(0,20),rep(1,20),rep(2,20),rep(0,20),rep(1,20),rep(2,20),
            rep(c(rep(0,5),rep(1,5),rep(2,5),rep(3,5)),6),
            rep(c(0,1,2,3,4), 24)), nrow= 120, ncol = 4)

########## Bernoulli Bandit For one single simulation  #########
#contrast matrix
#g_part = diag(rep(1,11))
#Contrast matrix
g_part = cbind(rep(1,11),
               rbind(rep(0,14),
                     c(-1,1,rep(0,12)),
                     cbind(matrix(0,nrow = 2, ncol = 2),rep(-1,2),diag(1,2),matrix(0,nrow = 2, ncol = 9)),
                     cbind(matrix(0,nrow = 3, ncol = 5),rep(-1,3),diag(1,3),matrix(0,nrow = 3, ncol = 5)),
                     cbind(matrix(0,nrow = 4, ncol = 9),rep(-1,4),diag(1,4))))

set.seed(300)
theta1 = rnorm(0,sd =sqrt(0.1),n = 10) 
mu = rnorm(mean = qnorm(p = 0.05), sd = sqrt(0.1), n = 1)
theta= matrix(c(mu,theta1),nrow = 11)
sort(pnorm(X%*%theta),decreasing = T)[1:20]

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


hist(pnorm(X%*%theta),xlim = c(0,1),breaks = 8)
mean(pnorm(X%*%theta)) #
which.max(pnorm(X%*%theta)) #13-0.7397591
mustar = X[which.max(pnorm(X%*%theta)),]%*%theta
sort(pnorm(X%*%theta),decreasing = T)[1:20]
order(pnorm(X%*%theta),decreasing = T)[1:20]
#############

funpaste = function(x){return(paste('V',x,sep = ""))}


N=30
arm_opt_mult_rpm = integer(N)
for(R in 1:N){
  
  theta_prior = matrix(rnorm(n = 11,mean=rep(0,11), sd = rep(1,11)),nrow = 11)
  #w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
  
  L3 = NULL
  y_acc = numeric(120)
  n_acc = numeric(120)
  
  ### phase 1: rough estimate of theta ####
  kk = 8
  xt = matrix(rep(t(X),kk),ncol = 11,byrow = T)
  yt = matrix(rbinom(n = kk*120,size = 1,prob = pnorm(xt%*%theta)),nrow =kk*120) #exp(xt%*%theta)/(1+exp(xt%*%theta))
  ind = 1:120
  
  yyt = as.data.frame(cbind(rep(ind,kk), yt))
  names(yyt ) = c("ind","V2")
  uni_yyt = aggregate(V2~ind,data = yyt, FUN = sum) #config:yt
  
  #y_acc conversions per congiguration(<=kk)
  y_acc[uni_yyt$ind] =  y_acc[uni_yyt$ind] + uni_yyt$V2
  
  uni_ind = sort(unique(ind))
  n = as.data.frame(table(yyt$ind)) # trials 
  names(n) = c("ind", "Freq")
  
  #n_acc trials per configuration(=kk)
  n_acc[as.numeric(levels(n$ind))] =  n_acc[as.numeric(levels(n$ind))] + n$Freq 
  
  
  #calculate theta: (similar probability matching method, to make it comparable to TS)
  theta_hat = apply(sim.post(nruns = 2000,ndraws = 500,n_acc,y_acc,X,theta_prior, last_N = 500),MARGIN = 1,FUN = mean)
  #data.frame(theta_hat,theta,theta_prior)
  #sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]
  #which.max(sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]) #1
  
  ####### update g_part start (use wilcox test  #####
  
  
  comp= data.frame(first = rep(0,4),second = rep(0,4),wilcox = rep(0,4), continue = rep(0,4), 
                   row.names = c("Factor1","Factor2","Factor3","Factor4"))
  ind_no = c(1)
  for (i in 1:4){
    
    if (i==1 ){
      if (theta_hat[2]>0){
        comp[1,1] = 2
        comp[1,2] = NA
        comp[1,3] = (wilcox.test(pnorm(X[X[,2]==1,]%*%theta_hat),pnorm(X[!X[,2]==1,]%*%theta_hat),alternative = 'greater'))$p.value
      }else{
        comp[1,2] = 2
        comp[1,1] = NA
        comp[1,3] = (wilcox.test(pnorm(X[!X[,2]==1,]%*%theta_hat),pnorm(X[X[,2]==1,]%*%theta_hat),alternative = 'greater'))$p.value
        
      }
      if (comp[1,3]>0.05) ind_no = c(ind_no,2)
    }else {
      if(i ==2){
        case = 3:4
        baseline = pnorm(X[!X[,3]*X[,4]==1,]%*%theta_hat)
      }else if(i ==3){
        case = 5:7
        baseline = pnorm(X[!X[,5]*X[,6]*X[,7]==1,]%*%theta_hat)
      }else {
        case =8:11
        baseline = pnorm(X[!X[,8]*X[,9]*X[,10]*X[,11]==1,]%*%theta_hat)
        
      }
      
      if(max(theta_hat[case])<0){
        comp[i,2] = (case)[order(theta_hat[case],decreasing = T)[1]]
        comp[i,1] = NA
        comp[i,3] = (wilcox.test(baseline,pnorm(X[X[,comp[i,2]]==1,]%*%theta_hat),alternative = 'greater'))$p.value
        
      }else if (theta_hat[(case)[order(theta_hat[case],decreasing = T)[2]]]<0){
        comp[i,1] = (case)[order(theta_hat[case],decreasing = T)[1]]
        comp[i,2] = NA
        comp[i,3] = (wilcox.test(pnorm(X[X[,comp[i,1]]==1,]%*%theta_hat),baseline,alternative = 'greater'))$p.value
        
      }else{
        comp[i,1:2] =  (case)[order(theta_hat[case],decreasing = T)[1:2]]
        comp[i,3] = (wilcox.test(pnorm(X[X[,comp[i,1]]==1,]%*%theta_hat),pnorm(X[X[,comp[i,2]]==1,]%*%theta_hat),alternative = 'greater'))$p.value
        
      }
      if (comp[i,3]>0.05) ind_no = c(ind_no,case)
    }
  }
  comp[,4] = (comp[,3]>0.05)
  
  g_part = matrix(0,nrow = length(ind_no), ncol = length(theta))
  for ( i in 1:NROW(g_part)){
    g_part[i,ind_no[i]] =1
    
  }
  g_part[,1] = rep(1,length(ind_no))
  
  
  
  ## update g_part end
  
  ###### Phase 2 : Design phase I ####
  
  space = X
  #g_part <- diag(rep(1,11))
  
  design_stg1 = xt
  theta_stg1 = theta_hat
  
  
  if(length(ind_no)>1){
    d1d <- opt_design(design_stg1,theta_stg1 ,n1=1600,pp=0,space)
    
    
    ####### Phase 3: execute the experiments using the generated design ####
    
    
    # #            use the probability matching to calculate the estimate of theta
    ##   data preparation
    yt = c()
    config = c()
    nt = c()
    for (i in 1:NROW(d1d$design_exact)){
      xt_ = matrix(d1d$design_exact[i,],ncol = 11) ###
      config[i] = row.match(as.numeric(xt_),X)
      
      yt[i] = rbinom(n=1, size= d1d$weight_exact[i],prob = pnorm(xt_%*%theta))
      nt[i] = d1d$weight_exact[i]
      
    }
    
    exp = data.frame(config = config, y = yt, n = nt)
    
    #replication may occur with  d1d$design_exact need to check
    if(sum(dim(d1d$design_exact)!=dim(unique(d1d$design_exact)))>0){
      df = aggregate(x = cbind(exp$y,exp$n), by = list(config),FUN = 'sum')
      names(df) = c("config", "y","n")
      exp = df
    }
    
    
    exp1 = data.frame(config = setdiff(1:120,exp$config),y=0,n=0)
    exp = rbind(exp,exp1)
    exp = exp[order(exp$config),]
    rownames(exp) <- NULL
    #> exp[1:5,]
    #config y  n
    #1      1 8 26
    #2      2 0  0
    #3      3 0  0
    
    ## (similar probability matching method, to make it comparable to TS)
    #theta_hat_pos1 = apply(sim.post(nruns = 2000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X[,ind_no],theta_stg1[ind_no], last_N = 500),MARGIN = 1,FUN = mean)
    theta_hat_pos1 = apply(sim.post(nruns = 3000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X,theta_stg1, last_N = 500),MARGIN = 1,FUN = mean)
    
    
    theta_hat_pos = rep(0,length(theta))
    #theta_hat_pos[ind_no]= theta_hat_pos1
    theta_hat_pos[ind_no]= theta_hat_pos1[ind_no]
    theta_hat_pos[setdiff(1:11,ind_no)] = theta_hat[setdiff(1:11,ind_no)]
    #order(pnorm(X%*%theta_hat_pos),decreasing = T)[1:10]
    #sort(pnorm(X%*%theta_hat_pos),decreasing = T)[1:20]
    #var(pnorm(X%*%theta_hat_pos))
  
    
  }else{
    #when stage 1 already selected the best levels for each factor, no need for stage2
    theta_hat_pos = theta_hat
  }
  
  
  arm_opt_mult_rpm[R]=which.max(pnorm(X%*%theta_hat_pos))
  
  #if(which.max(pnorm(X%*%theta_hat_pos))!=13) browser()
  
} 

table(arm_opt_mult_rpm)
#arm_opt_mult_rpm 1200+2400
#11 13 
#2 48 
#arm_opt_mult_rpm 960+600 kk=8 new g_part 4/5
#3 11 12 13 14 31 33 
#1  4  2 20  1  1  1 

#arm_opt_mult_rpm  960+800 kk=8 new g_part 4/5
#1  3 11 13 14 31 33 
#1  1 12 76  2  1  7 
#arm_opt_mult_rpm 1080+900 kk=9 new g_part 4/5
#11 13 33 
#9 40  1
#arm_opt_mult_rpm 840+900 kk=7 new g_part 4/5
#1 11 12 13 33 
#1  5  1 16  7 
#arm_opt_mult_rpm 960+900 kk=8 new g_part 4/5
#11 12 13 33 53 
#2  1 25  1  1
#arm_opt_mult_rpm 840 + 800 kk=7 new g_part 4/5
#3 11 13 33 53 
#1  4 23  1  1 
#arm_opt_mult_rpm 840 + 800 kk=7 new g_part 7/10
#11 12 13 31 33 
#6  2 19  1  2 
#arm_opt_mult_rpm 840 + 800 kk=7 new g_part 8/10
#11 12 13 14 33 53 
#11  1 15  1  1  1