# This is a multi-replicate version of code_22(multistage) and code_21(binomial bandit/fractional factorial)
# Goal: to see for three methods in the final convergence , 
#         1. the number of times the method will finally choose the actual best arm
#         2. the loss function


# Shared set up: model structure probit model with 120 config(4 factors, each of 2,3,4, and 5 levels)
#                the same set of generated model parameters(#22 being the best)
####           Data preparation        ####
set.seed(300)
library(truncnorm)
library(MASS)


N= matrix(c(rep(0,60),rep(1,60), 
            rep(0,20),rep(1,20),rep(2,20),rep(0,20),rep(1,20),rep(2,20),
            rep(c(rep(0,5),rep(1,5),rep(2,5),rep(3,5)),6),
            rep(c(0,1,2,3,4), 24)), nrow= 120, ncol = 4)


X = matrix(c(rep(1,120),#int
             rep(1,60),rep(0,60),#factor1: 1
             rep(0,60),rep(1,60),#factor1: 2
             rep(1,20),rep(0,40),rep(1,20),rep(0,40),#factor2: 1
             rep(0,20),rep(1,20), rep(0,20),rep(0,20),rep(1,20), rep(0,20),#factor2: 2
             rep(0,40),rep(1,20),rep(0,40),rep(1,20),#factor2: 3
             rep(c(rep(1,5),rep(0,15)),6),#factor3: 1
             rep(c(rep(0,5),rep(1,5),rep(0,10)),6),#factor3: 2
             rep(c(rep(0,10),rep(1,5),rep(0,5)),6),#factor3: 3
             rep(c(rep(0,15),rep(1,5)),6),#factor3: 4
             rep(c(1,0,0,0,0),24),#factor4: 1
             rep(c(0,1,0,0,0),24),#factor4: 2
             rep(c(0,0,1,0,0),24),#factor4: 3
             rep(c(0,0,0,1,0),24),#factor4: 4
             rep(c(0,0,0,0,1),24)), nrow=120, ncol=15)#factor4: 5


set.seed(300)
theta1 = rnorm(0,sd =sqrt(0.1),n = 14)
set.seed(300)
mu = rnorm(mean = qnorm(p = 0.05), sd = sqrt(0.1), n = 1)
theta= matrix(c(mu,theta1),nrow = 15)
#hist(pnorm(X%*%theta),xlim = c(0,0.2),breaks = 11)
#mean(pnorm(X%*%theta)) #
#which.max(pnorm(X%*%theta)) #22-0.6443
mustar = X[which.max(pnorm(X%*%theta)),]%*%theta
#sort(pnorm(X%*%theta),decreasing = T)[1:20]



#####       Binomial Bandit       #####
# Repeat 100 times. Each time stop when one arm is getting over 50 spots in the sample
N=100
arm_opt = c()
iter = c()
for (R in 1:N){
  w0 = matrix(c(rep(1,240)),nrow = 120) #parameters of beta distribution(alpha, beta).(1,1) are the priors
  L1 = NULL
  weight = NULL
  proceed = True
  t=1
  while(proceed){
    # sample weights from updated beta distribution
    for(i in 1:120){
      weight[i] = rbeta(n=1, shape1 = w0[i,1], shape2 = w0[i,2])
    }
    # Sample 100 arms
    if(t>1){
      y = rep(0,120)
      y[uni_yyt$ind] = uni_yyt$V2
      nn = rep(0,120)
      nn[as.numeric(levels(n$ind))] = n$Freq
      
      weight1 = compute.probopt(y=y,n=nn,param=w0)
    }else weight1 = weight
    
    
    ind = sample(x = c(1:120), prob = weight1,size = 100,replace = T) #w0[,1]/rowSums(w0)
    
    xt = matrix(X[ind,],nrow = 100)
    yt = matrix(rbinom(n = 100,size = 1,prob = exp(xt%*%theta)/(1+exp(xt%*%theta))),nrow = 100)
    yyt = as.data.frame(cbind(ind, yt))
    uni_yyt = aggregate(V2~ind,data = yyt, FUN = sum) # number of '1' response for each chosen config
    uni_ind = sort(unique(ind)) 
    
    n = as.data.frame(table(ind)) # times of being selected for chosen config
    # update
    w0[uni_ind,] = w0[uni_ind,] + cbind(uni_yyt$V2,n$Freq-uni_yyt$V2) 
    # Expected regret
    L1[t] = sum(matrix(n$Freq,nrow = 1)%*%
                  (expit(matrix(as.numeric(mustar),nrow = length(uni_ind)))-expit((X[uni_ind,]%*%theta))))
    #browser()#Q to quit
    proceed = (max(table(ind))<70)
    t = t + 1
  }
  #plot(L1,ylim = c(0,60),type = 'l')  
  arm_opt[R] = as.numeric(levels(n[which.max(n$Freq),]$ind)[n[which.max(n$Freq),]$ind])
  iter[R] = t
}  

#> arm_opt
#[1]  22  22   2  22  42   2  22  22  62  22  22  22  82  22  22  22  22  22  22  22  22  42   2  22  22
#[26]  22  22  22  22  22  22  22  22   2  22  22  22  22  22  22  22  22  22  22  22  22   2  22  22  22
#[51]  22  22  97  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22  22   2
#[76]  22  22  42  67  22  22  22 102  22  22  22   2  22  22   2  22   2   2  22  62  22  22  22  97  37
#> iter
#[1] 107  94 121  44  17 182 140 107  49  74 148  58  81  98 158  75 162  53 208  75  99  83  70 136 184
#[26] 171 155 142 100  89 180  88 194  47 224 152  59 128 229  67 314  88  68  85  89 154 103  63 123  89
#[51]  91  90  42 108 242  82 116 109  17 198 163  55 134 135 117 175 107  77 154 118 193 103 147 209  56
#[76] 110 164  31  22  75 228 125  31  80 162 114 163  75 237 155 111  91 117 157  72  86  20 164  76  29
  
#####       Fractional Binomial Bandit   #####

library(truncnorm)
library(MASS)
N=100
arm_opt_frac = c()
iter_frac = c()
for (R in 1:N){
  theta_prior = matrix(rnorm(n = 15,mean=rep(0,15), sd = rep(1,15)),nrow = 15)
  w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
  L2 = NULL
  y_acc = numeric(120)
  n_acc = numeric(120)
  proceed = True
  t=1
  while(proceed){
    # Sample 100 arms according to conjugacy priors
    ind = sample(x = c(1:120), prob = w1,size = 100,replace = T)#replacement?
    
    
    xt = matrix(X[ind,],nrow = 100)
    yt = matrix(rbinom(n = 100,size = 1,prob = pnorm(xt%*%theta)),nrow = 100) #exp(xt%*%theta)/(1+exp(xt%*%theta))
    
    yyt = as.data.frame(cbind(ind, yt))
    uni_yyt = aggregate(V2~ind,data = yyt, FUN = sum)# conversions 
    
    y_acc[uni_yyt$ind] =  y_acc[uni_yyt$ind] + uni_yyt$V2
    
    uni_ind = sort(unique(ind))
    
    n = as.data.frame(table(ind)) # trials 
    n_acc[as.numeric(levels(n$ind))] =  n_acc[as.numeric(levels(n$ind))] + n$Freq
    
    # Expected regret
    L2[t] = sum(matrix(n$Freq,nrow = 1)%*%
                  (pnorm(matrix(as.numeric(mustar),nrow = length(uni_ind)))-pnorm((X[uni_ind,]%*%theta))))
    
    
    
    #update w1
    
    # ppst = sim.post(nruns = 1000,ndraws= 50,X=X,xt=xt,uni_yyt = uni_yyt,n = n, theta_prior = theta_prior, last_N = 100)
    
    
    #mu_mat = (X%*%t(post1))
    
    #  tbl = as.data.frame(table(apply(mu_mat,MARGIN = 2,which.max)))# var1 . freq
    # tbl$weight = tbl$Freq/50
    #  w1 = rep(0,120)
    #  w1[as.numeric(levels(tbl$Var1))] = tbl$weight
    
    w10 =compute.win.prob(nruns = 1000,ndraws= 500,X=X,y_acc = y_acc,n_acc = n_acc, 
                          theta_prior = theta_prior, last_N = 500)
    #w1 = w1 + as.data.frame(w10)$Freq
    w1 = as.data.frame(w10)$Freq
  
    proceed = (max(table(ind))<70)
    t = t + 1
}
  arm_opt_frac[R] = as.numeric(levels(n[which.max(n$Freq),]$ind)[n[which.max(n$Freq),]$ind])
  iter_frac[R] = t
  
  #plot(L2,ylim = c(0,60),type='l')
  
}


##### Multistage fractional factorial #####
N=20
arm_opt_mult_glm = c()
arm_opt_mult_rpm =c()
for (R in 1:N){
  theta_prior = matrix(rnorm(n = 15,mean=rep(0,15), sd = rep(1,15)),nrow = 15)
  w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
  L3 = NULL
  y_acc = numeric(120)
  n_acc = numeric(120)
  
  
  ### phase 1: rough estimate of theta 
  
  xt = X
  yt = matrix(rbinom(n = 120,size = 1,prob = pnorm(xt%*%theta)),nrow = 120) #exp(xt%*%theta)/(1+exp(xt%*%theta))
  ind = 1:120
  
  yyt = as.data.frame(cbind(ind, yt))
  uni_yyt = aggregate(V2~ind,data = yyt, FUN = sum)# conversions 
  
  y_acc[uni_yyt$ind] =  y_acc[uni_yyt$ind] + uni_yyt$V2
  
  uni_ind = sort(unique(ind))
  
  n = as.data.frame(table(ind)) # trials 
  n_acc[as.numeric(levels(n$ind))] =  n_acc[as.numeric(levels(n$ind))] + n$Freq
  
  
  
  
  #calculate theta: (similar probability matching method, to make it comparable to TS)
  theta_hat = apply(sim.post(nruns = 1000,ndraws = 500,n_acc,y_acc,X,theta_prior, last_N = 500),MARGIN = 1,FUN = mean)
  
  #sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]
  #which.max(sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]) #1
  
  
  ###### Phase 2 : Design phase
  
  space = X
  g_part <- diag(rep(1,15))
  design_stg1 = X
  theta_stg1 = theta_hat
  
  invisible(d1d <- opt_design(design_stg1,theta_stg1 ,n=500,pp=0,space))
  
  #d1d$design_exact
  #d1d$weight-exact
  
  ####### Phase 3: execute the experiments using the generated design

  # #            use the probability matching to calculate the estimate of theta
  ##   data preparation
  yt = c()
  config = c()
  nt = c()
  for (i in 1:NROW(d1d$design_exact)){
    xt = matrix(d1d$design_exact[i,],ncol = 15)
    config[i] = prodlim::row.match(as.numeric(xt),X)
    yt[i] = rbinom(n=1, size= d1d$weight_exact[i],prob = pnorm(xt%*%theta))
    nt[i] = d1d$weight_exact[i]
  }
  exp = data.frame(config = config, y = yt, n = nt)
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
  theta_hat_pos = apply(sim.post(nruns = 1000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X,theta_stg1, last_N = 500),MARGIN = 1,FUN = mean)
  arm_opt_mult_rpm[R] = which.max(pnorm(X%*%theta_hat_pos))
  #[1] 22
  
  
  
  # #          ....or, use glm() to estimate theta
  ##   data preparation
  datX = data.frame(x1 = c(rep(1,60),rep(2,60)),
                    x2 = rep(c(rep(1,20),rep(2,20),rep(3,20)),2),
                    x3 = rep(c(rep(1,5),rep(2,5),rep(3,5),rep(4,5)),6),
                    x4 = rep(c(1:5),24))
  
  exp[exp$n!=0,]$config
  dat = as.data.frame(matrix(nrow = 1,ncol = 6))
  tmp = c()
  row=1
  for (i in 1:120){
    if (i%in% exp[exp$n!=0,]$config){
      tot = exp[i,]$y
      for(j in 1:exp[i,]$n){
        
        dat[row,] = c(as.numeric(tot>0),datX[i,],exp$config[i])
        row= row +1
        tot = max(tot-1,0)
      }
    }
  }
  names(dat) = c("y","x1","x2","x3","x4","config")
  
  fit = glm(y~factor(x1)+factor(x2)+factor(x3)+factor(x4),data = dat, family = binomial(link = "probit"))
  
  arm_opt_mult_glm[R] = dat[which.max(fitted.values(fit)),]$config # 77??
  
  
  
  
  # Expected regret for the experiment stage
  
  #L3= sum(matrix(exp[exp$n!=0,]$n,nrow = 1)%*%
  #          (pnorm(matrix(as.numeric(mustar),nrow =NROW(d1d$design_exact) ))-pnorm((X[exp[exp$n!=0,]$config,]%*%theta_hat_pos))))
  
}

