
#### 
#### version fin4: multistage g_part changed for comparing only the top five config; stage 2 adjusted accordingly; with weights
##  ##

####code _functions changed ####
sample_post = function(n_acc,y_acc,X1,theta_pri ,k0){
  z=c()
  
  for (k in 1:dim(X1)[1]){
    mean0 = as.numeric(matrix(X1[k,],nrow = 1)%*%theta_pri)+mu0[k] ##
    if(y_acc[k]==0 || n_acc[k]-y_acc[k] == 0){
      z[k] = 0#0.5
    }else{
      z1 = sum(rtruncnorm(n=y_acc[k], sd=1,a=0, mean = mean0))#a=0
      z2 = sum(rtruncnorm(n=n_acc[k]-y_acc[k], sd=1,b=0,mean = mean0))#b=0.5
      z[k] = z1+z2
    }
    
  }
  
  p = length(theta_pri)
  
  omg = solve(diag(rep(1,p))+t(X1)%*%diag(n_acc)%*%X1) 
  
  #theta_tilde = omg%*%(t(X1)%*%matriX1(z,nrow = 120))  ###
  theta_tilde = omg%*%(t(X1)%*%matrix(z,nrow = k0))
  
  samp = mvrnorm(n = 1, mu = theta_tilde, Sigma = omg) 
  return(samp)
}



sim.post<- function(nruns,ndraws,n_acc,y_acc,X1,theta_prior, last_N ,k0){
  
  # generate Markov chain of posterior distributions
  
  p = length(theta_prior)
  post = matrix(NA,nrow = nruns, ncol =p)###
  post[1,] = sample_post(n_acc,y_acc,X1,theta_pri = theta_prior ,k0)
  for (kk in 2:nruns){
    post [kk,] = sample_post(n_acc,y_acc,X1,theta_pri = post[kk-1,] ,k0)
  }
  
  #sample 100 theta from the (stabilzed) posterior distribution(assume the last 70 is stable)
  ind1 = sample(c((nruns-last_N+1):nruns),size = ndraws, replace = T)
  post1 = post[ind1,]
  
  return(t(post1))
}



prob.winner = function(post){
  k = dim(X1)[1]
  w = table(factor(max.col(t(post)%*%t(X1)),level = 1:k))
  return(w/sum(w))
}

compute.win.prob = function(nruns,ndraws,n_acc,y_acc,X1,theta_prior, last_N ,k0){
  return(prob.winner(sim.post(nruns,ndraws,n_acc,y_acc,X1,theta_prior, last_N ,k0)))
  
}

##### code_fucntions end


#### the iterations ####

N=30
arm_opt_mult_rpm = integer(N)
n1 = integer(N)
for(R in 1:N){
  
  theta_prior = matrix(rnorm(n = 11,mean=rep(0,11), sd = rep(1,11)),nrow = 11)
  #w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
  
  L3 = NULL
  y_acc = numeric(120)
  n_acc = numeric(120)
  
  ### phase 1: rough estimate of theta ####
  kk = 10
  xt = matrix(rep(t(X),kk),ncol = 11,byrow = T)
  yt = matrix(rbinom(n = kk*120,size = 1,prob = pnorm(xt%*%theta)),nrow =kk*120) #exp(xt%*%theta)/(1+exp(xt%*%theta))
  ind = 1:120
  
  yyt = as.data.frame(cbind(rep(ind,kk), yt))
  names(yyt ) = c("ind","V2")
  uni_yyt = aggregate(V2~ind,data = yyt, FUN = sum) #config:yt
  
  #y_acc conversions per congiguration(<=kk)
  y_acc[uni_yyt$ind] =  y_acc[uni_yyt$ind] + uni_yyt$V2
  
  #uni_ind = sort(unique(ind))
  n = as.data.frame(table(yyt$ind)) # trials 
  names(n) = c("ind", "Freq")
  
  #n_acc trials per configuration(=kk)
  n_acc[as.numeric(levels(n$ind))] =  n_acc[as.numeric(levels(n$ind))] + n$Freq 
  
  k0 = 120
  mu0 = rep(0,k0)
  
  
  #calculate theta: (similar probability matching method, to make it comparable to TS)
  theta_hat = apply(sim.post(nruns = 2000,ndraws = 500,n_acc,y_acc,X,theta_prior, last_N = 500,k0),MARGIN = 1,FUN = mean) ###
  #data.frame(theta_hat,theta,theta_prior)
  #sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]
  #which.max(sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]) #1
  
  ####### update g_part start - to compare top five effectively  #####
  
  
  # out = order(pnorm(X%*%theta_hat),decreasing = T)[1:5]
  out = order(y_acc, decreasing = T)[1:5]
  
  
  temp = X[out,]
  g_part = matrix(0, nrow = 4,ncol = 11)
  for (i in 2:NROW(temp)){
    g_part[i-1,] = temp[i,]-temp[1,]
  }
  
  ind_no = which(apply(g_part==0,MARGIN = 2,FUN = sum)<4) # find parameters that are not needed in further estimation
  
  ## update g_part end
  
  ###### Phase 2 : Design phase I ####
  
  space = X
  #g_part <- diag(rep(1,11))
  
  design_stg1 = xt
  theta_stg1 = theta_hat
  
  #total number of levels k0
  #k0 = ifelse(comp[1,]$continue==T,2,1)*ifelse(comp[2,]$continue==T,3,1)*ifelse(comp[3,]$continue==T,4,1)*ifelse(comp[4,]$continue==T,5,1)
  #n1[R] = floor(length(ind_no)*1200/11)
  d1d <- opt_design(design_stg1,theta_stg1 ,n1 = 1600,pp=0,space)
  
  
  ####### Phase 3: execute the experiments using the generated design ####
  
  
  # #            use the probability matching to calculate the estimate of theta
  ##   data preparation
  X1 = unique(X[,ind_no])
  #X1= d1d$design_exact[]
  
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
  
  #replication may occur with  d1d$design_exact need to check
  
  #if(sum(dim(d1d$design_exact)!=dim(unique(d1d$design_exact)))>0){
  #  df = aggregate(x = cbind(exp$y,exp$n), by = list(config),FUN = 'sum')
  #  names(df) = c("config", "y","n")
  #  exp = df
  #}
  
  if(length(exp$config)!=length(unique(exp$config))){
    df = aggregate(x = cbind(exp$y,exp$n), by = list(exp$config),FUN = 'sum')
    names(df) = c("config", "y","n")
    exp = df
  }
  
  
  #total number of levels k0
  #k0 = ifelse(comp[1,]$continue==T,2,1)*ifelse(comp[2,]$continue==T,3,1)*ifelse(comp[3,]$continue==T,4,1)*ifelse(comp[4,]$continue==T,5,1)
  k0 = NROW(X1)
  if (NROW(exp)<k0){
    exp1 = data.frame(config = setdiff(1:k0,exp$config),y=0,n=0)
    exp = rbind(exp,exp1)
    
  }
  exp = exp[order(exp$config),]
  rownames(exp) <- NULL
  # mu0 = exp$mu0
  mu0 = rep(0,k0)
  
  #theta_pri = c(theta_stg1[1]+)
  ## (similar probability matching method, to make it comparable to TS)
  #theta_hat_pos1 = apply(sim.post(nruns = 2000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X,theta_stg1, last_N = 500),MARGIN = 1,FUN = mean)
  theta_hat_pos1 = apply(sim.post(nruns = 3000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X1,theta_stg1[ind_no], last_N = 500,k0 ),MARGIN = 1,FUN = mean)
  
  post = sim.post(nruns = 3000,ndraws = 1000,n_acc =exp$n ,y_acc = exp$y,X1,theta_stg1[ind_no], last_N = 500,k0 )
  pri = matrix(rep(theta_stg1,dim(post)[2]),nrow = dim(post)[2],byrow = T)
  pri[,ind_no] = t(post)
  
  ## this
  k = 120
  w = table(factor(max.col(pri%*%t(X)),level = 1:k))
  w1 = (w/sum(w))
  #which.max(w1)
 
  ## or this
 
 # w = table(factor(max.col(t(post)%*%t(X1)),level = 1:k0))
 # row.match(as.numeric(X1[ which.max(w1),]),X[out,ind_no])
  
  
  
  #######phase 3 end ###
  
  arm_opt_mult_rpm[R] = which.max(w1)
  
  #if( arm_opt_mult_rpm[R]!=13) browser()
  
} 

table(arm_opt_mult_rpm)

### these are for theta(13)  ###

#sort(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 0.7397591 0.6983675 0.6976299 0.6534694 0.6494925 0.6448303 0.6443016 0.6436732 0.6029765 0.6022263 0.5981295 0.5969277 0.5967617 0.5961087 0.5824527 0.5582508 0.5541772
#[18] 0.5479152 0.5476154 0.5443221
#> order(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 13 33 11 31 12 14 53  3 32 15 34 23 51  1 73 18 35 21  8 52

#arm_opt_mult_rpm 10+1600
#1  3  8 11 12 13 14 18 28 31 33 34 38 53 
#1  2  1  1  6  4  3  3  1  1  1  1  1  4 


### these are for theta(26)  ###