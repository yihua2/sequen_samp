
#### 
#### version fin3: multistage g_part changed for comparing only the top five config; stage 2 adjusted accordingly
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
    
  
    #######phase 3 end ###
    
    if (max(g_part[,ind_no]%*%theta_hat_pos1)>0){
      arm_opt_mult_rpm[R] =out[ which.max(g_part[,ind_no]%*%theta_hat_pos1)+1]
    }else{
      arm_opt_mult_rpm[R] = out[1]
    }
    

  #if( arm_opt_mult_rpm[R]!=13) browser()
  
} 

table(arm_opt_mult_rpm)

### these are for theta(13)  ###


#arm_opt_mult_rpm 8+800
#11 13 14 15 32 53 
#1 24  2  1  1  1 

#arm_opt_mult_rpm 8+800
#11 13 14 33 
#4 21  2  3 


#arm_opt_mult_rpm 9+800
#11 13 15 31 33 51 
#5 80  2  1 11  1 

#arm_opt_mult_rpm 10+1600(2800)
#11 13 
#2 28 


#notes: why the different performance of TS and multistage? maybe because the difference between optimal and suboptimal is big/small??
#> sort(pnorm(X%*%theta),decreasing = T)[1:5]
#[1] 0.7397591 0.6983675 0.6976299 0.6534694 0.6494925
#> sort(X%*%theta,decreasing = T)[1:5]
#[1] 0.6426029 0.5197112 0.5175960 0.3947042 0.3839507

#tried theta = theta_store1, best arm 62, 2 77
#[,1]
#[1,] -1.37289875
#[2,] -0.06152554
#[3,]  0.59810961
#[4,]  0.14918034
#[5,]  0.06869351
#[6,] -0.27484273
#[7,] -0.10052901
#[8,]  0.07692614
#[9,]  0.33139563
#[10,] -0.43408377
#[11,] -0.74351460
#> sort(pnorm(X%*%theta),decreasing = T)[1:5] #smaller probabilities, but difference between 1st and 2nd
#[1] 0.3539418 0.3313366 0.3287406 0.3068078 0.2932474


# remember to record variance



### these are for theta(2)  ###

#sort(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 0.8208992 0.7382919 0.7233833 0.6598241 0.6337775 0.6225607 0.6188620 0.6062569 0.5856221 0.5522208 0.5343096 0.5243865 0.5182656 0.5154332 0.5086900 0.5063827
#[17] 0.4955589 0.4906728 0.4782836 0.4563721
#> order(pnorm(X%*%theta),decreasing = T)[1:20]
#[1]   2  12  42   4  62  52   5  22  17  14  44  72   7   3  15 102  32  45   1  57

#arm_opt_mult_rpm 10+1600
#2  5 12 17 42 
#25  1  1  1  2 


### these are for theta(32)  ###
# sort(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 0.7038422 0.6285020 0.5828967 0.5773997 0.5585128 0.5133815 0.4950735 0.4801765 0.4759205 0.4484183 0.4479160 0.4309146 0.4289790 0.3984715 0.3863868 0.3849043
#[17] 0.3795400 0.3679630 0.3613941 0.3535014
#> order(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 32 31 35 37 22 92 36 27 21 52 40 91 25 26 34 95 97 51 82 30



### these are for theta(27) ##
#> sort(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 0.6798228 0.6694615 0.6218070 0.6108186 0.5406304 0.5359089 0.5274656 0.5244620 0.5134859 0.5064892 0.5020082 0.4950094 0.4780887 0.4649094 0.4390600 0.4339999
#[17] 0.4277405 0.4227057 0.3916393 0.3789762
#> order(pnorm(X%*%theta),decreasing = T)[1:20]
#[1] 27 22 28 23 32 26 37 21 29 30 24 25 33 38 47 87 42 82 31 36

#arm_opt_mult_rpm #10+1600
#22 23 25 27 28 37 
#11  2  1 14  1  1 