
#######################
#### version fin1: multistage change stage 2 only estimate part of the parameters version 1 ####
########################

####code _functions changed ####
sample_post = function(n_acc,y_acc,X1,theta_pri ,k0){
  z=c()
  
  for (k in 1:dim(X1)[1]){
    mean0 = as.numeric(matrix(X1[k,],nrow = 1)%*%theta_pri) ##
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
N=100
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
  
  k0 = 120
  #calculate theta: (similar probability matching method, to make it comparable to TS)
  theta_hat = apply(sim.post(nruns = 2000,ndraws = 500,n_acc,y_acc,X,theta_prior, last_N = 500,k0),MARGIN = 1,FUN = mean) ###
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
X1 = unique(X[,ind_no])

yt = c()
config = c()
nt = c()
for (i in 1:NROW(d1d$design_exact)){
  xt_ = matrix(d1d$design_exact[i,],ncol = 11) ###
  #xt_ = matrix(d1d$design_exact[i,ind_no],ncol = length(ind_no)) ###
  
  config[i] = row.match(as.numeric(xt_[ind_no]),X1)
  
  #yt[i] = rbinom(n=1, size= d1d$weight_exact[i],prob = pnorm(xt_%*%theta[ind_no]))
  yt[i] = rbinom(n=1, size= d1d$weight_exact[i],prob = pnorm(xt_%*%theta))
  nt[i] = d1d$weight_exact[i]
  
}

exp = data.frame(config = config, y = yt, n = nt)

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
k0 = ifelse(comp[1,]$continue==T,2,1)*ifelse(comp[2,]$continue==T,3,1)*ifelse(comp[3,]$continue==T,4,1)*ifelse(comp[4,]$continue==T,5,1)
if (NROW(exp)<k0){
  exp1 = data.frame(config = setdiff(1:k0,exp$config),y=0,n=0)
  exp = rbind(exp,exp1)
  
}
exp = exp[order(exp$config),]
rownames(exp) <- NULL


#theta_pri = c(theta_stg1[1]+)
## (similar probability matching method, to make it comparable to TS)
#theta_hat_pos1 = apply(sim.post(nruns = 2000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X,theta_stg1, last_N = 500),MARGIN = 1,FUN = mean)
theta_hat_pos1 = apply(sim.post(nruns = 3000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X1,theta_stg1[ind_no], last_N = 500,k0 ),MARGIN = 1,FUN = mean)


#######phase 3 end ###



theta_hat_pos = rep(0,length(theta))
#theta_hat_pos[ind_no]= theta_hat_pos1
theta_hat_pos[ind_no]= theta_hat_pos1
theta_hat_pos[setdiff(1:length(theta),ind_no)] = theta_hat[setdiff(1:length(theta),ind_no)]
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

#arm_opt_mult_rpm 8+1600 
#11 13 31 33 
#9 88  1  2

#arm_opt_mult_rpm 8+800
#arm_opt_mult_rpm
#11 13 31 33 
#13 69  5 13 