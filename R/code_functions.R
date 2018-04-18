####### expit #####

expit<-function(x){
  exp(x)/(1+exp(x))
}

#######       compute sampling weights using quadrature --binomial bandit  ##########
compute.probopt = function(y,n,param){
  k = length(y)
  ans = numeric(k)
  for ( i in 1:k){
    indx = (1:k)[-i]
    f = function(x){
      r = dbeta(x,y[i]+param[i,1],n[i]-y[i]+param[i,2])
      for (j in indx) r = r*pbeta(x,y[j]+param[j,1],n[j]-y[j]+param[j,2])
      return(r)
    }
    ans[i] = integrate(f,0,1)$value
  }
  return(ans)
}

########        compute sampling weights using simulation --binomial bandit   ##########
sim.post<- function(y,n,ndraws){
  k <- length(y)
  ans <- matrix(nrow=ndraws, ncol = k)
  no=  n-y
  for (i in 1:k){
    ans[,i] = rbeta(ndraws,y[i]+1,no[i]+1)
    
  }
  return(ans)
}

prob.winner = function(post){
  k = ncol(y)
  w = table(factor(max.col(post),level = 1:k))
  return(w/sum(w))
}

compute.win.prob = function(y,n,ndraws){
  return(prob.winner(sim.post(y,n,ndraws)))
  
}

########       compute sampling weights using simulation --fractional factorial    ############

#calculate poster distribution of theta (the standard way) _______not modified(don't use)
#sample_post= function(xt,yt,theta_prior){
#  z=c()
#  for(k in 1:100){
#    mean0 = as.numeric(matrix(xt[k,],nrow = 1)%*%theta_prior)   
#    if (yt[k]>0) {
#      z[k] = rtruncnorm(n=1, a=0, mean=mean0, sd=1)    
#    }
#    if (yt[k]==0) {
#      z[k] = rtruncnorm(n=1, b=0, mean=mean0, sd=1)
#    }
#  }
#  omg = solve(diag(rep(1,15))+t(xt)%*%xt)
#  theta_tilde = omg%*% (t(xt)%*%matrix(z,nrow = 100))
  
#  samp = mvrnorm(n = 1, mu = theta_tilde, Sigma = omg) 
#  return(samp)
#} 


#the indicator way to calculate posterior distribution theta #
sample_post = function(n_acc,y_acc,X,theta_pri){
  z=c()
  
  for (k in 1:dim(X)[1]){
    mean0 = as.numeric(matrix(X[k,],nrow = 1)%*%theta_pri)
    if(y_acc[k]==0 || n_acc[k]-y_acc[k] == 0){
      z[k] = 0#0.5
    }else{
      z1 = sum(rtruncnorm(n=y_acc[k], sd=1,a=0, mean = mean0))#a=0
      z2 = sum(rtruncnorm(n=n_acc[k]-y_acc[k], sd=1,b=0,mean = mean0))#b=0.5
      z[k] = z1+z2
    }
   
  }
 # z = rep(0,120)
#  z[uni_yyt$ind] = z0
  
 # n_hat = rep(0,120)
#  n_hat[as.numeric(levels(n$ind))] = n$Freq
  p = length(theta_pri)
  omg = solve(diag(rep(1,p))+t(X)%*%diag(n_acc)%*%X)
  theta_tilde = omg%*%(t(X)%*%matrix(z,nrow = 120))
  
  samp = mvrnorm(n = 1, mu = theta_tilde, Sigma = omg) 
  return(samp)
}



sim.post<- function(nruns,ndraws,n_acc,y_acc,X,theta_prior, last_N){
  
  # generate Markov chain of posterior distributions:nruns
  p = length(theta_prior)
  post = matrix(NA,nrow = nruns, ncol =p)
  post[1,] = sample_post(n_acc,y_acc,X,theta_pri = theta_prior)
  for (kk in 2:nruns){
    post [kk,] = sample_post(n_acc,y_acc,X,theta_pri = post[kk-1,])
  }
  
  #sample ndraws theta from the (stabilzed) posterior distribution(assume the last last_N is stable)
  ind1 = sample(c((nruns-last_N+1):nruns),size = ndraws, replace = T)
  post1 = post[ind1,]
  
  return(t(post1))
}



prob.winner = function(post){
  k = dim(X)[1]
  w = table(factor(max.col(t(post)%*%t(X)),level = 1:k))
  return(w/sum(w))
}

compute.win.prob = function(nruns,ndraws,n_acc,y_acc,X,theta_prior, last_N){
  return(prob.winner(sim.post(nruns,ndraws,n_acc,y_acc,X,theta_prior, last_N)))
  
}