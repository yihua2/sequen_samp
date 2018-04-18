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
