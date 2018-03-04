## Bernoulli

#list of configurations
N= matrix(c(rep(0,60),rep(1,60), 
            rep(0,20),rep(1,20),rep(2,20),rep(0,20),rep(1,20),rep(2,20),
            rep(c(rep(0,5),rep(1,5),rep(2,5),rep(3,5)),6),
            rep(c(0,1,2,3,4), 24)), nrow= 120, ncol = 4)


X = matrix(c(rep(1,120),
             rep(1,60),rep(0,60),
             rep(0,60),rep(1,60),
             rep(1,20),rep(0,40),rep(1,20),rep(0,40),
             rep(0,20),rep(1,20), rep(0,20),rep(0,20),rep(1,20), rep(0,20),
             rep(0,40),rep(1,20),rep(0,40),rep(1,20),
             rep(c(rep(1,5),rep(0,15)),6),
             rep(c(rep(0,5),rep(1,5),rep(0,10)),6),
             rep(c(rep(0,10),rep(1,5),rep(0,5)),6),
             rep(c(rep(0,15),rep(1,5)),6),
             rep(c(1,0,0,0,0),24),
             rep(c(0,1,0,0,0),24),
             rep(c(0,0,1,0,0),24),
             rep(c(0,0,0,1,0),24),
             rep(c(0,0,0,0,1),24)), nrow=120, ncol=15)


########## Bernoulli Bandit For one single simulation  #########
set.seed(300)
theta1 = rnorm(0,sd =sqrt(0.1),n = 14)
mu = rnorm(mean = qnorm(p = 0.05), sd = sqrt(0.1), n = 1)
theta= matrix(c(mu,theta1),nrow = 15)
hist(pnorm(X%*%theta),xlim = c(0,0.2),breaks = 11)
mean(pnorm(X%*%theta))
which.max(pnorm(X%*%theta))
mustar = X[which.max(pnorm(X%*%theta)),]%*%theta


#   w/ replacement(makes more sense as the regret goes down faster)
w0 = matrix(c(rep(1,240)),nrow = 120) #parameters of beta distribution(alpha, beta).(1,1) are the priors

L1 = NULL
weight = NULL
for (t in 1:200){
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
}
plot(L1,ylim = c(0,60),type = 'l')  


###### Fractional Factorial for one simulation  #####

library(truncnorm)
library(MASS)

theta_prior = matrix(rnorm(n = 15,mean=rep(0,15), sd = rep(1,15)),nrow = 15)
w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
L2 = NULL
y_acc = numeric(120)
n_acc = numeric(120)

for (t in 1:100){
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
}


plot(L2,ylim = c(0,60),type='l')




