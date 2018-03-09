#######################
# version 1: multistage
########################
source("/Users/Vera/Dropbox/study/google/owea interference-copy.R")
set.seed(300)
library(truncnorm)
library(MASS)
library(prodlim)


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


########## Bernoulli Bandit For one single simulation  #########
set.seed(300)
theta1 = rnorm(0,sd =sqrt(0.1),n = 14)
mu = rnorm(mean = qnorm(p = 0.05), sd = sqrt(0.1), n = 1)
theta= matrix(c(mu,theta1),nrow = 15)
hist(pnorm(X%*%theta),xlim = c(0,0.2),breaks = 11)
mean(pnorm(X%*%theta)) #
which.max(pnorm(X%*%theta)) #22-0.6443
mustar = X[which.max(pnorm(X%*%theta)),]%*%theta
sort(pnorm(X%*%theta),decreasing = T)[1:20]


###### Fractional Factorial with multistage #####

N=100
arm_opt_mult_glm = c()
arm_opt_mult_rpm =c()
for (R in 1:N){
  
theta_prior = matrix(rnorm(n = 15,mean=rep(0,15), sd = rep(1,15)),nrow = 15)
w1 = pnorm(X%*%theta_prior)#matrix(c(rep(1/2,120)),nrow = 1)
L3 = NULL
y_acc = numeric(120)
n_acc = numeric(120)


### phase 1: rough estimate of theta 
kk = 50
xt = matrix(rep(t(X),kk),ncol = 15,byrow = T)
yt = matrix(rbinom(n = kk*120,size = 1,prob = pnorm(xt%*%theta)),nrow =kk*120) #exp(xt%*%theta)/(1+exp(xt%*%theta))
ind = 1:120

yyt = as.data.frame(cbind(rep(ind,kk), yt))
names(yyt ) = c("ind","V2")
uni_yyt = aggregate(V2~ind,data = yyt, FUN = sum)# conversions 

y_acc[uni_yyt$ind] =  y_acc[uni_yyt$ind] + uni_yyt$V2

uni_ind = sort(unique(ind))

n = as.data.frame(table(yyt$ind)) # trials 
names(n) = c("ind", "Freq")
n_acc[as.numeric(levels(n$ind))] =  n_acc[as.numeric(levels(n$ind))] + n$Freq




#calculate theta: (similar probability matching method, to make it comparable to TS)
theta_hat = apply(sim.post(nruns = 1000,ndraws = 500,n_acc,y_acc,X,theta_prior, last_N = 500),MARGIN = 1,FUN = mean)

#sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]
#which.max(sort(pnorm(X%*%theta_hat),decreasing = T)[1:20]) #1


###### Phase 2 : Design phase

space = X
g_part <- diag(rep(1,15))
design_stg1 = xt
theta_stg1 = theta_hat

d1d <- opt_design(design_stg1,theta_stg1 ,n=6000,pp=0,space)

#d1d$design_exact
#d1d$weight_exact

####### Phase 3: execute the experiments using the generated design


# #            use the probability matching to calculate the estimate of theta
##   data preparation
yt = c()
config = c()
nt = c()
for (i in 1:NROW(d1d$design_exact)){
  xt = matrix(d1d$design_exact[i,],ncol = 15)
  config[i] = row.match(as.numeric(xt),X)
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

#sort(pnorm(X%*%theta_hat_pos),decreasing = T)[1:10] # estimated ranking
#[1] 0.6010915 0.5690317 0.5610377 0.5521877 0.5340979 0.5284351 0.5195111 0.5114153 0.4932192 0.4842826 
#c(1:120)[order(pnorm(X%*%theta_hat_pos),decreasing = T)][1:20]
#[1]  22  37   2  82  27  17  97  62   7  87  42  77  32  23  57  67  25  12 102  38

#sort(pnorm(X%*%theta),decreasing = T)[1:10] # true ranking
#[1] 0.6443884 0.6172192 0.6016559 0.5825437 0.5736236 0.5542383 0.5527565 0.5483755 0.5381710 0.5241612
#c(1:120)[order(pnorm(X%*%theta),decreasing = T)][1:20]
#[1]  22   2  37  82  17  62  27  42  97   7  77  57  32  87 102  12  67  47 117  92




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
#168.1667 it's big..
# to conteract more obs in this experiments(n=500), 168.1667/5 = 33.63334 is not big


#total # of observations used to get best arm: 120+500 = 620
#total # of observations used to get best arm for fractional factorial: 16*100 = 1600
#total # of observations used to get best arm for binomial bandit: 80*100 = 8000
}

table(arm_opt_mult_rpm)
#arm_opt_mult_rpm  4800+4800
#1  2  7 17 21 22 25 27 32 37 42 97 
#1 22  1  1  2 60  1  2  1  7  1  1 
#arm_opt_mult_rpm  2400+2400
#2  17  21  22  27  32  37  45  52  57  62  67  71  82  87  97 102 
#16   3   1  48   3   2  15   1   1   2   1   2   1   1   1   1   1 


save.image("/Users/Vera/Documents/GitHub/sequen_samp/R/code_22_copyimg.RData")
