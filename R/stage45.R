###### Phase 4 : Design phase II #####
#combine stg1 and the previous design\
design_stg2 = xt
for (i in 1:NROW(d1d$design_exact)){
  design_stg2 = rbind(design_stg2,matrix(rep(d1d$design_exact[i,],d1d$weight_exact[i]),ncol = 11,byrow = T))
}

d1d2 <- opt_design(design_stg2,theta_hat_pos ,n1=1200,pp=0,space)

####### Phase 5: execute the experiments using the generated design for Phase 4 ####


# #            use the probability matching to calculate the estimate of theta
##   data preparation
yt = c()
config = c()
nt = c()
for (i in 1:NROW(d1d2$design_exact)){
  xt_ = matrix(d1d2$design_exact[i,],ncol = 11)
  config[i] = row.match(as.numeric(xt_),X)
  
  yt[i] = rbinom(n=1, size= d1d2$weight_exact[i],prob = pnorm(xt_%*%theta))
  nt[i] = d1d2$weight_exact[i]
  
}

exp = data.frame(config = config, y = yt, n = nt)

#replication may occur with  d1d$design_exact need to check
if(sum(dim(d1d2$design_exact)!=dim(unique(d1d2$design_exact)))>0){
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
theta_hat_pos2 = apply(sim.post(nruns = 1000,ndraws = 500,n_acc =exp$n ,y_acc = exp$y,X,theta_hat_pos, last_N = 500),MARGIN = 1,FUN = mean)


arm_opt_mult_rpm[R]=which.max(pnorm(X%*%theta_hat_pos2))
