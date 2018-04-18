####### use wilcox test  #####


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






####### use top 5 ranking #####


out = order(pnorm(X%*%theta_hat),decreasing = T)[1:5]


ind_no = c()
dat = X[out,]

i=2
while (i<= length(theta)){
  if (i==2) {
    case = 2
  }else if(i %in% 3:4){
    case = 3:4
  }else if(i %in% 5:7){
    case = 5:7
  }else {
    case =8:11
  }
  
  #if(NROW(unique(dat[,case]))<=2) ind_no = c(ind_no,case)
  if(max(count(as.data.frame(dat),vars = apply(as.matrix(case),MARGIN = 1,FUN = funpaste))$freq)>=4) {
    ind_no = c(ind_no,case)
  }
  i = max(case)+1 
}
#record conclusions 
theta_hat_no = theta_hat[ind_no]

ind_no = setdiff(1:length(theta),ind_no)#param that  need to estimate further
theta_hat_con = theta_hat[ind_no]

#browser()
g_part = matrix(0,nrow = length(ind_no), ncol = length(theta))
for ( i in 1:NROW(g_part)){
  g_part[i,ind_no[i]] =1
  
}
g_part[,1] = rep(1,length(ind_no))


####  use contrast of the mu's ####

out = order(pnorm(X%*%theta_hat),decreasing = T)[1:5]

dat = X[out,]
for (i in 1:NROW(dat)){
  
}



