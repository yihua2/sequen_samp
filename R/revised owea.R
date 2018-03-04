opt_design <- function(t,p,sigma,n,pp,space){
  point_space <- space
  index.initial <- sample(1:NROW(point_space),p) # initial selection
  design <- point_space[index.initial,] # initial design
  print(c('initial selection:'))
  print(index.initial)
  print(c('initial design is:'))
  print(design)
  # initialize parameters #
  d_p <- 1 
  weight_opt <- matrix(1/(p),ncol = 1, nrow = p)
  iter <- 0
  # initialization ends #
  
  
  
  time1 <- system.time({
    while(d_p>0.00001 & iter <= 200){
      iter <- iter + 1
      # x* obtained #
      tmp <- weight2(design,weight_opt,t,p,sigma,pp)
      # newton optimization #
      weight_opt <- tmp$weight
      design <- tmp$design
      I <- infor_design(design,weight_opt,t,p,sigma)
      I.inv <- MASS::ginv(I)
      covar <- g_part%*%I.inv%*%t(g_part)
      # find x* = argmax dp(x*,design) #
      dp_value <- rep(0,NROW(point_space))
      for (i in 1:NROW(point_space)){
        dp_value[i] <- dp(pp,point_space[i,],design,I.inv,covar,t,p,sigma,weight_opt)
        point_new <- point_space[which.max(dp_value),]
      }
      d_p <- max(dp_value)
      #print(d_p)
      if (d_p < 0.00001){
        print('optimal')
        break
      }else{
        design <- rbind(design,point_new)
        weight_opt <- rbind(weight_opt,0)
        next
      }
    }
  }) 
  
  # round to exact design
  time2=system.time({
    weight_opt1 <- round(weight_opt*n)
    weight_opt2 <- weight_opt*n-weight_opt1
    weight_opt3 <- matrix(0,nrow=NROW(weight_opt1),ncol=NCOL(weight_opt1))
    diffn <- n-sum(weight_opt1)
    if(abs(diffn) > 0){
      for (i in 1:abs(diffn)){
        if(diffn > 0){
          # need extra points
          tmp.ind <- which.max(weight_opt2)
          weight_opt3[tmp.ind] <- 1
          weight_opt2[tmp.ind] <- 0
        }
        if(diffn < 0){
          # need remove points
          tmp.ind <- which.min(weight_opt2)
          weight_opt3[tmp.ind] <- -1
          weight_opt2[tmp.ind] <- 0
        }
      }
    }
    
    
    # summarying exact design and weights #
    weight_exact <- weight_opt1+weight_opt3
    tmp.zero <- weight_exact==0
    weight_exact <- weight_exact[!tmp.zero]
    design_exact <- design[!tmp.zero,]
    
    # comprehensive rounding #
    for(i in 1:NROW(space)){
      phi.old <- phi(g_part%*%MASS::ginv(infor_design(design_exact,weight_exact,t,p,sigma))%*%t(g_part),pp)
      for(j in 1:NROW(design_exact)){
        design_tmp <- rbind(design_exact,space[i,])
        weight_tmp <- weight_exact
        weight_tmp[j] <- weight_tmp[j] - 1
        weight_tmp <- c(weight_tmp,1)
        tmp.zero <- weight_tmp==0
        weight_tmp <- weight_tmp[!tmp.zero]
        design_tmp <- design_tmp[!tmp.zero,]
        phi.new <- phi(g_part%*%MASS::ginv(infor_design(design_tmp,weight_tmp,t,p,sigma))%*%t(g_part),pp)
        if(phi.new - phi.old < 0){
          tmp.zero <- weight_tmp==0
          weight_exact <- weight_tmp[!tmp.zero]
          design_exact <- design_tmp[!tmp.zero,]
        }               
      }
    }
  })  
  out <- list(design_approx=design,weight_approx=weight_opt,design_exact=design_exact,weight_exact=weight_exact,dp=d_p,
              time_approx = time1, time_exact = time2)
  print(out)
  return(out)
}
