# interference model # 
#single-stage


rm(list=ls())

# non-circular design #
# information matrix #

#modified for probit model information matrix
#here "point" is a row of the design matrix eg.(0,1,0,1,0,0,...,0), 0's and 1's, each position corresponds to a 
#             level of a factor. 
infor_point <- function(point, #a row (or rows) from the design, 1 by number of parameters
                        #t, # treatments
                        #p, ## (YH number of periods?)
                        #sigma# block size #(YH:covariance between periods?)
                        theta_stg1 #tentatively estimated theta(from stage 1)
){
  point = as.matrix(point,nrow = length(theta_stg1))
  s = as.numeric(t(point)%*%theta_stg1_0)
  infor = (dnorm(s))^2/pnorm(s)/(1-pnorm(s))*point%*%t(point)
  return(infor)
}

infor_design <- function(design_stg1=NULL, #design in stage 1: n1 by number of parameters matrix
                         design, # new design generated and updated in stage 2
                         weight, # weights of support points in stage2, updating
                         theta_stg1) # estimate of theta from stage 1
                         #t,p,
                         #sigma)
  { 
  
  #calculate information from stage1: 
  #here we assume stage 1 each point is used exactly once to gain rough estimate
  t = dim(design)[2]
  if (!is.null(design_stg1)){
    out <- matrix(0,nrow=t,ncol=t)
    for (i in 1:dim(design_stg1)[1]){
      out = out+ infor_point(design_stg1[i,], theta_stg1)
    }
    out0 = out
    a1 = dim(design_stg1)[1]/(dim(design_stg1)[1]+n)
    a2 = n/(dim(design_stg1)[1]+n)
  }else{
    out0 = matrix(0,nrow=t,ncol=t)
    a1 = 1
    a2 = 1
  }
  
  
  out <- matrix(0,nrow=t,ncol=t)
  for(i in 1:NROW(design)){
    out <- out + weight[i]*infor_point(design[i,],theta_stg1)#weight needs to be adjusted n0>0
  }
  out1 = out
  out = a1*out1+a2*out0
  return(out)
}# can be used to compute the info_mat for both old design and new design, (and add up)

# requires partial_g function to proceed # now assume target is to estimate theta_stg1

# phi_p function #
phi <- function(covar,pp){
  #infor.treat <- infor_design(design,t,p,drop,weight)[(t+2):(2*t+1),(t+2):(2*t+1)] 
  #covar <- MASS::ginv(infor_design(design,t,p))
  if (pp == 0){
    out <- log(det(covar))
  }else{
    out <- sum(diag(covar))^pp
  }
  return(out)
}

# derivative of phi_p #
# derivatives of both 1st and 2nd order #
dw <- function(design,weight,I,covar,theta_stg1,pp){
  d1w <- matrix(0,nrow=NROW(design)-1)
  d2w <- matrix(0,ncol=NROW(design)-1,nrow=NROW(design)-1)
  I.inv <- MASS::ginv(I)
  Im <- infor_point(design[nrow(design),],theta_stg1) #????
  
  d1w_part <- list()
  # D-optimal when pp=0 #
  if(pp==0){
    covar.inv <- MASS::ginv(covar)
    for (i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],theta_stg1)-Im
      d1w_part[[i]] <- -(g_part%*%I.inv%*%(mid_i)%*%I.inv%*%t(g_part)) # A(g, xi)
      d1w[i,1] <- sum(diag(covar.inv%*%d1w_part[[i]]))}
    for(i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],theta_stg1)-Im
      for(j in 1:(NROW(design)-1)){
        mid_j <- infor_point(design[j,],theta_stg1)-Im
        mid_ij <- mid_i%*%I.inv%*%mid_j
        mid_ji <- mid_j%*%I.inv%*%mid_i
        d2w[i,j] <- sum(diag(covar.inv%*%g_part%*%(I.inv%*%(mid_ji+mid_ij)%*%I.inv)%*%t(g_part)))-
          sum(diag(covar.inv%*%d1w_part[[j]]%*%covar.inv%*%d1w_part[[i]]))
      }
    }
  }
  # A-optimal #
  if(pp==1){
    # first order d1w_part #
    for (i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],theta_stg1)-Im
      d1w_part[[i]] <- -(g_part%*%I.inv%*%(mid_i)%*%I.inv%*%t(g_part))
      d1w[i,1] <- sum(diag(d1w_part[[i]]))}
    # end of first order #
    for(i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],theta_stg1)-Im
      for(j in 1:(NROW(design)-1)){
        mid_j <- infor_point(design[j,],theta_stg1)-Im
        mid_ij <- mid_i%*%I.inv%*%mid_j
        mid_ji <- mid_j%*%I.inv%*%mid_i
        d2w[i,j] <- sum(diag(g_part%*%(I.inv%*%(mid_ji+mid_ij)%*%I.inv)%*%t(g_part)))
      }
    }
  }
  return(list(d1w=d1w,d2w=d2w))
}

# directional derivative #
dp <- function(pp,point_new,design,design_stg1,I.inv,covar,theta_stg1,weight){
  A <- g_part%*%(I.inv%*%(infor_point(point_new,theta_stg1)-
                            infor_design(design_stg1 = NULL,design,weight,theta_stg1))%*%I.inv)%*%t(g_part)
  if(pp==0){
    out <- sum(diag(MASS::ginv(covar)%*%A))
  }
  if(pp == 1){
    nu <- dim(g_part)[1]
    # calculating matrix power using eigen values #
    out <- (1/nu)*sum(diag(A))
  }
  return(out)
}

# weight 1 returns positive critical points of weight #
weight1 <- function(design,design_stg1,weight,theta_stg1,pp){
  weight_long <- matrix(weight,ncol=1)
  weight_short <- matrix(weight_long[-length(weight_long)],ncol=1)
  #initialization parameters #
  diff <- 1
  alpha <- 1
  # initialzation ends #
  while(diff>0.000000000000001){
    I <- infor_design(design_stg1 ,design,weight_long,theta_stg1) 
    covar <- g_part%*%MASS::ginv(I)%*%t(g_part) 
    tmp <- dw(design,weight_long,I,covar,theta_stg1,pp)
    #print(weight_short)
    # update the weight #
    weight_new <- weight_short - alpha*MASS::ginv(tmp$d2w)%*%tmp$d1w
    # new weight #
    weight_long <- matrix(c(weight_new,1-sum(weight_new)),ncol=1) 
    # check non-positive weights #
    if(min(weight_long)<0){
      if(alpha > 0.00000001){
        alpha <- alpha/2
      }else{
        break
      }
    }else{
      # all weights are positive #
      # tmp <- dw(design,weight_long,I,covar,t,p,drop,pp)
      # diff <- sum((tmp$d1w)^2)
      diff <- sum((weight_new-weight_short)^2)
      weight_short <- weight_new
      weight_long <- matrix(c(weight_new,1-sum(weight_new)),ncol=1)
    }
  }
  # returns the design with all postive weights
  return(list(design = as.matrix(design),weight=weight_long))
}

# critical points on boundary #
weight2 <- function(design,design_stg1,weight,theta_stg1,pp){
  weight_opt <- matrix(weight,ncol=1)
  tmp1 <- weight1(design,design_stg1,weight_opt,theta_stg1,pp)
  weight_opt <- tmp1$weight
  design <- tmp1$design
  while(NROW(design)>1 & min(weight_opt)<0.000001){
    #print(tmp1)
    if(min(weight_opt)<0.000001){
      #print('remove points at boundary')
      index.rm <- which.min(weight_opt)
      weight_opt <- weight_opt[-index.rm]
      design <- design[-index.rm,]
    }
    tmp1 <- weight1(design,design_stg1,weight_opt,theta_stg1,pp)
  }
  if(NROW(design)==1){weight_opt = 1}
  return(list(design=as.matrix(design),weight=as.matrix(weight_opt)))
}


opt_design <- function(design_stg1,theta_stg1,n,pp,space){
  p=length(theta_stg1)
  n<<-n
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
      tmp <- weight2(design,design_stg1,weight_opt,theta_stg1,pp)
      # newton optimization #
      weight_opt <- tmp$weight
      design <- tmp$design
      I <- infor_design(design_stg1,design,weight_opt,theta_stg1)
      I.inv <- MASS::ginv(I)
      covar <- g_part%*%I.inv%*%t(g_part)
      # find x* = argmax dp(x*,design) #
      dp_value <- rep(0,NROW(point_space))
      for (i in 1:NROW(point_space)){
        dp_value[i] <- dp(pp,point_space[i,],design,design_stg1,I.inv,covar,theta_stg1,weight_opt)
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
    
    # comprehensive rounding #  ##?????? ignore for now
#    for(i in 1:NROW(space)){
#      phi.old <- phi(g_part%*%MASS::ginv(infor_design(design_stg1,design_exact,weight_exact,theta_stg1))%*%t(g_part),pp)
#      for(j in 1:NROW(design_exact)){
#        design_tmp <- rbind(design_exact,space[i,])
#        weight_tmp <- weight_exact
#        weight_tmp[j] <- weight_tmp[j] - 1
#        weight_tmp <- c(weight_tmp,1)
#        tmp.zero <- weight_tmp==0
#        weight_tmp <- weight_tmp[!tmp.zero]
#        design_tmp <- design_tmp[!tmp.zero,]
 #       phi.new <- phi(g_part%*%MASS::ginv(infor_design(design_stg1,design_tmp,weight_tmp,theta_stg1))%*%t(g_part),pp)
#        if(phi.new - phi.old < 0){
#          tmp.zero <- weight_tmp==0
#          weight_exact <- weight_tmp[!tmp.zero]
#          design_exact <- design_tmp[!tmp.zero,]
#        }               
#      }
#    }
  })  
  out <- list(design_approx=design,weight_approx=weight_opt,design_exact=design_exact,weight_exact=weight_exact,dp=d_p,
              time_approx = time1, time_exact = time2)
  print(out)
  return(out)
}
