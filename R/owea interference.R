# interference model # 
rm(list=ls())

# non-circular design #
# information matrix #
infor_point <- function(point, # design point, row vector
                        t, # treatments
                        p,
                        sigma# block size
){
  X <- matrix(0,ncol=3*t,nrow=p)
  Td <- X[,1:t] <- t(sapply(point,function(x){as.numeric(1:t %in% x)}))
  X[2:p,t+1:t] <- Td[-p,]
  X[1:(p-1),2*t+1:t] <- Td[-1,]
  sig.inv <- MASS::ginv(sigma)
  ii <- matrix(1,ncol=1,nrow=p)
  mid <- sig.inv-sig.inv%*%ii%*%t(ii)%*%sig.inv/sum(sig.inv)
  out <- t(X)%*%mid%*%X
}

infor_design <- function(design,weight,t,p,sigma){
  out <- matrix(0,nrow=3*t,ncol=3*t)
  for(i in 1:NROW(design)){
    out <- out + weight[i]*infor_point(design[i,],t,p,sigma)
  }
  return(out)
}

# requires partial_g function to proceed #

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
dw <- function(design,weight,I,covar,t,p,sigma,pp){
  d1w <- matrix(0,nrow=NROW(design)-1)
  d2w <- matrix(0,ncol=NROW(design)-1,nrow=NROW(design)-1)
  I.inv <- MASS::ginv(I)
  Im <- infor_point(design[nrow(design),],t,p,sigma)
  d1w_part <- list()
  # D-optimal when pp=0 #
  if(pp==0){
    covar.inv <- MASS::ginv(covar)
    for (i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],t,p,sigma)-Im
      d1w_part[[i]] <- -(g_part%*%I.inv%*%(mid_i)%*%I.inv%*%t(g_part))
      d1w[i,1] <- sum(diag(covar.inv%*%d1w_part[[i]]))}
    for(i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],t,p,sigma)-Im
      for(j in 1:(NROW(design)-1)){
        mid_j <- infor_point(design[j,],t,p,sigma)-Im
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
      mid_i <- infor_point(design[i,],t,p,sigma)-Im
      d1w_part[[i]] <- -(g_part%*%I.inv%*%(mid_i)%*%I.inv%*%t(g_part))
      d1w[i,1] <- sum(diag(d1w_part[[i]]))}
    # end of first order #
    for(i in 1:(NROW(design)-1)){
      mid_i <- infor_point(design[i,],t,p,sigma)-Im
      for(j in 1:(NROW(design)-1)){
        mid_j <- infor_point(design[j,],t,p,sigma)-Im
        mid_ij <- mid_i%*%I.inv%*%mid_j
        mid_ji <- mid_j%*%I.inv%*%mid_i
        d2w[i,j] <- sum(diag(g_part%*%(I.inv%*%(mid_ji+mid_ij)%*%I.inv)%*%t(g_part)))
      }
    }
  }
  return(list(d1w=d1w,d2w=d2w))
}

# directional derivative #
dp <- function(pp,point_new,design,I.inv,covar,t,p,sigma,weight){
  A <- g_part%*%(I.inv%*%(infor_point(point_new,t,p,sigma)-
                            infor_design(design,weight,t,p,sigma))%*%I.inv)%*%t(g_part)
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
weight1 <- function(design,weight,t,p,sigma,pp){
  weight_long <- matrix(weight,ncol=1)
  weight_short <- matrix(weight_long[-length(weight_long)],ncol=1)
  #initialization parameters #
  diff <- 1
  alpha <- 1
  # initialzation ends #
  while(diff>0.000000000000001){
    I <- infor_design(design,weight_long,t,p,sigma) 
    covar <- g_part%*%MASS::ginv(I)%*%t(g_part) 
    tmp <- dw(design,weight_long,I,covar,t,p,sigma,pp)
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
weight2 <- function(design,weight,t,p,sigma,pp){
  weight_opt <- matrix(weight,ncol=1)
  tmp1 <- weight1(design,weight_opt,t,p,sigma,pp)
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
    tmp1 <- weight1(design,weight_opt,t,p,sigma,pp)
  }
  if(NROW(design)==1){weight_opt = 1}
  return(list(design=as.matrix(design),weight=as.matrix(weight_opt)))
}


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

# test example #
# example 1 D-optimal #
pp <- 0
t <- 4 
p <- 4 
n <- 10
g_part <- matrix(0,nrow = t-1, ncol = 3*t)
g_part[,1:t] <- cbind(1,diag(-1,t-1))
sigma <- diag(1,p)
space <- gtools::permutations(t,p,repeats.allowed = T)
d1d <- opt_design(t,p,sigma,n,pp,space)

