#functions
estimate_multi <- function(W){
  
  n = dim(W)[1]
  rp = dim(W)[2]
  W=W
  index <- t(combn(1:rp,3))
  # matrix of Wi*Wj
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  
  #matrix of Wi*wj*wk
  #Wijk<-as.matrix(apply(W,1,prod))
  #Wijk <- matrix(nrow=n,ncol=1)
  
  
  
  #matrix of Wi2*wj2*wk2
  #Wijk2<-as.matrix(apply(W^2,1,prod))
  
  
  
  # matrix of W1*W2
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  
  # matrix of W^2
  W2 <- W^2
  
  # matrix of W^3
  W3 <- W^3
  
  # matrix of W^4
  W4 <- W^4
  
  Wijk2 <- matrix(nrow=n,ncol=choose(rp,3))
  for(i in c(1:choose(rp,3))){
    Wijk2[,i] = as.matrix(apply(W2[,index[i,]],1,prod))
  }
  
  Wijk <- matrix(nrow=n,ncol=choose(rp,3))
  for(i in c(1:choose(rp,3))){
    Wijk[,i] = as.matrix(apply(W[,index[i,]],1,prod))
  }
  
  # matrix of Wi^2*Wj
  W2W <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W2W[,k] <- W2[,i]*W[,j]
      k=k+1
    }
  }
  
  # matrix of Wi^3*Wj
  W3W <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W3W[,k] <- W3[,i]*W[,j]
      k=k+1
    }
  }
  
  # matrix of Wi^3*Wj^2
  W3W2 <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W3W2[,k] <- W3[,i]*W2[,j]
      k=k+1
    }
  }
  
  # matrix of Wi^2*Wj^2
  W2W2 <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W2W2[,k] <- W2[,i]*W2[,j]
      k=k+1
    }
  }
  
  # matrix of Wi^4*Wj
  W4W <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W4W[,k] <- W4[,i]*W2[,j]
      k=k+1
    }
  }
  
  
  
  #matrix of  w1^2*W2*W3
  
  Wi2jk <- matrix(nrow=n,ncol=choose(rp,3))
  for(i in c(1:choose(rp,3))){
    Wi2jk[,i] = as.matrix(apply(W[,index[i,2:3]],1,prod))*W2[,index[i,1]]
  }
  
  
  #matrix of  w1^3*W2^2*W3
  
  
  Wi3j2k <- matrix(nrow=n,ncol=choose(rp,3))
  for(i in c(1:choose(rp,3))){
    Wi3j2k[,i] =W[,index[i,1]]*W2[,index[i,2]]*W3[,index[i,3]]
  }
  
  
  #matrix of  w1^2*W2^2*W3

  Wi2j2k <- matrix(nrow=n,ncol=choose(rp,3))
  for(i in c(1:choose(rp,3))){
    Wi2j2k[,i] = as.matrix(apply(W2[,index[i,2:3]],1,prod))*W[,index[i,1]]
  }
  
  Wi3jk <- matrix(nrow=n,ncol=choose(rp,3))
  for(i in c(1:choose(rp,3))){
    Wi3jk[,i] = as.matrix(apply(W[,index[i,2:3]],1,prod))*W3[,index[i,1]]
  }
  
  
  # estimate of a2
  a2 <- (mean(W2)*mean(W)-mean(W2W))/(mean(Wij)*mean(W)-mean(Wijk))
  
  # estimate of sig_mul
  sig_mul <- sqrt(2*log(sqrt(a2)))
  
  # estimate of sig_add_2^2
  sig_add_2 <- mean(W2)-a2*mean(Wij)
  
  # estimate of sig_add
  sig_add <- sqrt(mean(W2)-a2*mean(Wij))
  
  mu1 <- mean(W)/sqrt(a2)
  
  mu2 <- (mean(Wij))/(a2)
  
  result <- list(a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add, EWW=mean(Wij), EW=mean(W), W2=W2, W2W=W2W,n=n,Wij=Wij,W=W,Wijk=Wijk,Wijk2=Wijk2,W2=W2,W3=W3,W4=W4, W3W=W3W,
                 W3W2=W3W2,W2W2=W2W2,W4W=W4W,Wi2jk=Wi2jk,Wi3jk=Wi3jk,Wi2j2k=Wi2j2k,Wi3j2k=Wi3j2k,mu1=mu1,mu2=mu2)
  return(result)
}
estimate_triple <- function(W){
  
  n = dim(W)[1]
  rp = dim(W)[2]
  W=W
  
  # matrix of Wi*Wj
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  
  #matrix of Wi*wj*wk
  Wijk<-as.matrix(apply(W,1,prod))
  #Wijk <- matrix(nrow=n,ncol=1)
  
  #matrix of Wi2*wj2*wk2
  Wijk2<-as.matrix(apply(W^2,1,prod))
  
  
  # matrix of W1*W2
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  
  # matrix of W^2
  W2 <- W^2
  
  # matrix of W^3
  W3 <- W^3
  
  # matrix of W^4
  W4 <- W^4
  
  # matrix of Wi^2*Wj
  W2W <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W2W[,k] <- W2[,i]*W[,j]
      k=k+1
    }
  }
  
  # estimate of a2
  a2 <- (mean(W2)*mean(W)-mean(W2W))/(mean(Wij)*mean(W)-mean(apply(W,1,prod)))
  
  # estimate of sig_mul
  sig_mul <- sqrt(2*log(sqrt(a2)))
  
  # estimate of sig_add_2^2
  sig_add_2 <- mean(W2)-a2*mean(Wij)
  
  # estimate of sig_add
  sig_add <- sqrt(mean(W2)-a2*mean(Wij))
  
  result <- list(a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add,EWW=mean(Wij), EW=mean(W), W2=W2, W2W=W2W,n=n,Wij=Wij,W=W,Wijk=Wijk)
  return(result)
}
hypothesis_triple <- function(W,Type){
  
  result_est <- estimate_triple(W)
  
  
  a2_est <- result_est$a2
  
  sig_add_2_est <- result_est$sig_add_2
  
  EW<-result_est$EW
  
  EWW<-result_est$EWW
  
  W2<-result_est$W2
  
  W2W<-result_est$W2W
  
  Wij<-result_est$Wij
  
  Wijk<-result_est$Wijk
  
  n<-result_est$n
  
  Wijk=result_est$Wijk
  Wijk2=result_est$Wijk2
  
  W2=result_est$W2
  W3=result_est$W3
  W4=result_est$W4
  W3W=result_est$W3W
  W3W2=result_est$W3W2
  W2W2=result_est$W2W2
  W4W=result_est$W4W
  Wi2jk=result_est$Wi2jk
  Wi3jk=result_est$Wi3jk
  Wi2j2k=result_est$Wi2j2k
  Wi3j2k=result_est$Wi3j2k
  
  if(a2_est<1){a2_est <- 1}
  
  if(sig_add_2_est<0) {sig_add_2_est<-0}
  #W<-result_est$W
  
  innerfunction <- function(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk){
    
    #####################
    #change the variance#
    #####################
    
    g1_i <- sig_add^2*apply(W, 1, mean) + a2*Wijk - apply(W2W,1,mean)
    
    g2_i <- sig_add^2 - apply(W2,1,mean) + a2*apply(Wij,1,mean)
    
    v1 <- cov(data.frame(g1=g1_i,g2=g2_i))
    
    v<-matrix(nrow=2, ncol=2)
  
    
    ########################################
    # Note there are four ways to set gamma# 
    ########################################
    
    gamma <- matrix(c(mean(Wijk),EW,EWW,1),nrow = 2,byrow = TRUE)
    #gamma <- matrix(c(2*sqrt(a2)*mean(Wijk),EW,2*sqrt(a2)*EWW,1),nrow = 2,byrow = TRUE)
    #gamma <- matrix(c(2*sqrt(a2)*mean(Wijk),2*sig_add*EW,2*sqrt(a2)*EWW,1),nrow = 2,byrow = TRUE)
    #gamma <- matrix(c(2*sqrt(a2)*mean(Wijk)+sig_add^2*EW,EW*sqrt(a2),2*sqrt(a2)*EWW,1),nrow = 2,byrow = TRUE)
    #cov_est <- ginv(gamma)%*%v2%*%t(ginv(gamma))/n
    cov_est <- ginv(gamma)%*%v1%*%t(ginv(gamma))/n
    #cov_est <- solve(gamma)%*%v2%*%t(solve(gamma))/n
    return(cov_est)
  }
  
  # Type 1 HP
  if(Type==1){
    
    a2 <- 1
    
    sig_add <- 0
    
    
    cov_est <-innerfunction(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk)
    d2 <- t(c(a2_est-1, sig_add_2_est))%*%ginv(cov_est)%*%c(a2_est-1, sig_add_2_est)
    p_value <- pchisq(d2,2,lower.tail = FALSE)/4
    
  }
  
  if(Type==2){
    
    
    a2 <- 1
    
    sig_add <- sqrt(sig_add_2_est)
    
    
    cov_est <-innerfunction(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk)
    p_value <- 1 - pnorm(a2_est,1,sqrt(cov_est[1,1]))
    #p_value <- 1 - phalfnorm((sqrt(a2_est)-1)/sqrt(cov_est[1,1]))
  }
  if(Type==3){
    a2 <- a2_est
    sig_add <-0
    cov_est <-innerfunction(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk)
    p_value <- 1 - pnorm(sig_add_2_est/sqrt(cov_est[2,2]))
  }
  
  return(p_value)
}
hypothesis_multi <- function(W,Type){
  
  result_est <- estimate_multi(W)
  
  
  a2_est <- result_est$a2
  
  sig_add_2_est <- result_est$sig_add_2
  
  EW<-result_est$EW
  
  EWW<-result_est$EWW
  
  W2<-result_est$W2
  
  W2W<-result_est$W2W
  
  Wij<-result_est$Wij
  
  Wijk<-result_est$Wijk
  
  n<-result_est$n
  
  Wijk=result_est$Wijk
  Wijk2=result_est$Wijk2
  
  W2=result_est$W2
  W3=result_est$W3
  W4=result_est$W4
  W3W=result_est$W3W
  W3W2=result_est$W3W2
  W2W2=result_est$W2W2
  W4W=result_est$W4W
  Wi2jk=result_est$Wi2jk
  Wi3jk=result_est$Wi3jk
  Wi2j2k=result_est$Wi2j2k
  Wi3j2k=result_est$Wi3j2k
  
  if(a2_est<1){a2_est <- 1}
  
  if(sig_add_2_est<0) {sig_add_2_est<-0}
  #W<-result_est$W
  
  innerfunction <- function(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk){
    
    #####################
    #change the variance#
    #####################
    
    g1_i <- sig_add^2*apply(W, 1, mean) + a2*apply(Wijk,1,mean) - apply(W2W,1,mean)
    
    g2_i <- sig_add^2 - apply(W2,1,mean) + a2*apply(Wij,1,mean)
    
    v1 <- cov(data.frame(g1=g1_i,g2=g2_i))
    
    v<-matrix(nrow=2, ncol=2)
    
    
    ########################################
    # Note there are four ways to set gamma# 
    ########################################
    
    gamma <- matrix(c(mean(Wijk),EW,EWW,1),nrow = 2,byrow = TRUE)
    #gamma <- matrix(c(2*sqrt(a2)*mean(Wijk),EW,2*sqrt(a2)*EWW,1),nrow = 2,byrow = TRUE)
    #gamma <- matrix(c(2*sqrt(a2)*mean(Wijk),2*sig_add*EW,2*sqrt(a2)*EWW,1),nrow = 2,byrow = TRUE)
    #gamma <- matrix(c(2*sqrt(a2)*mean(Wijk)+sig_add^2*EW,EW*sqrt(a2),2*sqrt(a2)*EWW,1),nrow = 2,byrow = TRUE)
    #cov_est <- ginv(gamma)%*%v2%*%t(ginv(gamma))/n
    cov_est <- ginv(gamma)%*%v1%*%t(ginv(gamma))/n
    #cov_est <- solve(gamma)%*%v2%*%t(solve(gamma))/n
    return(cov_est)
  }
  
  # Type 1 HP
  if(Type==1){
    
    a2 <- 1
    
    sig_add <- 0
    
    
    cov_est <-innerfunction(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk)
    d2 <- t(c(a2_est-1, sig_add_2_est))%*%ginv(cov_est)%*%c(a2_est-1, sig_add_2_est)
    p_value <- pchisq(d2,2,lower.tail = FALSE)/4
    
  }
  
  if(Type==2){
    
    
    a2 <- 1
    
    sig_add <- sqrt(sig_add_2_est)
    
    
    cov_est <-innerfunction(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk)
    p_value <- 1 - pnorm(a2_est,1,sqrt(cov_est[1,1]))
    #p_value <- 1 - phalfnorm((sqrt(a2_est)-1)/sqrt(cov_est[1,1]))
  }
  if(Type==3){
    a2 <- a2_est
    sig_add <-0
    cov_est <-innerfunction(a2,sig_add,EW,EWW,W2,W2W,Wij,n,W,Wijk)
    p_value <- 1 - pnorm(sig_add_2_est/sqrt(cov_est[2,2]))
  }
  
  return(p_value)
}
estimate_regression <- function(W,Z,Y){
  n <- dim(W)[1]
  rp = dim(W)[2]
  estimate <- estimate_double(W)
  W_m <- apply(W,1,mean)
  alpha1 <- (sqrt(estimate$a2))*(estimate$mu2-estimate$mu1^2)/(var(W_m))
  alpha0 <- estimate$mu1-mean(W_m)*alpha1
  x_hat <- alpha1*W_m+alpha0
  if(is.na(Z)) data <- data.frame(x_hat) else data <- data.frame(x_hat,Z)
  fit1 <- lm(Y~.,data = data)
  beta1 <- fit1$coefficients[2:length(fit1$coefficients)]
  beta0 <- fit1$coefficients[1]
  # matrix of Wi*Wj
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  
  # matrix of W^2
  W2 <- W^2
  
  # matrix of W^3
  W3 <- W^3
  
  # matrix of Wi^2*Wj
  W2W <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W2W[,k] <- W2[,i]*W[,j]
      k=k+1
    }
  }
  phi1 <- Y-beta0-as.matrix(data)%*%beta1
  phi2 <- as.matrix(data)*c(Y-beta0-as.matrix(data)%*%beta1)
  phi3 <- estimate$a2^2*(W_m*estimate$sig_add_2-apply(W2W,1,mean))-3*W_m*estimate$sig_add_2+apply(W3,1,mean)
  phi4 <- estimate$sig_add_2-apply(W2,1,mean)+estimate$a2*apply(Wij,1,mean)
  phi <- cbind(phi1,phi2,phi3,phi4)
  Bn <- t(phi)%*%phi/n
  #Bn11 <- mean((Y-beta0-beta1*x_hat)^2)
  #Bn12 <- mean((Y-beta0-beta1*x_hat)^2*x_hat)
  #Bn22 <- mean((Y-beta0-beta1*x_hat)^2*x_hat^2)
  #Bn <- matrix(c(Bn11,Bn12,Bn12,Bn22),byrow = TRUE,nrow = 2)
  da2_alpha_1 <- ((mean(W2)-estimate$sig_add_2)*(-1.5)*(estimate$a2^(-2.5))-mean(W)^2*(-0.5)*(estimate$a2^(-1.5)))/(var(W_m))
  da2_x_hat <- (W_m-mean(W_m))*da2_alpha_1-0.5*mean(W_m)*(estimate$a2^(-1.5))
  dsig_alpha_1 <- -estimate$a2^(-1.5)/(var(W_m))
  dsig_x_hat <- (W_m-mean(W_m))*dsig_alpha_1
  # An12 <- -mean(x_hat)
  # An22 <- -mean(x_hat^2)
  An211 <- -mean(beta1[1]*da2_x_hat)
  An212 <- -mean(beta1[1]*dsig_x_hat)
  An221 <- mean((Y-beta0-as.matrix(data)%*%beta1)*da2_x_hat-beta1[1]*da2_x_hat*data$x_hat)
  An222 <- mean((Y-beta0-as.matrix(data)%*%beta1)*dsig_x_hat-beta1[1]*dsig_x_hat*data$x_hat)
  if(is.na(Z[1])){
    An2 <- matrix(c(An211,An212,An221,An222),byrow = TRUE,nrow = 2)
    
  }else{
    An231 <- mean(-beta1[1]*da2_x_hat*Z)
    An232 <- mean(-beta1[1]*dsig_x_hat*Z)
    An2 <- matrix(c(An211,An212,An221,An222,An231,An232),byrow = TRUE,nrow = 3)
  }
  An1 <- -t(as.matrix(cbind(1,data)))%*%as.matrix(cbind(1,data))/n
  An411 <- 2*estimate$a2*(mean(W)*estimate$sig_add_2-mean(W2W))
  An412 <- estimate$a2^2*mean(W)-3*mean(W)
  An421 <- mean(Wij)
  An422 <- 1
  # An <- matrix(c(-1,An12,An13,An14,An12,An22,An23,An24,rep(0,2),An33,An34,rep(0,2),An43,An44),byrow = TRUE,nrow = 4)
  An4 <- matrix(c(An411,An412,An421,An422),byrow = TRUE,nrow = 2)
  An3 <- matrix(0,dim(An4)[1],dim(An1)[2])
  An <- rbind(cbind(An1,An2),cbind(An3,An4))
  cov_est <- solve(An)%*%Bn%*%t(solve(An))/n
  return(list(beta1=beta1,beta0=beta0,cov_est=cov_est,alpha0=alpha0,alpha1=alpha1))
} 
estimate_regression_multi <- function(W,Z,Y){
  n <- dim(W)[1]
  rp = dim(W)[2]
  estimate <- estimate_multi(W)
  W_m <- apply(W,1,mean)
  alpha1 <- (sqrt(estimate$a2))*(estimate$mu2-estimate$mu1^2)/(var(W_m))
  alpha0 <- estimate$mu1-mean(W_m)*alpha1
  x_hat <- alpha1*W_m+alpha0
  if(is.na(Z)) data <- data.frame(x_hat) else data <- data.frame(x_hat,Z)
  fit1 <- lm(Y~.,data = data)
  beta1 <- fit1$coefficients[2:length(fit1$coefficients)]
  beta0 <- fit1$coefficients[1]
  # matrix of Wi*Wj
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  
  # matrix of W^2
  W2 <- W^2
  
  # matrix of W^3
  W3 <- W^3
  
  # matrix of Wi^2*Wj
  W2W <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  for(i in c(1:(rp))){
    for(j in c(1:rp)[-i]){
      W2W[,k] <- W2[,i]*W[,j]
      k=k+1
    }
  }
  phi1 <- Y-beta0-as.matrix(data)%*%beta1
  phi2 <- as.matrix(data)*c(Y-beta0-as.matrix(data)%*%beta1)
  phi3 <- estimate$a2^2*(W_m*estimate$sig_add_2-apply(W2W,1,mean))-3*W_m*estimate$sig_add_2+apply(W3,1,mean)
  phi4 <- estimate$sig_add_2-apply(W2,1,mean)+estimate$a2*apply(Wij,1,mean)
  phi <- cbind(phi1,phi2,phi3,phi4)
  Bn <- t(phi)%*%phi/n
  #Bn11 <- mean((Y-beta0-beta1*x_hat)^2)
  #Bn12 <- mean((Y-beta0-beta1*x_hat)^2*x_hat)
  #Bn22 <- mean((Y-beta0-beta1*x_hat)^2*x_hat^2)
  #Bn <- matrix(c(Bn11,Bn12,Bn12,Bn22),byrow = TRUE,nrow = 2)
  da2_alpha_1 <- ((mean(W2)-estimate$sig_add_2)*(-1.5)*(estimate$a2^(-2.5))-mean(W)^2*(-0.5)*(estimate$a2^(-1.5)))/(var(W_m))
  da2_x_hat <- (W_m-mean(W_m))*da2_alpha_1-0.5*mean(W_m)*(estimate$a2^(-1.5))
  dsig_alpha_1 <- -estimate$a2^(-1.5)/(var(W_m))
  dsig_x_hat <- (W_m-mean(W_m))*dsig_alpha_1
  # An12 <- -mean(x_hat)
  # An22 <- -mean(x_hat^2)
  An211 <- -mean(beta1[1]*da2_x_hat)
  An212 <- -mean(beta1[1]*dsig_x_hat)
  An221 <- mean((Y-beta0-as.matrix(data)%*%beta1)*da2_x_hat-beta1[1]*da2_x_hat*data$x_hat)
  An222 <- mean((Y-beta0-as.matrix(data)%*%beta1)*dsig_x_hat-beta1[1]*dsig_x_hat*data$x_hat)
  if(is.na(Z[1])){
    An2 <- matrix(c(An211,An212,An221,An222),byrow = TRUE,nrow = 2)
    
  }else{
    An231 <- mean(-beta1[1]*da2_x_hat*Z)
    An232 <- mean(-beta1[1]*dsig_x_hat*Z)
    An2 <- matrix(c(An211,An212,An221,An222,An231,An232),byrow = TRUE,nrow = 3)
  }
  An1 <- -t(as.matrix(cbind(1,data)))%*%as.matrix(cbind(1,data))/n
  An411 <- 2*estimate$a2*(mean(W)*estimate$sig_add_2-mean(W2W))
  An412 <- estimate$a2^2*mean(W)-3*mean(W)
  An421 <- mean(Wij)
  An422 <- 1
  # An <- matrix(c(-1,An12,An13,An14,An12,An22,An23,An24,rep(0,2),An33,An34,rep(0,2),An43,An44),byrow = TRUE,nrow = 4)
  An4 <- matrix(c(An411,An412,An421,An422),byrow = TRUE,nrow = 2)
  An3 <- matrix(0,dim(An4)[1],dim(An1)[2])
  An <- rbind(cbind(An1,An2),cbind(An3,An4))
  cov_est <- solve(An)%*%Bn%*%t(solve(An))/n
  return(list(beta1=beta1,beta0=beta0,cov_est=cov_est,alpha0=alpha0,alpha1=alpha1))
} 

# data('bloodpressure')
# bloodpressure
# na_flag<-apply(!is.na(bloodpressure),1,sum)
# bloodpressure <- bloodpressure[-which(na_flag<6),]
# W <- as.matrix(bloodpressure[,c(3:6)])
# est<-estimate_multi(W)

# sbp <- read.delim2("/Users/yuxiangzong/Downloads/Dryad.txt")
# sbp_data <- sbp[,c('Age','SBP_30','SBP_60','SBP_90','SBP_120','Creatinine')]
# na_flag<-apply(!is.na(sbp_data),1,sum)
# sbp_data <- sbp_data[-which(na_flag<6),]
# W <- as.matrix(sbp_data[,2:5])
# hypothesis_triple(W,1)
# hypothesis_triple(W,2)
# hypothesis_triple(W,3)
# W <- as.matrix(sbp_data[,2:5])
# est<-estimate_multi(log(W))
# hypothesis_multi(log(W),1)
# hypothesis_multi(log(W),2)
# hypothesis_multi(log(W),3)

#install.packages("augSIMEX")
# Generic data
library(augSIMEX)
data(GeneRepeat)
data("GeneUni")
data0 <- GeneRepeat$Main
library(ggplot2)
summary(data0[,2:7])
library(tidyr)
data_new <- gather(data0[,2:7],key="rep",value="Totdist")
data_new$time <- rep(1:672,6)
data_new$mean <- rep(apply(data0[,2:7],1,mean),6)
#Statistical Analysis
p <- ggplot(data_new, aes(x=Totdist,color=rep)) + 
  geom_density()
p

p1 <- ggplot(data_new, aes(y=Totdist,x=time,color=rep)) + 
  geom_line()
p1

p1 <- ggplot(data_new, aes(y=Totdist,x=mean,color=rep)) + 
  geom_line()
p1
library(kurtosis)
apply(data0[,2:7],2,mean)
apply(data0[,2:7],2,var)
apply(data0[,2:7],2,sd)
apply(data0[,2:7],2,kurtosis)
kurtosis

mystats <- function(x,na.omit=FALSE){
  if (na.omit)
    x <- x[!is.na(x)]
  m <- mean(x)
  n <- length(x)
  s <- sd(x)
  var <- var(x)
  skew <- sum((x-m)^3/s^3)/n
  kurt <- sum((x-m)^4/s^4)/n - 3
  return(c(n=n,mean=m,stdev=s,var=var,skew=skew,kurtosis=kurt))
}
round(sapply(data0[,2:7],mystats),3)


for(i in c(2:5)){
  W <-as.matrix(data0[,i:(i+2)])
  est <- estimate_triple(W)
  print(round(c(i,i+1,i+2,est$sig_mul,est$sig_add,hypothesis_triple(W,1),hypothesis_triple(W,2),hypothesis_triple(W,3)),3))
  results_3rep <- rbind(results_3rep,c(round(c(est$sig_mul,est$sig_add,hypothesis_triple(W,1),hypothesis_triple(W,2),hypothesis_triple(W,3)),3)))
}

results_4rep <- data.frame(matrix(ncol = 5,nrow = 0))
colnames(results_4rep) <- c('sigma1','sigma2','HT1','HT2','HT3')

for(i in c(3:4)){
  W <-as.matrix(data0[,i:(i+3)])
  est <- estimate_multi(W)
  print(round(c(i,i+1,i+2,i+3,est$sig_mul,est$sig_add,hypothesis_multi(W,1),hypothesis_multi(W,2),hypothesis_multi(W,3)),3))
  results_4rep <- rbind(results_4rep,c(round(c(est$sig_mul,est$sig_add,hypothesis_multi(W,1),hypothesis_multi(W,2),hypothesis_multi(W,3)),3)))
}

results_5rep <- data.frame(matrix(ncol = 5,nrow = 0))
colnames(results_5rep) <- c('sigma1','sigma2','HT1','HT2','HT3')


for(i in c(2:3)){
  W <-as.matrix(data0[,i:(i+4)])
  est <- estimate_multi(W)
  print(round(c(i,i+1,i+2,i+3,i+4,est$sig_mul,est$sig_add,hypothesis_multi(W,1),hypothesis_multi(W,2),hypothesis_multi(W,3)),3))
  results_5rep <- rbind(results_5rep,c(round(c(est$sig_mul,est$sig_add,hypothesis_multi(W,1),hypothesis_multi(W,2),hypothesis_multi(W,3)),3)))
}

#Regression Analysis
W <-as.matrix(data0[,5:7])
est <- estimate_triple(W)

Y <- data0$rs223979909
RG1 <- estimate_regression(W,data0$Weight,Y)
RG2 <- estimate_regression(W,NA,Y)

fit1 <- lm(Y~apply(W,1,mean))
summary(fit1)
fit2 <- lm(Y~apply(W,1,mean)+data0$Weight)
summary(fit2)
fit3 <- lm(Y~apply(W,1,mean)+data0$Weight+as.factor(data0$Batch_O))
summary(fit3)

###########################################
#Framingham dataset
devtools::install_github("timothyhyndman/deconvolve")
fdata <- deconvolve::framingham
W <- as.matrix(log(fdata[,3:4]-50))
est <- estimate_double(W)
round(c(est$sig_mul,est$sig_add,hypothsis_double(W,1),hypothsis_double(W,2),hypothsis_double(W,3)),3)


W <- as.matrix(log(fdata[,5:6]-50))
est <- estimate_double(W)
round(c(est$sig_mul,est$sig_add,hypothsis_double(W,1),hypothsis_double(W,2),hypothsis_double(W,3)),3)


W <- as.matrix(log(fdata[,8:9]))
est <- estimate_double(W)
round(c(est$sig_mul,est$sig_add,hypothsis_double(W,1),hypothsis_double(W,2),hypothsis_double(W,3)),3)


fdata_new <- data.frame(SBP1=apply(fdata[,3:4],1,mean),SBP2=apply(fdata[,5:6],1,mean))
W <- as.matrix(log(fdata_new-50))
est <- estimate_double(W)
round(c(est$sig_mul,est$sig_add,hypothsis_double(W,1),hypothsis_double(W,2),hypothsis_double(W,3)),3)

