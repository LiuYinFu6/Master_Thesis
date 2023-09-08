estimate_double <- function(W){
  # X: input real value
  # sig_mul0: real standard deviation of multiplicative error
  # sig_add0: real standard deviation of additive error
  n = dim(W)[1]
  rp = dim(W)[2]
  
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
  
  # parameters of cubic function 
  # A*x^3 + B*x^2 + C*x + D 
  A <- mean(W)*mean(Wij)
  
  B <- mean(W2W) - mean(W2)*mean(W)
  
  C <- -3*mean(W)*mean(Wij)
  
  D <- 3*mean(W)*mean(W2)-mean(W3)
  
  #discriminant delta
  delta <- (4*(B^2-3*A*C)^3 - (2*B^3-9*A*B*C+27*A^2*D)^2)/(27*A^2)
  
  A <- as.complex(A)
  
  B <- as.complex(B)
  
  C <- as.complex(C)
  
  D <- as.complex(D)
  
  # calculate the roots through analytical way
  x1 <- (((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3) - (- B^2/(9*A^2) + C/(3*A))/(((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3) - B/(3*A)
  
  x2 <- (- B^2/(9*A^2) + C/(3*A))/(2*(((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3)) - B/(3*A) - (3^(1/2)*((- B^2/(9*A^2) + C/(3*A))/(((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3) + (((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3))*1i)/2 - (((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3)/2
  
  x3 <- (3^(1/2)*((- B^2/(9*A^2) + C/(3*A))/(((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3) + (((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3))*1i)/2 - B/(3*A) + (- B^2/(9*A^2) + C/(3*A))/(2*(((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3)) - (((D/(2*A) + B^3/(27*A^3) - (B*C)/(6*A^2))^2 + (- B^2/(9*A^2) + C/(3*A))^3)^(1/2) - D/(2*A) - B^3/(27*A^3) + (B*C)/(6*A^2))^(1/3)/2
  
  #root0 <- polyroot(z = c(D,C,B,A))
  
  roots <- c(x1,x2,x3)
  
  reals <- sapply(roots, function(i) isTRUE(all.equal(Im(i),0)))
  
  a2 <- max(Re(roots[reals]))
  
  #a2 <- max(Re(root))
  
  sig_mul <- sqrt(2*log(sqrt(a2)))
  
  sig_add_2 <- mean(W2)-a2*mean(Wij)
  
  sig_add <- sqrt(mean(W2)-a2*mean(Wij))
  
  mu1 <- mean(W)/sqrt(a2)
  
  mu2 <- (mean(Wij))/(a2)
  
  result <- list(A=A,B=B,C=C,D=D,delta=delta,roots=roots,a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add,mu1=mu1,mu2=mu2)
  
  return(result)
}

#hypothesis test function
hypothsis_double <- function(W,Type,v){
  result_est <- estimate_double(W)
  
  a2_est <- result_est$a2
  
  sig_add_2_est <- result_est$sig_add_2
  
  # X: input real value
  # sig_mul0: real standard deviation of multiplicative error
  # sig_add0: real standard deviation of additive error
  n = dim(W)[1]
  rp = dim(W)[2]
  
  if(a2_est<1){a2_est <- 1}
  
  if(sig_add_2_est<0) sig_add_2_est<-0
  
  innerfunction <- function(a2,sig_add,W,v){
    n = dim(W)[1]
    rp = dim(W)[2]
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
    
    EWW <- mean(Wij)
    
    EW <- mean(W)
    
    W4 <- W^4
    
    W5 <- W^5
    
    W6 <- W^6
    
    W4W2 <- matrix(nrow=n,ncol=rp*(rp-1))
    
    k=1
    for(i in c(1:(rp))){
      for(j in c(1:rp)[-i]){
        W4W2[,k] <- W4[,i]*W2[,j]
        k=k+1
      }
    }
    
    W5W1 <- matrix(nrow=n,ncol=rp*(rp-1))
    
    k=1
    for(i in c(1:(rp))){
      for(j in c(1:rp)[-i]){
        W5W1[,k] <- W5[,i]*W[,j]
        k=k+1
      }
    }
    
    W4W1 <- matrix(nrow=n,ncol=rp*(rp-1))
    
    k=1
    for(i in c(1:(rp))){
      for(j in c(1:rp)[-i]){
        W4W1[,k] <- W4[,i]*W[,j]
        k=k+1
      }
    }
    
    W3W2 <- matrix(nrow=n,ncol=rp*(rp-1))
    
    k=1
    for(i in c(1:(rp))){
      for(j in c(1:rp)[-i]){
        W3W2[,k] <- W3[,i]*W2[,j]
        k=k+1
      }
    }
    
    W3W1 <- matrix(nrow=n,ncol=rp*(rp-1))
    
    k=1
    for(i in c(1:(rp))){
      for(j in c(1:rp)[-i]){
        W3W1[,k] <- W3[,i]*W[,j]
        k=k+1
      }
    }
    
    W2W2 <- matrix(nrow=n,ncol=rp*(rp-1))
    
    k=1
    for(i in c(1:(rp))){
      for(j in c(1:rp)[-i]){
        W2W2[,k] <- W2[,i]*W2[,j]
        k=k+1
      }
    }
    
    g1_i <- a2^2*apply(W,1,mean)*sig_add^2 - 3*apply(W,1,mean)*sig_add^2 - a2^2*apply(W2W,1,mean) + apply(W3,1,mean)
    
    g2_i <- sig_add^2 - apply(W2,1,mean) + a2*apply(Wij,1,mean)
    
    # v1
    v1 <- cov(data.frame(g1=g1_i,g2=g2_i))
    
    # v2
    v11 <- mean(W6 -2*a2^2*W5W1 + a2^4*W4W2 + (2*sig_add^2*a2^2 - 6*sig_add^2)*W4 + 2*a2^2*(- sig_add^2*a2^2 + 3*sig_add^2)*W3W1 + (- sig_add^2*a2^2 + 3*sig_add^2)^2*W2)
    
    v12 <- mean(- W5 + (a2^2 + a2)*W4W1 - a2^3*W3W2 + (- sig_add^2*a2^2 + 4*sig_add^2)*W3 + (- a2^2*sig_add^2*W2W) - (a2*(3*sig_add^2 - a2^2*sig_add^2))*W2W - sig_add^2*(3*sig_add^2 - a2^2*sig_add^2)*W)
    
    v22 <- mean(W4 -2*a2*W3W1 + a2^2*W2W2 - 2*sig_add^2*W2) + mean(2*a2*sig_add^2*Wij) + sig_add^4
    
    v2_1 <- matrix(c(v11,v12,v12,v22),nrow = 2,byrow = TRUE)
    
    v2_2 <- c(mean(g1_i),mean(g2_i))%*%t(c(mean(g1_i),mean(g2_i)))
    
    v2 <- v2_1-v2_2
    
    if(v==1) V <- v1
    else if(v==2) V <- v2
    else return("error")
    
    gamma <- matrix(c(2*a2*sig_add^2*EW-2*a2*mean(W2W),a2^2*EW-3*EW,EWW,1),nrow = 2,byrow = TRUE)
    
    cov_est <- solve(gamma)%*%V%*%t(solve(gamma))/n
    
    return(cov_est)
  }
  
  
  # Type 1 HP
  if(Type==1){
    
    a2 <- 1
    
    sig_add <- 0
    
    cov_est <- innerfunction(a2,sig_add,W,v)
    
    #d2 <- t(c(a2_est-1, sig_add_2_est))%*%ginv(cov_est)%*%c(a2_est-1, sig_add_2_est)
    
    #if(is.na(d2)){
    #d2=0
    #}
    
    #eigens <- eigen(cov_est)
    
    d2 <- t(c(a2_est-1, sig_add_2_est))%*%ginv(cov_est)%*%c(a2_est-1, sig_add_2_est)
    
    p_value <- pchisq(d2,2,lower.tail = FALSE)/4 + ifelse(d2==0,1,0)*0.75
    
    # p value calculated based on 1/4 multivariate normal distribution
    #p_value <- 1 - ptmvnorm(c(chol(ginv(cov_est),pivot=TRUE)%*%c(a2_est-1, sig_add_2_est)),mu=c(0,0),diag(1,2),c(0,0),c(Inf,Inf))
    
    #p_value <- 1 - ptmvnorm(c(a2_est-1, sig_add_2_est),mu=c(0,0),cov_est,c(0,0),c(Inf,Inf))
  }
  
  else if(Type==2){
    
    a2 <- 1
    
    sig_add <- sqrt(sig_add_2_est)
    
    cov_est<-innerfunction(a2,sig_add,W,v)
    
    p_value <- 1 - pnorm(a2_est,1,sqrt(cov_est[1,1])) + ifelse(a2_est==1,1,0)*0.5
    
    # p-value calculated based on half normal distribution
    
    #p_value <- 1 - phalfnorm((a2_est-1)/sqrt(cov_est[1,1]))
    
    #d2 <- t(c(a2_est-1, sig_add_est))%*%solve(cov_est)%*%c(a2_est-1, sig_add_est)
    
    #p_value <- 1-pchisq(d2,1)
  }
  else{
    
    a2 <- a2_est
    
    sig_add <- 0
    
    cov_est<-innerfunction(a2,sig_add,W,v)
    
    p_value <- 1 - pnorm(sig_add_2_est/sqrt(cov_est[2,2])) + ifelse(sig_add_2_est==0,1,0)*0.5
    
    #d2 <- t(c(a2_est-1, sig_add_est))%*%solve(cov_est)%*%c(a2_est-1, sig_add_est)
    
    #p_value <- 1-pchisq(d2,1)
    
    # p-value calculated based on half normal distribution
    
    #p_value <- 1 - phalfnorm(sig_add_2_est/sqrt(cov_est[2,2]))
  }
  return(p_value)
}

#Simulation Function
#HP_Simulation_Double uses same X when doing N hypothesis tests
HP_Simulation_Double <- function(X,sig_mul_vec,sig_add_vec,N,v){
  # X: real data matrix
  # sig_mul_vec: the vector of the standard deviations of multiplicative errors
  # sig_add_vec: the vector of the standard deviations of additive errors
  # N: the number of hypothesis tests
  num_mul <- length(sig_mul_vec)
  num_add <- length(sig_add_vec)
  HP_1_m <- HP_2_m <- HP_3_m <- matrix(nrow = num_mul,ncol = num_add)
  
  for(i in c(1:num_mul)){
    for(j in c(1:num_add)){
      n <- dim(X)[1]
      rp <- dim(X)[2]
      a2_vec <- c()
      sig_mul_rv <- c()
      sig_add_rv <- c()
      HP_1_rv <- HP_2_rv <- HP_3_rv <- c()
      
      for(k in c(1:N)){
        W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec[i]),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec[j]),n,rp)
        HP_1_rv <- c(HP_1_rv,hypothsis_double(W,1,v))
        HP_2_rv <- c(HP_2_rv,hypothsis_double(W,2,v))
        HP_3_rv <- c(HP_3_rv,hypothsis_double(W,3,v))
      }
      
      HP_1_m[i,j] <- sum(HP_1_rv < 0.05)/length(HP_1_rv)
      HP_2_m[i,j] <- sum(HP_2_rv < 0.05)/length(HP_2_rv)
      HP_3_m[i,j] <- sum(HP_3_rv < 0.05)/length(HP_3_rv)
      print(i,j)
    }
  }
  result <- list(HP_1_m=HP_1_m,HP_2_m=HP_2_m,HP_3_m=HP_3_m)
}


#HP_Simulation_Double_2 generates different X when doing N hypothesis tests
HP_Simulation_Double_2 <- function(command,sig_mul_vec,sig_add_vec,N,v){
  # commond: example "X <- rlnorm(500,0,0.5)"
  # sig_mul_vec: the vector of the standard deviations of multiplicative errors
  # sig_add_vec: the vector of the standard deviations of additive errors
  # N: the number of hypothesis tests
  # v: 1 = v1; 2 = v2
  num_mul <- length(sig_mul_vec)
  num_add <- length(sig_add_vec)
  HP_1_m <- HP_2_m <- HP_3_m <- matrix(nrow = num_mul,ncol = num_add)
  
  for(i in c(1:num_mul)){
    for(j in c(1:num_add)){
      
      a2_vec <- c()
      sig_mul_rv <- c()
      sig_add_rv <- c()
      HP_1_rv <- HP_2_rv <- HP_3_rv <- c()
      
      for(k in c(1:N)){
        eval(parse(text=command))
        X <- matrix(rep(X,2),nrow=length(X),2)
        n <- dim(X)[1]
        rp <- dim(X)[2]
        W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec[i]),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec[j]),n,rp)
        HP_1_rv <- c(HP_1_rv,hypothsis_double(W,1,v))
        HP_2_rv <- c(HP_2_rv,hypothsis_double(W,2,v))
        HP_3_rv <- c(HP_3_rv,hypothsis_double(W,3,v))
      }
      
      HP_1_m[i,j] <- sum(HP_1_rv < 0.05)/length(HP_1_rv)
      HP_2_m[i,j] <- sum(HP_2_rv < 0.05)/length(HP_2_rv)
      HP_3_m[i,j] <- sum(HP_3_rv < 0.05)/length(HP_3_rv)
      print(i,j)
    }
  }
  result <- list(HP_1_m=HP_1_m,HP_2_m=HP_2_m,HP_3_m=HP_3_m)
}


###################################
# Conduct Simulation 
# v: v1, v2
# i: sample size 500, 2000
# j: lognormal distribution, exponential distribution
N <- 100
for(v in c(1,2)){
  for(i in c(500,2000)){
    #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
      eval(parse(text=j))
      sig_add_vec <- sd(X)*seq(0,0.8,0.1)
      sig_mul_vec <- sd(X)*seq(0,0.8,0.1)
      name <- paste0("simu_HP_",substr(j,7,9),"_v_",v,"_",i)
      result <- HP_Simulation_Double_2(j,sig_mul_vec,sig_add_vec,N,v)
      row.names(result$HP_1_m) <- round(sig_mul_vec,digits = 3)
      colnames(result$HP_1_m) <- round(sig_add_vec,digits = 3) 
      row.names(result$HP_2_m) <- round(sig_mul_vec,digits = 3)
      colnames(result$HP_2_m) <- round(sig_add_vec,digits = 3)
      row.names(result$HP_3_m) <- round(sig_mul_vec,digits = 3)
      colnames(result$HP_3_m) <- round(sig_add_vec,digits = 3)
      assign(name,result)
    }
  }
}