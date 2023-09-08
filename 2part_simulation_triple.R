options(scipen = 200,digits = 4, nsmall=3)
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
  Wijk<-apply(W,1,prod)
  #Wijk <- matrix(nrow=n,ncol=1)
  
  #k=2
  #W[k,1]
  #for(k in 1:n){
  #     Wijk[k,1] <- W[k,1]*W[k,2]*W[k,3]      #k <- k+1
  # }
  
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
  
  result <- list(a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add, EWW=mean(Wij), EW=mean(W), W2=W2, W2W=W2W,n=n,Wij=Wij,W=W,Wijk=Wijk)
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
  
  innerfunction <- function(a2, sig_add, EW, EWW, W2, W2W, Wij, n, W, Wijk,W3, W4, W3W, W3W2, W2W2, W4W, Wi2jk, Wi3jk,Wi2j2k,Wi3j2k){
    
    #####################
    #change the variance#
    #####################
    
    g1_i <- sig_add^2*apply(W, 1, mean) + a2*Wijk - apply(W2W,1,mean)
    
    g2_i <- sig_add^2 - apply(W2,1,mean) + a2*apply(Wij,1,mean)
    
    v1 <- cov(data.frame(g1=g1_i,g2=g2_i))
    
    v<-matrix(nrow=2, ncol=2)
    
    
    #v[1,1]<-a2^2*apply(Wijk2,2,mean)+a2*sig_add^2*mean(Wi2jk)-a2*mean(Wi3j2k)+a2*sig_add^2*apply(Wi2jk,2,mean)+sig_add^4*mean(W2)-sig_add^2*mean(W3W)-a2*apply(Wi3jk,2,mean)-sig_add^2*mean(W3W)+mean(W4W)
    v[1,1]<-a2^2*mean(Wijk2)+a2*sig_add^2*mean(Wi2jk)-a2*mean(Wi3j2k)+a2*sig_add^2*mean(Wi2jk)+sig_add^4*mean(W2)-sig_add^2*mean(W3W)-a2*mean(Wi3jk)-sig_add^2*mean(W3W)+mean(W4W)
    v[1,2]<-mean(a2*sig_add^2*Wijk)-mean(a2*Wi3jk)+mean(a2^2*Wi2j2k)+mean(sig_add^4*W-sig_add^2*W3)+mean((sig_add^2*a2-sig_add^2)*W2W)-mean(W4W+a2*W3W2)
    v[2,1]<-v[1,2]
    v[2,2]<-sig_add^4-sig_add^2*mean(W2)+a2*sig_add^2*mean(Wij)-sig_add^2*mean(W2)+mean(W4)-mean(a2*W3W)+mean(a2*sig_add^2*Wij)+mean(a2*W3W)+a2^2*mean(W2W2)
    
    v2_2 <- c(mean(g1_i),mean(g2_i))%*%t(c(mean(g1_i),mean(g2_i)))
    
    v2 <- v-v2_2
    
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
    
    
    cov_est <-innerfunction(a2,sig_add,EW, EWW, W2, W2W, Wij, n, W, Wijk,W3, W4, W3W, W3W2, W2W2, W4W, Wi2jk, Wi3jk,Wi2j2k,Wi3j2k)
    d2 <- t(c(a2_est-1, sig_add_2_est))%*%ginv(cov_est)%*%c(a2_est-1, sig_add_2_est)
    p_value <- pchisq(d2,2,lower.tail = FALSE)/4
    
  }
  
  if(Type==2){
    
    
    a2 <- 1
    
    sig_add <- sqrt(sig_add_2_est)
    
    
    cov_est <-innerfunction(a2,sig_add,EW, EWW, W2, W2W, Wij, n, W, Wijk,W3, W4, W3W, W3W2, W2W2, W4W, Wi2jk, Wi3jk,Wi2j2k,Wi3j2k)
    p_value <- 1 - pnorm(a2_est,1,sqrt(cov_est[1,1]))
    #p_value <- 1 - phalfnorm((sqrt(a2_est)-1)/sqrt(cov_est[1,1]))
  }
  if(Type==3){
    a2 <- a2_est
    sig_add <-0
    cov_est <-innerfunction(a2,sig_add,EW, EWW, W2, W2W, Wij, n, W, Wijk,W3, W4, W3W, W3W2, W2W2, W4W, Wi2jk, Wi3jk,Wi2j2k,Wi3j2k)
    p_value <- 1 - pnorm(sig_add_2_est/sqrt(cov_est[2,2]))
  }
  
  return(p_value)
}

HP_Simulation_Triple <- function(command,sig_mul_vec,sig_add_vec,N){
  # commond: example "X <- rlnorm(500,0,0.5)"
  # sig_mul_vec: the vector of the standard deviations of multiplicative errors
  # sig_add_vec: the vector of the standard deviations of additive errors
  # N: the number of hypothesis tests
  # v: 1 = v1; 2 = v2
  num_mul <- length(sig_mul_vec)
  num_add <- length(sig_add_vec)
  HP_1_m <- HP_2_m <- HP_3_m <-HP_1_m_1 <- HP_2_m_1 <- HP_3_m_1 <- HP_1_m_2 <- HP_2_m_2 <- HP_3_m_2 <- matrix(nrow = num_mul,ncol = num_add)
  
  for(i in c(1:num_mul)){
    for(j in c(1:num_add)){
      
      a2_vec <- c()
      sig_mul_rv <- c()
      sig_add_rv <- c()
      HP_1_rv <- HP_2_rv <- HP_3_rv <- c()
      
      for(k in c(1:N)){
        eval(parse(text=command))
        X <- matrix(rep(X,3),nrow=length(X),3)
        n <- dim(X)[1]
        rp <- dim(X)[2]
        W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec[i]),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec[j]),n,rp)
        HP_1_rv <- c(HP_1_rv,hypothesis_triple(W,1))
        HP_2_rv <- c(HP_2_rv,hypothesis_triple(W,2))
        HP_3_rv <- c(HP_3_rv,hypothesis_triple(W,3))
      }
      
      HP_1_m[i,j] <- sum(HP_1_rv < 0.01)/length(HP_1_rv)
      HP_2_m[i,j] <- sum(HP_2_rv < 0.01)/length(HP_2_rv)
      HP_3_m[i,j] <- sum(HP_3_rv < 0.01)/length(HP_3_rv)
      HP_1_m_1[i,j] <- sum(HP_1_rv < 0.05)/length(HP_1_rv)
      HP_2_m_1[i,j] <- sum(HP_2_rv < 0.05)/length(HP_2_rv)
      HP_3_m_1[i,j] <- sum(HP_3_rv < 0.05)/length(HP_3_rv)
      HP_1_m_2[i,j] <- sum(HP_1_rv < 0.1)/length(HP_1_rv)
      HP_2_m_2[i,j] <- sum(HP_2_rv < 0.1)/length(HP_2_rv)
      HP_3_m_2[i,j] <- sum(HP_3_rv < 0.1)/length(HP_3_rv)
      print(c(i,j))
    }
  }
  result <- data.frame(sig_mul = rep(sig_mul_vec,each=num_add), sig_add = rep(sig_add_vec, num_mul),
                       HP1_0.01 = c(t(HP_1_m)),HP2_0.01 = c(t(HP_2_m)),HP3_0.01 = c(t(HP_3_m)),
                       HP1_0.05 = c(t(HP_1_m_1)),HP2_0.05 = c(t(HP_2_m_1)),HP3_0.05 = c(t(HP_3_m_1)),
                       HP1_0.1 = c(t(HP_1_m_2)),HP2_0.1 = c(t(HP_2_m_2)),HP3_0.1 = c(t(HP_3_m_2)))
}

name_vec <- c()
N <- 500
sd_vec <- c(sqrt(exp(0.25)*(exp(0.25)-1)),sqrt(exp(1)*(exp(1)-1)),2,sqrt(6))
for(i in c(300,500,1000)){
  dis_vec <- c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rlnorm(",i,",0,1)"),
               paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))
  #dis_vec <- c(paste0("X <- rlnorm(",i,",0,0.5)"))
  for(j in 1:length(dis_vec)){
    command <- dis_vec[j]
    eval(parse(text=command))
    sig_add_vec <- sd_vec[j]*seq(0,0.75,0.25)
    sig_mul_vec <- sd_vec[j]*seq(0,0.75,0.25)
    #var_xe <- mean(X^2)*(exp(2*sig_mul_vec^2))-mean(X)^2*exp(sig_mul_vec^2)
    #sig_add_vec <- sqrt(var_xe)*seq(0,0.8,0.2)
    name <- paste0("simu_HPT_",substr(command,7,9),"_",j,"_",i)
    name_vec <- c(name_vec,name)
    result <- HP_Simulation_Triple(command,sig_mul_vec,sig_add_vec,N)
    assign(name,result)
    print(c(i,j))
  }
}

table1 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1:4)],collapse=","),")"))),3)
table2 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(5:8)],collapse=","),")"))),3)
table3 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(9:12)],collapse=","),")"))),3)

dis_name <- c("LogNorm(0,0.5)","LogNorm(0,1)","Exp(0.5)","\\chi_3^2")
dis <- c()
for(i in 1:length(dis_name)){dis <- c(dis,c(dis_name[i],rep(" ",15)))}
table1 <- cbind(data.frame(dis_name = dis),table1)
table2 <- cbind(data.frame(dis_name = dis),table2)
table3 <- cbind(data.frame(dis_name = dis),table3)

xtable1 <- xtable(table1,digits = 3)
xtable2 <- xtable(table2,digits = 3)
xtable3 <- xtable(table3,digits = 3)

print(xtable1, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable2, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable3, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
