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
  
  Wi2jk <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  #for (k in 1:n){
  for (i in c(1:rp)){
    for (j in c(1:rp)[-i]){
      for (m in c(1:rp)[-c(j,i)]){
        Wi2jk[,k] <-(W2[,i])*W[,j]*W[,m] 
        k=k+1
      }
    }
  }
  
  
  #matrix of  w1^3*W2^2*W3
  
  Wi3j2k <- matrix(nrow=n,ncol=6)
  
  k=1
  #for (k in 1:n){
  for (i in c(1:rp)){
    for (j in c(1:rp)[-i]){
      for (m in c(1:rp)[-c(j,i)]){
        Wi3j2k[,k] <-(W3[,i])*W[,j]*W[,m] 
        k=k+1
      }
    }
  }
  
  #matrix of  w1^2*W2^2*W3
  Wi2j2k <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  
  for (i in c(1:rp)){
    for (j in c(1:rp)[-i]){
      for (m in c(1:rp)[-c(j,i)]){
        Wi2j2k[,k] <-(W2[,i])*W2[,j]*W[,m] 
        k=k+1
      }
    }
  }
  
  
  #matrix of  w1^3*W2*W3
  Wi3jk <- matrix(nrow=n,ncol=rp*(rp-1))
  
  k=1
  #for (k in 1:n){
  for (i in c(1:rp)){
    for (j in c(1:rp)[-i]){
      for (m in c(1:rp)[-c(j,i)]){
        Wi3jk[,k] <-(W3[,i])*W[,j]*W[,m] 
        k=k+1
      }
    }
  }
  #}
  
  
  # estimate of a2
  a2 <- (mean(W2)*mean(W)-mean(W2W))/(mean(Wij)*mean(W)-mean(apply(W,1,prod)))
  
  # estimate of sig_mul
  sig_mul <- sqrt(2*log(sqrt(a2)))
  
  # estimate of sig_add_2^2
  sig_add_2 <- mean(W2)-a2*mean(Wij)
  
  # estimate of sig_add
  sig_add <- sqrt(mean(W2)-a2*mean(Wij))
  
  result <- list(a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add)
  return(result)
}

Diff_Simulation_Triple <- function(command,sig_mul_vec,sig_add_vec,N){
  num_mul <- length(sig_mul_vec)
  num_add <- length(sig_add_vec)
  
  a2_mean_m <- matrix(nrow = num_mul,ncol = num_add)
  a2_sd_m <- matrix(nrow = num_mul,ncol = num_add)
  sig_mul_mean_m <- matrix(nrow = num_mul,ncol = num_add)
  sig_add_mean_m <- matrix(nrow = num_mul,ncol = num_add)
  sig_mul_sd_m <- matrix(nrow = num_mul,ncol = num_add)
  sig_add_sd_m <- matrix(nrow = num_mul,ncol = num_add)
  a2_bias_m <- matrix(nrow = num_mul,ncol = num_add)
  a2_mse_m <- matrix(nrow = num_mul,ncol = num_add)
  sig_mul_bias_m <- sig_mul_mse_m <- sig_add_bias_m <- sig_add_mse_m <-matrix(nrow = num_mul,ncol = num_add)
  
  for(i in c(1:num_mul)){
    for(j in c(1:num_add)){
      a2_vec <- c()
      sig_mul_rv <- c()
      sig_add_rv <- c()
      
      for(k in c(1:N)){
        eval(parse(text=command))
        X <- matrix(rep(X,3),ncol=3)
        n <- dim(X)[1]
        rp <- dim(X)[2]
        W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec[i]),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec[j]),n,rp)
        model <- estimate_triple(W)
        if(is.na(model$sig_mul)) model$sig_mul <- 0
        if(is.na(model$sig_add)) model$sig_add <- 0
        a2_vec <- c(a2_vec,model$a2)
        sig_mul_rv <- c(sig_mul_rv,model$sig_mul)
        sig_add_rv <- c(sig_add_rv,model$sig_add)
      }
      a2_mean_m[i,j] <- mean(a2_vec)
      a2_bias_m[i,j] <- mean((a2_vec-exp(sig_mul_vec[i]^2/2)))
      a2_sd_m[i,j] <- sd(a2_vec)
      a2_mse_m[i,j] <- mean((a2_vec-exp(sig_mul_vec[i]^2/2))^2)
      sig_mul_mean_m[i,j] <- mean(sig_mul_rv)
      sig_mul_bias_m[i,j] <- mean(sig_mul_rv-sig_mul_vec[i])
      sig_mul_sd_m[i,j] <- sd(sig_mul_rv)
      sig_mul_mse_m[i,j] <- mean((sig_mul_rv-sig_mul_vec[i])^2)
      sig_add_mean_m[i,j] <- mean(sig_add_rv)
      sig_add_bias_m[i,j] <- mean(sig_add_rv-sig_add_vec[j])
      sig_add_sd_m[i,j] <- sd(sig_add_rv)
      sig_add_mse_m[i,j] <- mean((sig_add_rv-sig_add_vec[j])^2)
    }
  }
  result <- data.frame(sig_mul = rep(sig_mul_vec,each=num_add), sig_add = rep(sig_add_vec,num_mul),
                       #a2_bias=c(t(a2_bias_m)), a2_sd = c(t(a2_sd_m)),  
                       #a2_mse = c(t(a2_mse_m)),
                       sig_mul_bias = c(t(sig_mul_bias_m)),
                       sig_mul_RB = c(t(sig_mul_bias_m))/rep(sig_mul_vec,each=num_add),
                       sig_mul_sd = c(t(sig_mul_sd_m)),
                       sig_mul_mse = c(t(sig_mul_mse_m)),sig_add_bias = c(t(sig_add_bias_m)),sig_add_RB = c(t(sig_add_bias_m))/rep(sig_add_vec,num_mul),
                       sig_add_sd = c(t(sig_add_sd_m)),
                       sig_add_mse = c(t(sig_add_mse_m)))
}

name_vec <- c()
sd_vec <- c(sqrt(exp(0.25)*(exp(0.25)-1)),sqrt(exp(1)*(exp(1)-1)),2,sqrt(8),2,sqrt(6),2*sqrt(3*4/((3+4)^2*(3+4+1))),2*sqrt(1*2/((1+2)^2*(1+2+1))))
N <- 500
for(i in c(300,500,1000)){
  dis_vec <- c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rlnorm(",i,",0,1)"),
               paste0("X <- rexp(",i,",0.5)"),paste0("X <- rgamma(",i,",2,1/2)"),
               paste0("X <- rchisq(",i,",2)"),paste0("X <- rchisq(",i,",3)"),
               paste0("X <- 2*rbeta(",i,",3,4)+1"),paste0("X <- 2*rbeta(",i,",1,2)+1"))
  #dis_vec <- c(paste0("X <- 2*rbeta(",i,",3,4)+1"),paste0("3*rbeta(",i,",2,1)+2"))
  for(j in 1:length(dis_vec)){
    command <- dis_vec[j]
    eval(parse(text=command))
    sig_add_vec <- sd_vec[j]*seq(0,0.75,0.25)
    sig_mul_vec <- sd_vec[j]*seq(0,0.75,0.25)
    #var_xe <- mean(X^2)*(exp(2*sig_mul_vec^2))-mean(X)^2*exp(sig_mul_vec^2)
    #sig_add_vec <- sqrt(var_xe)*seq(0,0.8,0.2)
    name <- paste0("simu_ET_",j,"_",i)
    name_vec <- c(name_vec,name)
    result <- Diff_Simulation_Triple(command,sig_mul_vec,sig_add_vec,N)
    assign(name,result)
    print(c(i,j))
  }
}

table1 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1:4)],collapse=","),")"))),3)
table2 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(5:8)],collapse=","),")"))),3)
table3 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(9:12)],collapse=","),")"))),3)
table4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(13:16)],collapse=","),")"))),3)
table5 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(17:20)],collapse=","),")"))),3)
table6 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(21:24)],collapse=","),")"))),3)

dis_name <- c("LogNorm(0,0.5)","LogNorm(0,1)","Exp(0.5)","Gamma(2,1/2)","\\chi_2^2","\\chi_3^2","2*Beta(3,4)+1","2*Beta(1,2)+1")
dis <- c()
for(i in 1:length(dis_name)){dis <- c(dis,c(dis_name[i],rep(" ",15)))}
table1 <- cbind(data.frame(dis_name = dis[1:64]),table1)
table2 <- cbind(data.frame(dis_name = dis[65:128]),table2)
table3 <- cbind(data.frame(dis_name = dis[1:64]),table3)
table4 <- cbind(data.frame(dis_name = dis[65:128]),table4)
table5 <- cbind(data.frame(dis_name = dis[1:64]),table5)
table6 <- cbind(data.frame(dis_name = dis[65:128]),table6)

xtable1 <- xtable(table1,digits = 3)
xtable2 <- xtable(table2,digits = 3)
xtable3 <- xtable(table3,digits = 3)
xtable4 <- xtable(table4,digits = 3)
xtable5 <- xtable(table5,digits = 3)
xtable6 <- xtable(table6,digits = 3)

print(xtable1, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable2, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable3, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable4, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable5, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable6, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

