library(Matrix)
library(xtable)


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

estimate_regression_real <- function(W,Z,Y,sig_mul,sig_add){
  n <- dim(W)[1]
  rp = dim(W)[2]
  Wij <- matrix(nrow=n,ncol=rp*(rp-1)/2)
  
  k=1
  for(i in c(1:(rp-1))){
    for(j in c((i+1):rp)){
      Wij[,k] <- W[,i]*W[,j]
      k <- k+1
    }
  }
  estimate <- list()
  W_m <- apply(W,1,mean)
  estimate$a2 <- exp(sig_mul^2)
  estimate$mu1 <- mean(W_m)/exp(sig_mul^2/2)
  #estimate$mu2 <- (mean(W^2)-sig_add^2)/exp(sig_mul^2*2)
  estimate$mu2 <- (mean(Wij))/exp(sig_mul^2)
  alpha1 <- (sqrt(estimate$a2))*(estimate$mu2-estimate$mu1^2)/(var(W_m))
  alpha0 <- estimate$mu1-mean(W_m)*alpha1
  x_hat <- alpha1*W_m+alpha0
  if(is.na(Z)) data <- data.frame(x_hat) else data <- data.frame(x_hat,Z)
  fit1 <- lm(Y~.,data = data)
  beta1 <- fit1$coefficients[2:length(fit1$coefficients)]
  beta0 <- fit1$coefficients[1]
# Bn11 <- mean((Y-beta0-beta1*x_hat)^2)
# Bn12 <- mean((Y-beta0-beta1*x_hat)^2*x_hat)
# Bn22 <- mean((Y-beta0-beta1*x_hat)^2*x_hat^2)
# Bn <- matrix(c(Bn11,Bn12,Bn12,Bn22),byrow = TRUE,nrow = 2)
  phi1 <- Y-beta0-as.matrix(data)%*%beta1
  phi2 <- as.matrix(data)*c(Y-beta0-as.matrix(data)%*%beta1)
  phi <- cbind(phi1,phi2)
  Bn <- t(phi)%*%phi/n
  An <- -t(as.matrix(cbind(1,data)))%*%as.matrix(cbind(1,data))/n
  #An12 <- -mean(x_hat)
  #An22 <- -mean(x_hat^2)
  #An <- matrix(c(-1,An12,An12,An22),byrow = TRUE,nrow = 2)
  cov_est <- solve(An)%*%Bn%*%t(solve(An))/n
  return(list(beta1=beta1,beta0=beta0,cov_est=cov_est,alpha0=alpha0,alpha1=alpha1))
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

Regression_calib_simu <- function(command,N,sig_mul_vec,sig_add_vec,beta0,beta1,sig_e){
  beta1_vec <- c()
  beta0_vec <- c()
  sd_beta1_vec <- c()
  sd_beta0_vec <- c()
  beta1_real_vec <- c()
  beta0_real_vec <- c()
  sd_beta1_real_vec <- c()
  sd_beta0_real_vec <- c()
  beta1_ls_vec <- c()
  beta0_ls_vec <- c()
  sd_beta1_ls_vec <- c()
  sd_beta0_ls_vec <- c()
  for(k in 1:N){
    eval(parse(text=command))
    X <- matrix(rep(X,2),ncol=2)
    n <- dim(X)[1]
    rp <- dim(X)[2]
    W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
    eps <- rnorm(n,0,sig_e)
    y <- beta0 + beta1*X[,1] + eps
    est <- estimate_regression(W,NA,y)
    beta1_vec <- c(beta1_vec,est$beta1)
    beta0_vec <- c(beta0_vec,est$beta0)
    sd_beta1_vec <- c(sd_beta1_vec,sqrt(est$cov_est[2,2]))
    sd_beta0_vec <- c(sd_beta0_vec,sqrt(est$cov_est[1,1]))
    est_real <- estimate_regression_real(W,NA,y,sig_add_vec,sig_mul_vec)
    beta1_real_vec <- c(beta1_real_vec,est_real$beta1)
    beta0_real_vec <- c(beta0_real_vec,est_real$beta0)
    sd_beta1_real_vec <- c(sd_beta1_real_vec,sqrt(est_real$cov_est[2,2]))
    sd_beta0_real_vec <- c(sd_beta0_real_vec,sqrt(est_real$cov_est[1,1]))
    est_ls <- lm(y~apply(W,1,mean))
    est_ls_out <- summary(est_ls)
    beta1_ls_vec <- c(beta1_ls_vec,est_ls$coefficients[2])
    beta0_ls_vec <- c(beta0_ls_vec,est_ls$coefficients[1])
    sd_beta1_ls_vec <- c(sd_beta1_ls_vec,sqrt(est_ls_out$coefficients[2,2]))
    sd_beta0_ls_vec <- c(sd_beta0_ls_vec,sqrt(est_ls_out$coefficients[1,2]))
  }
  result <- data.frame(beta0_ls = c(mean(beta0_ls_vec)-beta0,sd(beta0_ls_vec),mean(sd_beta0_ls_vec),sd(sd_beta0_ls_vec)),
                       beta1_ls = c(mean(beta1_ls_vec)-beta1,sd(beta1_ls_vec),mean(sd_beta1_ls_vec),sd(sd_beta1_ls_vec)),
                       beta0_real = c(mean(beta0_real_vec)-beta0,sd(beta0_real_vec),mean(sd_beta0_real_vec),sd(sd_beta0_real_vec)),
                       beta1_real = c(mean(beta1_real_vec)-beta1,sd(beta1_real_vec),mean(sd_beta1_real_vec),sd(sd_beta1_real_vec)),
                       beta0 = c(mean(beta0_vec)-beta0,sd(beta0_vec),mean(sd_beta0_vec),sd(sd_beta0_vec)),
                       beta1 = c(mean(beta1_vec)-beta1,sd(beta1_vec),mean(sd_beta1_vec),sd(sd_beta1_vec))
  )
}

Regression_calib_Z_simu <- function(commandX,N,commandZ,sig_mul_vec,sig_add_vec,beta0,beta1,betaz,sig_e){
  beta1_vec <- beta0_vec <- betaz_vec <- c()
  sd_beta1_vec <- sd_beta0_vec <- sd_betaz_vec <- c()
  beta1_real_vec <- beta0_real_vec <- betaz_real_vec <- c()
  sd_beta1_real_vec <- sd_beta0_real_vec <- sd_betaz_real_vec <- c()
  beta1_ls_vec <- beta0_ls_vec <- betaz_ls_vec <- c()
  sd_beta1_ls_vec <- sd_beta0_ls_vec <- sd_betaz_ls_vec <- c()
  for(k in 1:N){
    eval(parse(text=commandX))
    eval(parse(text=commandZ))
    X <- matrix(rep(X,2),ncol=2)
    n <- dim(X)[1]
    rp <- dim(X)[2]
    W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
    eps <- rnorm(n,0,sig_e)
    y <- beta0 + beta1*X[,1] + betaz*Z + eps
    est <- estimate_regression(W,Z,y)
    beta1_vec <- c(beta1_vec,est$beta1[1])
    beta0_vec <- c(beta0_vec,est$beta0)
    betaz_vec <- c(betaz_vec,est$beta1[2])
    sd_beta1_vec <- c(sd_beta1_vec,sqrt(est$cov_est[2,2]))
    sd_betaz_vec <- c(sd_betaz_vec,sqrt(est$cov_est[3,3]))
    sd_beta0_vec <- c(sd_beta0_vec,sqrt(est$cov_est[1,1]))
    est_real <- estimate_regression_real(W,Z,y,sig_add_vec,sig_mul_vec)
    beta1_real_vec <- c(beta1_real_vec,est_real$beta1[1])
    beta0_real_vec <- c(beta0_real_vec,est_real$beta0)
    betaz_real_vec <- c(betaz_real_vec,est_real$beta1[2])
    sd_beta1_real_vec <- c(sd_beta1_real_vec,sqrt(est_real$cov_est[2,2]))
    sd_betaz_real_vec <- c(sd_betaz_real_vec,sqrt(est_real$cov_est[3,3]))
    sd_beta0_real_vec <- c(sd_beta0_real_vec,sqrt(est_real$cov_est[1,1]))
    est_ls <- lm(y~apply(W,1,mean)+Z)
    est_ls_out <- summary(est_ls)
    beta1_ls_vec <- c(beta1_ls_vec,est_ls$coefficients[2])
    beta0_ls_vec <- c(beta0_ls_vec,est_ls$coefficients[1])
    betaz_ls_vec <- c(betaz_ls_vec,est_ls$coefficients[3])
    sd_beta1_ls_vec <- c(sd_beta1_ls_vec,sqrt(est_ls_out$coefficients[2,2]))
    sd_beta0_ls_vec <- c(sd_beta0_ls_vec,sqrt(est_ls_out$coefficients[1,2]))
    sd_betaz_ls_vec <- c(sd_betaz_ls_vec,sqrt(est_ls_out$coefficients[3,2]))
  }
  result <- data.frame(beta0_ls = c(mean(beta0_ls_vec)-beta0,sd(beta0_ls_vec),mean(sd_beta0_ls_vec),sd(sd_beta0_ls_vec)),
                       beta1_ls = c(mean(beta1_ls_vec)-beta1,sd(beta1_ls_vec),mean(sd_beta1_ls_vec),sd(sd_beta1_ls_vec)),
                       betaz_ls = c(mean(betaz_ls_vec)-betaz,sd(betaz_ls_vec),mean(sd_betaz_ls_vec),sd(sd_betaz_ls_vec)),
                       beta0_real = c(mean(beta0_real_vec)-beta0,sd(beta0_real_vec),mean(sd_beta0_real_vec),sd(sd_beta0_real_vec)),
                       beta1_real = c(mean(beta1_real_vec)-beta1,sd(beta1_real_vec),mean(sd_beta1_real_vec),sd(sd_beta1_real_vec)),
                       betaz_real = c(mean(betaz_real_vec)-betaz,sd(betaz_real_vec),mean(sd_betaz_real_vec),sd(sd_betaz_real_vec)),
                       beta0 = c(mean(beta0_vec)-beta0,sd(beta0_vec),mean(sd_beta0_vec),sd(sd_beta0_vec)),
                       beta1 = c(mean(beta1_vec)-beta1,sd(beta1_vec),mean(sd_beta1_vec),sd(sd_beta1_vec)),
                       betaz = c(mean(betaz_vec)-betaz,sd(betaz_vec),mean(sd_betaz_vec),sd(sd_betaz_vec))
  )
}


#result <- Regression_calib_Z_simu(paste0(paste0("X <- rexp(",i,",1)")),N,paste0(paste0("Z <- rnorm(",i,",0,1)")),sig_add_vec,sig_add_vec,1,0.8,0.6,0.04)

# write.csv(result_chi,file = "result_chi.csv")
# write.csv(result_exp,file = "result_exp.csv")
# write.csv(result_lno,file = "result_lno.csv")


lab <- data.frame(label=rep(c("Bias","S.D.","Mean($\\widehat{S.E.}$)","S.D.($\\widehat{S.E.}$)"),9))
size <- data.frame(size=rep(c("n=100","","","","n=500","","","","n=2000","","",""),3))
dist <- data.frame(dist=rep(c("LogNorm(0,0.5)",rep("",11),"Exp(0.5)",rep("",11),"Chisq(3)",rep("",11))))

name_vec <- c()
N <- 1000
for(i in c(100,500,2000)){
  for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(k in c(0.04,0.4)){
    #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
    eval(parse(text=j))
    sig_add_vec <- sd(X)*0.3
    sig_mul_vec <- sd(X)*0.3
    name <- paste0("simu_RC_0.3_",substr(j,7,9),"_",i,"_",k)
    name_vec <- c(name_vec,name)
    result <- Regression_calib_simu(j,N,sig_add_vec,sig_add_vec,1,0.8,k)
    assign(name,result)
  }
}
}

result_lno_0.3_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1,7,13)],collapse = ","),")"))),3)
result_lno_0.3_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(2,8,14)],collapse = ","),")"))),3)

result_exp_0.3_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)],collapse = ","),")"))),3)
result_exp_0.3_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)],collapse = ","),")"))),3)

result_chi_0.3_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)+2],collapse = ","),")"))),3)
result_chi_0.3_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)+2],collapse = ","),")"))),3)

result_0.3_0.04 <- rbind(result_lno_0.3_0.04,result_exp_0.3_0.04,result_chi_0.3_0.04)
result_0.3_0.4 <- rbind(result_lno_0.3_0.4,result_exp_0.3_0.4,result_chi_0.3_0.4)

result_0.3_0.04_latex <- xtable(cbind(dist,size,lab,result_0.3_0.04),digits = 3)
result_0.3_0.4_latex <- xtable(cbind(dist,size,lab,result_0.3_0.4),digits = 3)
# 输出LaTeX代码
print(result_0.3_0.04_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(result_0.3_0.4_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

# without z
name_vec <- c()
N <- 1000
for(i in c(100,500,2000)){
  for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(k in c(0.04,0.4)){
      #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
      eval(parse(text=j))
    sig_add_vec <- sd(X)*0.6
    sig_mul_vec <- sd(X)*0.6
    name <- paste0("simu_RC_0.6_",substr(j,7,9),"_",i,"_",k)
    name_vec <- c(name_vec,name)
    result <- Regression_calib_simu(j,N,sig_add_vec,sig_add_vec,1,0.8,k)
    assign(name,result)
    print(c(i,j,k))
  }
}
}

result_lno_0.6_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1,7,13)],collapse = ","),")"))),3)
result_lno_0.6_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(2,8,14)],collapse = ","),")"))),3)

result_exp_0.6_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)],collapse = ","),")"))),3)
result_exp_0.6_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)],collapse = ","),")"))),3)

result_chi_0.6_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)+2],collapse = ","),")"))),3)
result_chi_0.6_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)+2],collapse = ","),")"))),3)

result_0.6_0.04 <- rbind(result_lno_0.6_0.04,result_exp_0.6_0.04,result_chi_0.6_0.04)
result_0.6_0.4 <- rbind(result_lno_0.6_0.4,result_exp_0.6_0.4,result_chi_0.6_0.4)


result_0.6_0.04_latex <- xtable(cbind(dist,size,lab,result_0.6_0.04),digits = 3)
result_0.6_0.4_latex <- xtable(cbind(dist,size,lab,result_0.6_0.4),digits = 3)
# 输出LaTeX代码
print(result_0.6_0.04_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(result_0.6_0.4_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

# result_exp <- round(rbind(simu_RC_0.3_exp_100,simu_RC_0.3_exp_500,simu_RC_0.3_exp_1000,simu_RC_0.3_exp_2000,simu_RC_0.3_exp_5000),3)
# result_chi <- round(rbind(simu_RC_0.3_chi_100,simu_RC_0.3_chi_500,simu_RC_0.3_chi_1000,simu_RC_0.3_chi_2000,simu_RC_0.3_chi_5000),3)
# result_lno <- round(rbind(simu_RC_0.3_lno_100,simu_RC_0.3_lno_500,simu_RC_0.3_lno_1000,simu_RC_0.3_lno_2000,simu_RC_0.3_lno_5000),3)

name_vec <- c()
N <- 1000
for(i in c(100,500,2000)){
  for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(k in c(0.04,0.4)){
      #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
    eval(parse(text=j))
    sig_add_vec <- sd(X)*0.6
    sig_mul_vec <- sd(X)*0.3
    name <- paste0("simu_RC_36_",substr(j,7,9),"_",i,"_",k)
    name_vec <- c(name_vec,name)
    result <- Regression_calib_simu(j,N,sig_add_vec,sig_add_vec,1,0.8,k)
    assign(name,result)
    print(c(i,j,k))
  }
}
}

result_lno_36_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1,7,13)],collapse = ","),")"))),3)
result_lno_36_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(2,8,14)],collapse = ","),")"))),3)

result_exp_36_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)],collapse = ","),")"))),3)
result_exp_36_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)],collapse = ","),")"))),3)

result_chi_36_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)+2],collapse = ","),")"))),3)
result_chi_36_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)+2],collapse = ","),")"))),3)

result_36_0.04 <- rbind(result_lno_36_0.04,result_exp_36_0.04,result_chi_36_0.04)
result_36_0.4 <- rbind(result_lno_36_0.4,result_exp_36_0.4,result_chi_36_0.4)


result_36_0.04_latex <- xtable(cbind(dist,size,lab,result_36_0.04),digits = 3)
result_36_0.4_latex <- xtable(cbind(dist,size,lab,result_36_0.4),digits = 3)

# 输出LaTeX代码
print(result_36_0.04_latex , tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(result_36_0.4_latex , tabular.environment = "longtable", floating = FALSE,include.rownames = F)

################################################
name_vec <- c()
N <- 1000
for(i in c(100,500,1000,2000)){
  for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(k in c(0.04,0.4)){
      #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
      eval(parse(text=j))
      sig_add_vec <- sd(X)*0.6
      sig_mul_vec <- sd(X)*0.6
      name <- paste0("simu_Z_RC_0.6_",substr(j,7,9),"_",i,"_",k)
      name_vec <- c(name_vec,name)
      result <- Regression_calib_Z_simu(j,i,paste0(paste0("Z <- rnorm(",i,",0,1)")),sig_add_vec,sig_add_vec,1,0.8,0.6,k)
      assign(name,result)
      print(c(i,j,k))
    }
  }
}

result_z_lno_0.6_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1,7,13)],collapse = ","),")"))),3)
result_z_lno_0.6_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(2,8,14)],collapse = ","),")"))),3)

result_z_exp_0.6_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)],collapse = ","),")"))),3)
result_z_exp_0.6_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)],collapse = ","),")"))),3)

result_z_chi_0.6_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)+2],collapse = ","),")"))),3)
result_z_chi_0.6_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)+2],collapse = ","),")"))),3)

result_z_0.6_0.04 <- rbind(result_z_lno_0.6_0.04,result_z_exp_0.6_0.04,result_z_chi_0.6_0.04)
result_z_0.6_0.4 <- rbind(result_z_lno_0.6_0.4,result_z_exp_0.6_0.4,result_z_chi_0.6_0.4)

result1 <- cbind(dist,size,lab,result_z_0.6_0.04)
result2 <- cbind(dist,size,lab,result_z_0.6_0.4)

result_z_0.6_0.04_latex <- xtable(cbind(dist,size,lab,result_z_0.6_0.04),digits = 3)
result_z_0.6_0.4_latex <- xtable(cbind(dist,size,lab,result_z_0.6_0.4),digits = 3)
# 输出LaTeX代码
print(result_z_0.6_0.04_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(result_z_0.6_0.4_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

name_vec <- c()
N <- 1000
for(i in c(100,500,2000)){
  for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(k in c(0.04,0.4)){
      #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
      eval(parse(text=j))
      sig_add_vec <- sd(X)*0.3
      sig_mul_vec <- sd(X)*0.3
      name <- paste0("simu_Z_RC_3_",substr(j,7,9),"_",i,"_",k)
      name_vec <- c(name_vec,name)
      result <- Regression_calib_Z_simu(j,i,paste0(paste0("Z <- rnorm(",i,",0,1)")),sig_add_vec,sig_add_vec,1,0.8,0.6,k)
      assign(name,result)
      print(c(i,j,k))
    }
  }
}
result_z_lno_0.3_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1,7,13)],collapse = ","),")"))),3)
result_z_lno_0.3_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(2,8,14)],collapse = ","),")"))),3)

result_z_exp_0.3_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)],collapse = ","),")"))),3)
result_z_exp_0.3_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)],collapse = ","),")"))),3)

result_z_chi_0.3_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)+2],collapse = ","),")"))),3)
result_z_chi_0.3_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)+2],collapse = ","),")"))),3)

result_z_0.3_0.04 <- rbind(result_z_lno_0.3_0.04,result_z_exp_0.3_0.04,result_z_chi_0.3_0.04)
result_z_0.3_0.4 <- rbind(result_z_lno_0.3_0.4,result_z_exp_0.3_0.4,result_z_chi_0.3_0.4)

result_z_0.3_0.04_latex <- xtable(cbind(dist,size,lab,result_z_0.3_0.04),digits = 3)
result_z_0.3_0.4_latex <- xtable(cbind(dist,size,lab,result_z_0.3_0.4),digits = 3)
# 输出LaTeX代码
print(result_z_0.3_0.04_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(result_z_0.3_0.4_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)


name_vec <- c()
N <- 1000
for(i in c(100,500,1000,2000)){
  for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"),paste0("X <- rchisq(",i,",3)"))){
    for(k in c(0.04,0.4)){
      #for(j in c(paste0("X <- rlnorm(",i,",0,0.5)"),paste0("X <- rexp(",i,",0.5)"))){
      eval(parse(text=j))
      sig_add_vec <- sd(X)*0.6
      sig_mul_vec <- sd(X)*0.3
      name <- paste0("simu_Z_RC_36_",substr(j,7,9),"_",i,"_",k)
      name_vec <- c(name_vec,name)
      result <- Regression_calib_Z_simu(j,i,paste0(paste0("Z <- rnorm(",i,",0,1)")),sig_add_vec,sig_add_vec,1,0.8,0.6,k)
      assign(name,result)
      print(c(i,j,k))
    }
  }
}
result_z_lno_36_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(1,7,13)],collapse = ","),")"))),3)
result_z_lno_36_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(2,8,14)],collapse = ","),")"))),3)

result_z_exp_36_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)],collapse = ","),")"))),3)
result_z_exp_36_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)],collapse = ","),")"))),3)

result_z_chi_36_0.04 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(3,9,15)+2],collapse = ","),")"))),3)
result_z_chi_36_0.4 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(4,10,16)+2],collapse = ","),")"))),3)

result_z_36_0.04 <- rbind(result_z_lno_36_0.04,result_z_exp_36_0.04,result_z_chi_36_0.04)
result_z_36_0.4 <- rbind(result_z_lno_36_0.4,result_z_exp_36_0.4,result_z_chi_36_0.4)

result_z_36_0.04_latex <- xtable(cbind(dist,size,lab,result_z_36_0.04),digits = 3)
result_z_36_0.4_latex <- xtable(cbind(dist,size,lab,result_z_36_0.4),digits = 3)

print(result_z_36_0.04_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(result_z_36_0.4_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

X <- rchisq(300,3)
sd(X)
X <- rlnorm(300,0,0.5)
sd(X)