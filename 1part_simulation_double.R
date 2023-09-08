## Section1
options(scipen = 200,digits = 4, nsmall=3)
### Double Simulation
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
  
  result <- list(A=A,B=B,C=C,D=D,delta=delta,roots=roots,a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add)
  
  return(result)
}

Diff_Simulation_Double <- function(command,sig_mul_vec,sig_add_vec,N){
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
        X <- matrix(rep(X,2),ncol=2)
        n <- dim(X)[1]
        rp <- dim(X)[2]
        W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec[i]),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec[j]),n,rp)
        model <- estimate_double(W)
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
  #dis_vec <- c(paste0("X <- rlnorm(",i,",0,1)"))
  for(j in 1:length(dis_vec)){
      command <- dis_vec[j]
      eval(parse(text=command))
      sig_add_vec <- sd_vec[j]*seq(0,0.75,0.25)
      sig_mul_vec <- sd_vec[j]*seq(0,0.75,0.25)
      #var_xe <- mean(X^2)*(exp(2*sig_mul_vec^2))-mean(X)^2*exp(sig_mul_vec^2)
      #sig_add_vec <- sqrt(var_xe)*seq(0,0.8,0.2)
      name <- paste0("simu_ED_",j,"_",i)
      name_vec <- c(name_vec,name)
      result <- Diff_Simulation_Double(command,sig_mul_vec,sig_add_vec,N)
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
table7 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(25:28)],collapse=","),")"))),3)
table8 <- round(eval(parse(text=paste0("rbind(",paste(name_vec[c(29:32)],collapse=","),")"))),3)
dis_name <- c("LogNorm(0,0.5)","LogNorm(0,1)","Exp(0.5)","Gamma(2,1/2)","$\\chisq$(2)","$\\chisq$(3)","2*Beta(3,4)+1","2*Beta(1,2)+1")
dis <- c()
for(i in 1:length(dis_name)){dis <- c(dis,c(dis_name[i],rep(" ",15)))}
table1 <- cbind(data.frame(dis_name = dis[1:64]),table1)
table2 <- cbind(data.frame(dis_name = dis[65:128]),table2)
table3 <- cbind(data.frame(dis_name = dis[1:64]),table3)
table4 <- cbind(data.frame(dis_name = dis[65:128]),table4)
table5 <- cbind(data.frame(dis_name = dis[1:64]),table5)
table6 <- cbind(data.frame(dis_name = dis[65:128]),table6)
table7 <- cbind(data.frame(dis_name = dis[1:64]),table7)
table8 <- cbind(data.frame(dis_name = dis[65:128]),table8)
# cols <- colnames(table1)
# cols <- c(cols[1:8],cols[10],cols[9],cols[11])
# table1 <- table1[,cols]
# table2 <- table2[,cols]
# table3 <- table3[,cols]
# table4 <- table4[,cols]
# table5 <- table5[,cols]
# table6 <- table6[,cols]
# table7 <- table7[,cols]
# table8 <- table8[,cols]
library(xtable)
xtable1 <- xtable(table1,digits = 3)
xtable2 <- xtable(table2,digits = 3)
xtable3 <- xtable(table3,digits = 3)
xtable4 <- xtable(table4,digits = 3)
xtable5 <- xtable(table5,digits = 3)
xtable6 <- xtable(table6,digits = 3)
xtable7 <- xtable(table7,digits = 3)
xtable8 <- xtable(table8,digits = 3)

print(xtable1, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable2, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable3, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable4, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable5, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable6, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable7, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
print(xtable8, tabular.environment = "longtable", floating = FALSE,include.rownames = F)





