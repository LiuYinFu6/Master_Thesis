simex_without_replicate<-function(B,sig_mul_vec,sig_add_vec,W,XI,y1){
  i=length(W)
  #
  #B<-500
  Ub<-matrix(nrow=B,ncol=i)
  for(k in 1:B){
    Ub[k,]<-rnorm(i,0,sig_mul_vec)}
  
  Vb<-matrix(nrow=B,ncol=i)
  for(k in 1:B){
    Vb[k,]<-rnorm(i,0,sig_add_vec)}
  
  l<-length(XI)
  Wbi<-array(0,dim=c(B,i,l))
  m=200
  n=20
  for(k in 1:l){
    for(m in 1:B){
      for(n in 1:i){
        Wbi[m,n,k]<-W[n]*(exp(sqrt(XI[k])*Ub[m,n]))+sqrt(XI[k])*Vb[m,n]*exp(sqrt(XI[k])*Ub[m,n])
        #  Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]
        #   Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]*exp(log(sqrt(XI[k])*Ub[m,n]))
      }
    }
  }
  
  # compute the coefficient of the beta
  Theta_1<-matrix(nrow=B,ncol=l)
  Theta_0<-matrix(nrow=B,ncol=l)
  # Theta_Z<-matrix(nrow=B,ncol=l)
  var_est_1<-matrix(nrow=B,ncol=l)
  var_est_0<-matrix(nrow=B,ncol=l)
  # var_est_Z<-matrix(nrow=B,ncol=l)
  delta_Theta_0<-matrix(nrow=B,ncol=l)
  delta_Theta_1<-matrix(nrow=B,ncol=l)
  delta_Theta_Z<-matrix(nrow=B,ncol=l)
  
  for (l in 1:l){
    for (k in 1:B){
      lm_result<-summary(lm(y1 ~ Wbi[k,,l]))
      Theta_1[k,l]<-lm_result$coefficients[2,1]
      var_est_1[k,l]<-(lm_result$coefficients[2,2])^2
      Theta_0[k,l]<-lm_result$coefficients[1,1]
      var_est_0[k,l]<-(lm_result$coefficients[1,2])^2
      #   Theta_Z<-lm_result$coefficients[3,1]
    }
  }
  
  # for (l in 1:l){
  #  for (k in 1:B){
  #     Theta_0[k,l]<-summary(lm(y1 ~ Wbi[k,,l]))$coefficients[1,1]
  #     var_est_0[k,l]<-(summary(lm(y1 ~ Wbi[k,,l]))$coefficients[1,2])^2
  #   }
  #  }
  Theta_1_mean<- colMeans(Theta_1,na.rm=T)
  var_est_1_mean<-colMeans(var_est_1,na.rm=T)
  Theta_0_mean<- colMeans(Theta_0,na.rm=T)
  var_est_0_mean<-colMeans(var_est_0,na.rm=T)
  
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_1[k,l]<-(Theta_1[k,l]-Theta_1_mean[l])^2
    }
  }
  
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_0[k,l]<-(Theta_0[k,l]-Theta_0_mean[l])^2
    }
  }
  
  #s^2
  delta_Theta_1_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_1_sqr[1,l]<-sum(delta_Theta_1[,l])/(B-1)
    }
  }
  #s^2 for the first coefficient
  delta_Theta_0_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_0_sqr[1,l]<-sum(delta_Theta_0[,l])/(B-1)
    }
  }
  
  #compute the difference
  diff_1=var_est_1_mean-delta_Theta_1_sqr
  simex_var_1<-Exploration(as.vector(diff_1),XI)$beta_simex
  simex_beta_1<-Exploration(as.vector(Theta_1_mean),XI)$beta_simex
  
  #compute the difference
  diff_0=var_est_0_mean-delta_Theta_1_sqr
  simex_var_0<-Exploration(as.vector(diff_0),XI)$beta_simex
  simex_beta_0<-Exploration(as.vector(Theta_0_mean),XI)$beta_simex
  
  #t value
  t_1<-simex_beta_1/sqrt(simex_var_0)
  
  #t value
  t_0<-simex_beta_0/sqrt(simex_var_1)
  
  #p value
  p_1<- 2*pt(q=t_1, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_0,Theta_1_mean)
  
  #p value
  p_0<- 2*pt(q=t_0, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_1,Theta_1_mean)
  # return<-list(simex_beta=simex_beta,ylab=ylab,beta0=beta0,beta1=beta1,beta2=beta2,tvalue=t,pvalue=p)
  # return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,tvalue=t,pvalue=p,Theta_1_mean=Theta_1_mean,diff=diff,simex_var,stderr=sqrt(simex_var))
  return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,Theta_1_mean=Theta_1_mean,simex_beta_1=simex_beta_1,SD_0=sqrt(simex_var_0),SD_1=sqrt(simex_var_1))
}


simex_Z_without_replicate<-function(B,sig_mul_vec,sig_add_vec,W,XI,y1,Z){
  i=length(W)
  #
  #B<-500
  Ub<-matrix(nrow=B,ncol=i)
  for(k in 1:B){
    Ub[k,]<-rnorm(i,0,sig_mul_vec)}
  
  Vb<-matrix(nrow=B,ncol=i)
  for(k in 1:B){
    Vb[k,]<-rnorm(i,0,sig_add_vec)}
  
  l<-length(XI)
  Wbi<-array(0,dim=c(B,i,l))
  m=200
  n=20
  for(k in 1:l){
    for(m in 1:B){
      for(n in 1:i){
        Wbi[m,n,k]<-W[n]*(exp(sqrt(XI[k])*Ub[m,n]))+sqrt(XI[k])*Vb[m,n]*exp(sqrt(XI[k])*Ub[m,n])
        #  Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]
        #   Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]*exp(log(sqrt(XI[k])*Ub[m,n]))
      }
    }
  }
  
  # compute the coefficient of the beta
  Theta_1<-matrix(nrow=B,ncol=l)
  Theta_0<-matrix(nrow=B,ncol=l)
  Theta_Z<-matrix(nrow=B,ncol=l)
  var_est_1<-matrix(nrow=B,ncol=l)
  var_est_0<-matrix(nrow=B,ncol=l)
  var_est_Z<-matrix(nrow=B,ncol=l)
  delta_Theta_0<-matrix(nrow=B,ncol=l)
  delta_Theta_1<-matrix(nrow=B,ncol=l)
  delta_Theta_Z<-matrix(nrow=B,ncol=l)
  
  for (l in 1:l){
    for (k in 1:B){
      lm_result<-summary(lm(y1 ~ Wbi[k,,l]+Z))
      Theta_1[k,l]<-lm_result$coefficients[2,1]
      var_est_1[k,l]<-(lm_result$coefficients[2,2])^2
      Theta_0[k,l]<-lm_result$coefficients[1,1]
      var_est_0[k,l]<-(lm_result$coefficients[1,2])^2
      Theta_Z[k,l]<-lm_result$coefficients[3,1]
      var_est_Z[k,l]<-(lm_result$coefficients[3,2])^2
    }
  }
  
  # for (l in 1:l){
  #  for (k in 1:B){
  #     Theta_0[k,l]<-summary(lm(y1 ~ Wbi[k,,l]))$coefficients[1,1]
  #     var_est_0[k,l]<-(summary(lm(y1 ~ Wbi[k,,l]))$coefficients[1,2])^2
  #   }
  #  }
  Theta_1_mean<- colMeans(Theta_1,na.rm=T)
  var_est_1_mean<-colMeans(var_est_1,na.rm=T)
  Theta_0_mean<- colMeans(Theta_0,na.rm=T)
  var_est_0_mean<-colMeans(var_est_0,na.rm=T)
  Theta_Z_mean<- colMeans(Theta_Z,na.rm=T)
  var_est_Z_mean<-colMeans(var_est_Z,na.rm=T)
  
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_1[k,l]<-(Theta_1[k,l]-Theta_1_mean[l])^2
    }
  }
  
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_0[k,l]<-(Theta_0[k,l]-Theta_0_mean[l])^2
    }
  }
  
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_Z[k,l]<-(Theta_Z[k,l]-Theta_Z_mean[l])^2
    }
  }
  
  #s^2
  delta_Theta_1_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_1_sqr[1,l]<-sum(delta_Theta_1[,l])/(B-1)
    }
  }
  #s^2 for the first coefficient
  delta_Theta_0_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_0_sqr[1,l]<-sum(delta_Theta_0[,l])/(B-1)
    }
  }
  
  delta_Theta_Z_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_Z_sqr[1,l]<-sum(delta_Theta_Z[,l])/(B-1)
    }
  }
  
  
  #compute the difference
  diff_1=var_est_1_mean-delta_Theta_1_sqr
  simex_var_1<-Exploration(as.vector(diff_1),XI)$beta_simex
  simex_beta_1<-Exploration(as.vector(Theta_1_mean),XI)$beta_simex
  
  #compute the difference
  diff_0=var_est_0_mean-delta_Theta_0_sqr
  simex_var_0<-Exploration(as.vector(diff_0),XI)$beta_simex
  simex_beta_0<-Exploration(as.vector(Theta_0_mean),XI)$beta_simex
  
  #compute the difference
  diff_Z=var_est_Z_mean-delta_Theta_Z_sqr
  simex_var_Z<-Exploration(as.vector(diff_Z),XI)$beta_simex
  simex_beta_Z<-Exploration(as.vector(Theta_Z_mean),XI)$beta_simex
  
  #t value
  t_1<-simex_beta_1/sqrt(simex_var_1)
  
  #t value
  t_0<-simex_beta_0/sqrt(simex_var_0)
  
  #t value
  t_Z<-simex_beta_Z/sqrt(simex_var_Z)
  
  #p value
  p_1<- 2*pt(q=t_1, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_0,Theta_1_mean)
  
  #p value
  p_0<- 2*pt(q=t_0, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_1,Theta_1_mean)
  # return<-list(simex_beta=simex_beta,ylab=ylab,beta0=beta0,beta1=beta1,beta2=beta2,tvalue=t,pvalue=p)
  # return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,tvalue=t,pvalue=p,Theta_1_mean=Theta_1_mean,diff=diff,simex_var,stderr=sqrt(simex_var))
  return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,Theta_0_mean=Theta_0_mean,Theta_1_mean=Theta_1_mean,Theta_Z_mean=Theta_Z_mean,simex_beta_Z=simex_beta_Z,SD_0=sqrt(simex_var_0),SD_1=sqrt(simex_var_1),SD_Z=sqrt(simex_var_Z))
}


simex_plot_without_replicate<-function(XI,data,simex_beta){
  XI_new<-c(-1, XI)
  ylab<-c(simex_beta,data)
  beta0<-Exploration(as.vector(data),XI)$beta0
  beta1<-Exploration(as.vector(data),XI)$beta1
  beta2<-Exploration(as.vector(data),XI)$beta2
  Exploration_function <- function(x) beta0+beta1*(x)+beta2*x^2
  data<-data.frame(XI_new,ylab)
  df1<-data.frame(x=seq(-1,0,length.out = 100),
                  y=Exploration_function(seq(-1,0,length.out = 100)))
  df2<-data.frame(x=seq(0,2,length.out = 100),
                  y=Exploration_function(seq(0,2,length.out = 100)))
  p<-ggplot() + geom_point(data=data,aes(x=XI_new,y=ylab))+
    
    # first line
    geom_line(data = df1, aes(x = x, y = y), color = "blue", size = 1.5) +
    
    # second line
    geom_line(data = df2, aes(x = x, y = y), color = "red", size = 1.5, linetype = "dashed")
  return(p)
}

Exploration_without_replicate<- function(data,XI){
  EXmodel<-lm(as.vector(data)~XI+I(XI^2))
  beta0<-summary(EXmodel)$coefficients[1,1]
  beta1<-summary(EXmodel)$coefficients[2,1]
  beta2<-summary(EXmodel)$coefficients[3,1]
  beta_simex<-beta0+beta1*(-1)+beta2*1
  return<-list(beta_simex=beta_simex,beta0=beta0,beta1=beta1,beta2=beta2) 
  # return(beta_simex)
}
Simultation_without_replicate<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set){
  #number_repitition_matrix<-c(10,500,1000)
  sample_number_matrix<-c(100,200,500)
  Bias_matrix<-matrix(nrow=3,ncol=3)
  SD_Mean_matrix<-matrix(nrow=3,ncol=3)
  SD_SE_matrix<-matrix(nrow=3,ncol=3)
  SDE_matrix<-matrix(nrow=3,ncol=3)
  #Full_table<-matrix(nrow=12,ncol=3)
  for(q in 1:3){
    sample_number<-sample_number_matrix[q]
    simexresult_matrix_1<-matrix(nrow=50,ncol=1)
    simexresult_matrix_Z<-matrix(nrow=50,ncol=1)
    simexresult_matrix_0<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_Z_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      XI<-seq(0,2,0.4)
      B=500
      #i=100
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W <- x_real1*U+V
      Z<-rnorm(sample_number,0,1)
      y1 <- 1+0.8*x_real1 +0.6*Z+ rnorm(sample_number, sd = 0.04)
      simexresult<-simex_Z(B,sig_mul_vec,sig_add_vec,W,XI,y1,Z)
      simexresult_matrix_1[index,1]<-simexresult$simex_beta_1
      simexresult_matrix_0[index,1]<-simexresult$simex_beta_0
      simexresult_matrix_Z[index,1]<-simexresult$simex_beta_Z
      SD_0_matrix[index,1]<-simexresult$SD_0
      SD_1_matrix[index,1]<-simexresult$SD_1
      SD_Z_matrix[index,1]<-simexresult$SD_Z
    }
    
    Bias_beta_1<-mean(simexresult_matrix_1)-0.8
    SDE_1<-sd(simexresult_matrix_1)
    SE_1_mean<-mean(SD_1_matrix)
    SE_1_sd<-sd(SD_1_matrix)
    Bias_beta_0<-mean(simexresult_matrix_0)-1
    SDE_0<-sd(simexresult_matrix_0)
    SE_0_mean<-mean(SD_0_matrix)
    SE_0_sd<-sd(SD_0_matrix)
    
    Bias_beta_Z<-mean(simexresult_matrix_Z)-0.6
    SDE_Z<-sd(simexresult_matrix_Z)
    SE_Z_mean<-mean(SD_Z_matrix)
    SE_Z_sd<-sd(SD_Z_matrix)
    
    
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    Bias_matrix[q,3]<-Bias_beta_Z
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SDE_matrix[q,3]<-SDE_Z
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_Mean_matrix[q,3]<-SE_Z_mean
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    SD_SE_matrix[q,3]<-SE_Z_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
    
    
    
  }
  Full_table<-matrix(nrow=12,ncol=3)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
}

Simultation_without_replicates_Z<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set){
  #number_repitition_matrix<-c(10,500,1000)
  sample_number_matrix<-c(100,500,1000)
  Bias_matrix<-matrix(nrow=3,ncol=3)
  SD_Mean_matrix<-matrix(nrow=3,ncol=3)
  SD_SE_matrix<-matrix(nrow=3,ncol=3)
  SDE_matrix<-matrix(nrow=3,ncol=3)
  for(q in 1:3){
    sample_number<-sample_number_matrix[q]
    simexresult_matrix_1<-matrix(nrow=50,ncol=1)
    simexresult_matrix_0<-matrix(nrow=50,ncol=1)
    simexresult_matrix_Z<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    SD_Z_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      XI<-seq(0,2,0.4)
      B=500
      #i=100
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      X <- matrix(rep(x_real1,3),nrow=length(x_real1),3)
      n <- dim(X)[1]
      rp <- dim(X)[2]
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
      Z<-rnorm(sample_number,0,1)
      y1 <- 1+0.8*x_real1 +0.6*Z+ rnorm(sample_number, sd = 0.04)
      est<-estimate_triple(W)
      sig_mul_est<-est$sig_mul
      sig_add_est<-est$sig_add
      #W_new<-rowMeans(W)
      #W_new<-rowMedians(W)
      W_new<-W[,1]
      simexresult<-simex_Z(B,sig_mul_est,sig_add_est,W_new,XI,y1,Z)
      simexresult_matrix_1[index,1]<-simexresult$simex_beta_1
      simexresult_matrix_0[index,1]<-simexresult$simex_beta_0
      simexresult_matrix_Z[index,1]<-simexresult$simex_beta_Z
      SD_0_matrix[index,1]<-simexresult$SD_0
      SD_1_matrix[index,1]<-simexresult$SD_1
      SD_Z_matrix[index,1]<-simexresult$SD_Z
    }
    
    Bias_beta_1<-mean(simexresult_matrix_1)-0.8
    SDE_1<-sd(simexresult_matrix_1,na.rm = TRUE)
    SE_1_mean<-mean(SD_1_matrix,na.rm = TRUE)
    SE_1_sd<-sd(SD_1_matrix,na.rm = TRUE)
    
    Bias_beta_0<-mean(simexresult_matrix_0,na.rm = TRUE)-1
    SDE_0<-sd(simexresult_matrix_0,na.rm = TRUE)
    SE_0_mean<-mean(SD_0_matrix,na.rm = TRUE)
    SE_0_sd<-sd(SD_0_matrix,na.rm = TRUE)
    
    Bias_beta_Z<-mean(simexresult_matrix_Z,na.rm = TRUE)-0.6
    SDE_Z<-sd(simexresult_matrix_Z,na.rm = TRUE)
    SE_Z_mean<-mean(SD_Z_matrix,na.rm = TRUE)
    SE_Z_sd<-sd(SD_Z_matrix,na.rm = TRUE)
    
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    Bias_matrix[q,3]<-Bias_beta_Z
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SDE_matrix[q,3]<-SDE_Z
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_Mean_matrix[q,3]<-SE_Z_mean
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    SD_SE_matrix[q,3]<-SE_Z_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=3)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
}

simex_with_replicates<-function(B,sig_mul_vec,sig_add_vec,W_rep,XI,y1,numrep){
  i=dim(W_rep)
  i=i[1]
  #k is the number of replicates 
  Theta_0_Full<-matrix(nrow=length(XI),ncol=numrep)
  Theta_1_Full<-matrix(nrow=length(XI),ncol=numrep)
  Var_0_Full<-matrix(nrow=length(XI),ncol=numrep)
  Var_1_Full<-matrix(nrow=length(XI),ncol=numrep)
  #This for loop creat the remeasured data for k replicates 
for (nrep in 1:numrep) {
  W=W_rep[,nrep]
  #B<-500
  Ub<-matrix(nrow=B,ncol=i)
  for(k in 1:B){
    Ub[k,]<-rnorm(i,0,sig_mul_vec)}
  
  Vb<-matrix(nrow=B,ncol=i)
  for(k in 1:B){
    Vb[k,]<-rnorm(i,0,sig_add_vec)}
  
  l<-length(XI)
  Wbi<-array(0,dim=c(B,i,l))
  m=200
  n=20
  for(k in 1:l){
    for(m in 1:B){
      for(n in 1:i){
        Wbi[m,n,k]<-W[n]*(exp(sqrt(XI[k])*Ub[m,n]))+sqrt(XI[k])*Vb[m,n]*exp(sqrt(XI[k])*Ub[m,n])
        #  Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]
        #   Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]*exp(log(sqrt(XI[k])*Ub[m,n]))
      }
    }
  }
  
  # compute the coefficient of the beta
  Theta_1<-matrix(nrow=B,ncol=l)
  Theta_0<-matrix(nrow=B,ncol=l)
  # Theta_Z<-matrix(nrow=B,ncol=l)
  var_est_1<-matrix(nrow=B,ncol=l)
  var_est_0<-matrix(nrow=B,ncol=l)
  # var_est_Z<-matrix(nrow=B,ncol=l)
  delta_Theta_0<-matrix(nrow=B,ncol=l)
  delta_Theta_1<-matrix(nrow=B,ncol=l)
  delta_Theta_Z<-matrix(nrow=B,ncol=l)
  
  for (l in 1:l){
    for (k in 1:B){
      lm_result<-summary(lm(y1 ~ Wbi[k,,l]))
      Theta_1[k,l]<-lm_result$coefficients[2,1]
      var_est_1[k,l]<-(lm_result$coefficients[2,2])^2
      Theta_0[k,l]<-lm_result$coefficients[1,1]
      var_est_0[k,l]<-(lm_result$coefficients[1,2])^2
      #   Theta_Z<-lm_result$coefficients[3,1]
    }
  }
  
  Theta_1_mean<- colMeans(Theta_1,na.rm=T)
  var_est_1_mean<-colMeans(var_est_1,na.rm=T)
  Theta_0_mean<- colMeans(Theta_0,na.rm=T)
  var_est_0_mean<-colMeans(var_est_0,na.rm=T)
  
  Theta_1_Full[,nrep ]<-Theta_1_mean
  Theta_0_Full[,nrep ]<-Theta_0_mean
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_1[k,l]<-(Theta_1[k,l]-Theta_1_mean[l])^2
    }
  }
  
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_0[k,l]<-(Theta_0[k,l]-Theta_0_mean[l])^2
    }
  }
  
  #s^2
  delta_Theta_1_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_1_sqr[1,l]<-sum(delta_Theta_1[,l])/(B-1)
    }
  }
  
  #s^2 for the first coefficient
  delta_Theta_0_sqr<-matrix(nrow=1,ncol=l)
  for (l in 1:l){
    for (k in 1:B){
      delta_Theta_0_sqr[1,l]<-sum(delta_Theta_0[,l])/(B-1)
    }
     } 
    diff_1=var_est_1_mean-delta_Theta_1_sqr
    diff_0=var_est_0_mean-delta_Theta_1_sqr
    Var_1_Full[,nrep]= t(diff_1)
    Var_0_Full[,nrep]= t(diff_0)
  }
  
  #compute the difference
  #diff_1=var_est_1_mean-delta_Theta_1_sqr
  simex_var_1<-Exploration(as.vector(Var_1_Full),XI,numrep)$beta_simex
  simex_beta_1<-Exploration(as.vector(Theta_1_Full),XI,numrep)$beta_simex
  
  #compute the difference
  diff_0=var_est_0_mean-delta_Theta_1_sqr
  simex_var_0<-Exploration(as.vector(Var_0_Full),XI,numrep)$beta_simex
  simex_beta_0<-Exploration(as.vector(Theta_0_Full),XI,numrep)$beta_simex
  
  #t value
  t_1<-simex_beta_1/sqrt(simex_var_0)
  
  #t value
  t_0<-simex_beta_0/sqrt(simex_var_1)
  
  #p value
  p_1<- 2*pt(q=t_1, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_0,Theta_1_mean)
  
  #p value
  p_0<- 2*pt(q=t_0, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_1,Theta_1_mean)
  # return<-list(simex_beta=simex_beta,ylab=ylab,beta0=beta0,beta1=beta1,beta2=beta2,tvalue=t,pvalue=p)
  # return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,tvalue=t,pvalue=p,Theta_1_mean=Theta_1_mean,diff=diff,simex_var,stderr=sqrt(simex_var))
  return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,Theta_1_mean=Theta_1_mean,simex_beta_1=simex_beta_1,SD_0=sqrt(simex_var_0),SD_1=sqrt(simex_var_1),Theta_0_Full=Theta_0_Full,Theta_1_Full=Theta_1_Full,Var_0_Full=Var_0_Full,Var_1_Full=Var_1_Full)
}


simex_with_replicates_Z<-function(B,sig_mul_vec,sig_add_vec,W_rep,XI,y1,numrep,Z){
  i=dim(W_rep)
  i=i[1]
  #k is the number of replicates 
  Theta_0_Full<-matrix(nrow=length(XI),ncol=numrep)
  Theta_1_Full<-matrix(nrow=length(XI),ncol=numrep)
  Var_0_Full<-matrix(nrow=length(XI),ncol=numrep)
  Var_1_Full<-matrix(nrow=length(XI),ncol=numrep)
  Theta_Z_Full<-matrix(nrow=length(XI),ncol=numrep)
  Var_Z_Full<-matrix(nrow=length(XI),ncol=numrep)
  #This for loop creat the remeasured data for k replicates 
  for (nrep in 1:numrep) {
    W=W_rep[,nrep]
    #B<-500
    Ub<-matrix(nrow=B,ncol=i)
    for(k in 1:B){
      Ub[k,]<-rnorm(i,0,sig_mul_vec)}
    
    Vb<-matrix(nrow=B,ncol=i)
    for(k in 1:B){
      Vb[k,]<-rnorm(i,0,sig_add_vec)}
    
    l<-length(XI)
    Wbi<-array(0,dim=c(B,i,l))
    m=200
    n=20
    for(k in 1:l){
      for(m in 1:B){
        for(n in 1:i){
          Wbi[m,n,k]<-W[n]*(exp(sqrt(XI[k])*Ub[m,n]))+sqrt(XI[k])*Vb[m,n]*exp(sqrt(XI[k])*Ub[m,n])
          #  Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]
          #   Wbi[m,n,k]<-(exp(log(W[n])+(sqrt(XI[k])*Ub[m,n])))+sqrt(XI[k])*Vb[m,n]*exp(log(sqrt(XI[k])*Ub[m,n]))
        }
      }
    }
    
    # compute the coefficient of the beta
    Theta_1<-matrix(nrow=B,ncol=l)
    Theta_0<-matrix(nrow=B,ncol=l)
    Theta_Z<-matrix(nrow=B,ncol=l)
    # Theta_Z<-matrix(nrow=B,ncol=l)
    var_est_1<-matrix(nrow=B,ncol=l)
    var_est_0<-matrix(nrow=B,ncol=l)
    var_est_Z<-matrix(nrow=B,ncol=l)
    # var_est_Z<-matrix(nrow=B,ncol=l)
    delta_Theta_0<-matrix(nrow=B,ncol=l)
    delta_Theta_1<-matrix(nrow=B,ncol=l)
    delta_Theta_Z<-matrix(nrow=B,ncol=l)
    
    for (l in 1:l){
      for (k in 1:B){
        lm_result<-summary(lm(y1 ~ Wbi[k,,l]+Z))
        Theta_1[k,l]<-lm_result$coefficients[2,1]
        var_est_1[k,l]<-(lm_result$coefficients[2,2])^2
        Theta_0[k,l]<-lm_result$coefficients[1,1]
        var_est_0[k,l]<-(lm_result$coefficients[1,2])^2
        Theta_Z[k,l]<-lm_result$coefficients[3,1]
        var_est_Z[k,l]<-(lm_result$coefficients[3,2])^2
      }
    }
    
    Theta_1_mean<- colMeans(Theta_1,na.rm=T)
    var_est_1_mean<-colMeans(var_est_1,na.rm=T)
    Theta_0_mean<- colMeans(Theta_0,na.rm=T)
    var_est_0_mean<-colMeans(var_est_0,na.rm=T)
    Theta_Z_mean<- colMeans(Theta_Z,na.rm=T)
    var_est_Z_mean<-colMeans(var_est_Z,na.rm=T)
    
    Theta_1_Full[,nrep ]<-Theta_1_mean
    Theta_0_Full[,nrep ]<-Theta_0_mean
    Theta_Z_Full[,nrep ]<-Theta_Z_mean
    for (l in 1:l){
      for (k in 1:B){
        delta_Theta_1[k,l]<-(Theta_1[k,l]-Theta_1_mean[l])^2
      }
    }
    
    for (l in 1:l){
      for (k in 1:B){
        delta_Theta_0[k,l]<-(Theta_0[k,l]-Theta_0_mean[l])^2
      }
    }
    
    for (l in 1:l){
      for (k in 1:B){
        delta_Theta_Z[k,l]<-(Theta_Z[k,l]-Theta_Z_mean[l])^2
      }
    }
    
    #s^2
    delta_Theta_1_sqr<-matrix(nrow=1,ncol=l)
    for (l in 1:l){
      for (k in 1:B){
        delta_Theta_1_sqr[1,l]<-sum(delta_Theta_1[,l])/(B-1)
      }
    }
    
    #s^2 for the first coefficient
    delta_Theta_0_sqr<-matrix(nrow=1,ncol=l)
    for (l in 1:l){
      for (k in 1:B){
        delta_Theta_0_sqr[1,l]<-sum(delta_Theta_0[,l])/(B-1)
      }
    } 
    
    delta_Theta_Z_sqr<-matrix(nrow=1,ncol=l)
    for (l in 1:l){
      for (k in 1:B){
        delta_Theta_Z_sqr[1,l]<-sum(delta_Theta_Z[,l])/(B-1)
      }
    }
    diff_1=var_est_1_mean-delta_Theta_1_sqr
    diff_0=var_est_0_mean-delta_Theta_0_sqr
    diff_Z=var_est_Z_mean-delta_Theta_Z_sqr
    Var_1_Full[,nrep]= t(diff_1)
    Var_0_Full[,nrep]= t(diff_0)
    Var_Z_Full[,nrep]= t(diff_Z)
  }
  
  #compute the difference
  #diff_1=var_est_1_mean-delta_Theta_1_sqr
  simex_var_1<-Exploration(as.vector(Var_1_Full),XI,numrep)$beta_simex
  simex_beta_1<-Exploration(as.vector(Theta_1_Full),XI,numrep)$beta_simex
  
  #compute the difference
 # diff_0=var_est_0_mean-delta_Theta_1_sqr
  simex_var_0<-Exploration(as.vector(Var_0_Full),XI,numrep)$beta_simex
  simex_beta_0<-Exploration(as.vector(Theta_0_Full),XI,numrep)$beta_simex
  
  #compute the difference
 # diff_0=var_est_0_mean-delta_Theta_1_sqr
  simex_var_Z<-Exploration(as.vector(Var_Z_Full),XI,numrep)$beta_simex
  simex_beta_Z<-Exploration(as.vector(Theta_Z_Full),XI,numrep)$beta_simex
  
  #t value
  t_1<-simex_beta_1/sqrt(simex_var_0)
  
  #t value
  t_0<-simex_beta_0/sqrt(simex_var_0)
  
  #t value
  t_Z<-simex_beta_Z/sqrt(simex_var_Z)
  
  #p value
  p_1<- 2*pt(q=t_1, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_0,Theta_1_mean)
  
  #p value
  p_0<- 2*pt(q=t_0, df=i-2, lower.tail=FALSE)
  ylab<-c(simex_beta_1,Theta_1_mean)
  # return<-list(simex_beta=simex_beta,ylab=ylab,beta0=beta0,beta1=beta1,beta2=beta2,tvalue=t,pvalue=p)
  # return<-list(simex_beta_1=simex_beta_1,simex_beta_0=simex_beta_0,tvalue=t,pvalue=p,Theta_1_mean=Theta_1_mean,diff=diff,simex_var,stderr=sqrt(simex_var))
  return<-list(simex_beta_1=simex_beta_1,simex_beta_Z=simex_beta_Z,simex_beta_0=simex_beta_0,Theta_1_mean=Theta_1_mean,simex_beta_1=simex_beta_1,SD_0=sqrt(simex_var_0),SD_1=sqrt(simex_var_1),SD_Z=sqrt(simex_var_Z))
}

Exploration<- function(data,XI,numrep){
 # EXmodel<-lm(as.vector(data)~XI+I(XI^2))
  Response<-c(data)
  full_XI<-c(matrix(rep(XI, numrep), nrow = 1))
  EXmodel<-lm((Response)~full_XI+I(full_XI^2))
  beta0<-summary(EXmodel)$coefficients[1,1]
  beta1<-summary(EXmodel)$coefficients[2,1]
  beta2<-summary(EXmodel)$coefficients[3,1]
  beta_simex<-beta0+beta1*(-1)+beta2*1
  return<-list(beta_simex=beta_simex,beta0=beta0,beta1=beta1,beta2=beta2) 
  # return(beta_simex)
}

simex_plot<-function(XI,data,beta_simex,numrep){
  #XI_new<-c(-1, XI)
 data=test$Theta_1_Full
  beta_simex=test$simex_beta_1
  #full_XI<-c(matrix(rep(XI, numrep), nrow = 1))
  Response<-c(data)
  full_XI<-c(matrix(rep(XI, numrep), nrow = 1))
  XI_new<-c(-1, full_XI)
  ylab<-c(beta_simex,Response)
  beta0<-Exploration(as.vector(data),XI,numrep)$beta0
  beta1<-Exploration(as.vector(data),XI,numrep)$beta1
  beta2<-Exploration(as.vector(data),XI,numrep)$beta2
  Exploration_function <- function(x) beta0+beta1*(x)+beta2*x^2
  data<-data.frame(XI_new,ylab)
  df1<-data.frame(x=seq(-1,0,length.out = 100),
                  y=Exploration_function(seq(-1,0,length.out = 100)))
  df2<-data.frame(x=seq(0,2,length.out = 100),
                  y=Exploration_function(seq(0,2,length.out = 100)))
  p<-ggplot() + geom_point(data=data,aes(x=XI_new,y=ylab))+
    
    # first line
    geom_line(data = df1, aes(x = x, y = y), color = "blue", size = 1.5) +
    
    # second line
    geom_line(data = df2, aes(x = x, y = y), color = "red", size = 1.5, linetype = "dashed")+
  scale_y_continuous(name='')+scale_x_continuous(name=TeX("$\\xi$"))
  return(p)
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
  if(a2<1){a2_est <- 1}
  if(a2<1){a2 <- 1}
  
  
  # estimate of sig_mul
  sig_mul <- sqrt(2*log(sqrt(a2)))
  
  # estimate of sig_add_2^2
  sig_add_2 <- mean(W2)-a2*mean(Wij)
  if(sig_add_2<0) {sig_add_2<-0}
  # estimate of sig_add
  sig_add <- sqrt(sig_add_2)
  
  result <- list(a2 = a2,sig_mul = sig_mul,sig_add_2 = sig_add_2,sig_add = sig_add, EWW=mean(Wij), EW=mean(W), W2=W2, W2W=W2W,n=n,Wij=Wij,W=W,Wijk=Wijk,Wijk2=Wijk2,W2=W2,W3=W3,W4=W4, W3W=W3W,
                 W3W2=W3W2,W2W2=W2W2,W4W=W4W,Wi2jk=Wi2jk,Wi3jk=Wi3jk,Wi2j2k=Wi2j2k,Wi3j2k=Wi3j2k)
  return(result)
}

Change_Distribution<- function(distribution,sample_number){
  if (distribution ==1 ){
    x_real1<-rlnorm(sample_number,0,0.5)
  } else if ( distribution == 2 ) {
    x_real1<-rexp(sample_number,0.5)
  } else {
    x_real1 <- rchisq(sample_number,2)
  }
  return<- x_real1
}

#simulation method for simex without Z but with estimated value 
Simultation_with_estimate<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set,numrep){
  #number_repitition_matrix<-c(10,500,1000)
  sample_number_matrix<-c(100,500,2000)
  Bias_matrix<-matrix(nrow=3,ncol=2)
  SD_Mean_matrix<-matrix(nrow=3,ncol=2)
  SD_SE_matrix<-matrix(nrow=3,ncol=2)
  SDE_matrix<-matrix(nrow=3,ncol=2)
  #numrep1=numrep
  for(q in 1:length(sample_number_matrix)){
    sample_number<-sample_number_matrix[q]
    simexresult_matrix_1<-matrix(nrow=50,ncol=1)
    simexresult_matrix_0<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      XI<-seq(0,2,0.4)
      B=50
      #i=100
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      X <- matrix(rep(x_real1,numrep),nrow=length(x_real1),numrep)
      n <- dim(X)[1]
      rp <- dim(X)[2]
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W_rep <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
      y1 <- 1+0.8*x_real1 + rnorm(sample_number, sd = 0.04)
      est<-estimate_triple(W_rep)
      sig_mul_est<-est$sig_mul
      sig_add_est<-est$sig_add
      #W_new<-rowMeans(W)
      #W_new<-rowMedians(W)
      #W_new<-W[,1]
      simexresult<-simex_with_replicates(B,sig_mul_est, sig_add_est,W_rep,XI,y1,numrep)
      simexresult_matrix_1[index,1]<-simexresult$simex_beta_1
      simexresult_matrix_0[index,1]<-simexresult$simex_beta_0
      SD_0_matrix[index,1]<-simexresult$SD_0
      SD_1_matrix[index,1]<-simexresult$SD_1
    }
    
    Bias_beta_1<-mean(simexresult_matrix_1)-0.8
    SDE_1<-sd(simexresult_matrix_1,na.rm = TRUE)
    SE_1_mean<-mean(SD_1_matrix,na.rm = TRUE)
    SE_1_sd<-sd(SD_1_matrix,na.rm = TRUE)
    Bias_beta_0<-mean(simexresult_matrix_0,na.rm = TRUE)-1
    SDE_0<-sd(simexresult_matrix_0,na.rm = TRUE)
    SE_0_mean<-mean(SD_0_matrix,na.rm = TRUE)
    SE_0_sd<-sd(SD_0_matrix,na.rm = TRUE)
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=2)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
 # result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix)
}

#simulation method for simex without Z but with true variance
Simultation_replicate_trueVariance<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set,numrep){
  #number_repitition_matrix<-c(10,500,1000)
  sample_number_matrix<-c(100,500,2000)
  Bias_matrix<-matrix(nrow=3,ncol=2)
  SD_Mean_matrix<-matrix(nrow=3,ncol=2)
  SD_SE_matrix<-matrix(nrow=3,ncol=2)
  SDE_matrix<-matrix(nrow=3,ncol=2)
  #numrep1=numrep
  for(q in 1:length(sample_number_matrix)){
    sample_number<-sample_number_matrix[q]
    simexresult_matrix_1<-matrix(nrow=50,ncol=1)
    simexresult_matrix_0<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      XI<-seq(0,2,0.4)
      B=50
      #i=100
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      X <- matrix(rep(x_real1,numrep),nrow=length(x_real1),numrep)
      n <- dim(X)[1]
      rp <- dim(X)[2]
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W_rep <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
      y1 <- 1+0.8*x_real1 + rnorm(sample_number, sd = 0.04)
      est<-estimate_triple(W_rep)
      #sig_mul_est<-est$sig_mul
      #sig_add_est<-est$sig_add
      #W_new<-rowMeans(W)
      #W_new<-rowMedians(W)
      #W_new<-W[,1]
      simexresult<-simex_with_replicates(B,sig_mul_vec, sig_add_vec,W_rep,XI,y1,numrep)
      simexresult_matrix_1[index,1]<-simexresult$simex_beta_1
      simexresult_matrix_0[index,1]<-simexresult$simex_beta_0
      SD_0_matrix[index,1]<-simexresult$SD_0
      SD_1_matrix[index,1]<-simexresult$SD_1
    }
    
    Bias_beta_1<-mean(simexresult_matrix_1)-0.8
    SDE_1<-sd(simexresult_matrix_1,na.rm = TRUE)
    SE_1_mean<-mean(SD_1_matrix,na.rm = TRUE)
    SE_1_sd<-sd(SD_1_matrix,na.rm = TRUE)
    Bias_beta_0<-mean(simexresult_matrix_0,na.rm = TRUE)-1
    SDE_0<-sd(simexresult_matrix_0,na.rm = TRUE)
    SE_0_mean<-mean(SD_0_matrix,na.rm = TRUE)
    SE_0_sd<-sd(SD_0_matrix,na.rm = TRUE)
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=2)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
  # result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix)
}

Naive_result<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set){
  sample_number_matrix<-c(100,500,2000)
  Bias_matrix<-matrix(nrow=3,ncol=2)
  SD_Mean_matrix<-matrix(nrow=3,ncol=2)
  SD_SE_matrix<-matrix(nrow=3,ncol=2)
  SDE_matrix<-matrix(nrow=3,ncol=2)
  for(q in 1:3){
    sample_number<-sample_number_matrix[q]
    naive_matrix_1<-matrix(nrow=50,ncol=1)
    naive_matrix_0<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W <- x_real1*U+V
      y1 <- 1+0.8*x_real1 + rnorm(sample_number, sd = 0.04)
      naive_lm_result<-summary(lm(y1~W))
      naive_matrix_1[index,1]<-naive_lm_result$coefficients[2,1]
      naive_matrix_0[index,1]<-naive_lm_result$coefficients[1,1]
      SD_0_matrix[index,1]<-naive_lm_result$coefficients[1,2]
      SD_1_matrix[index,1]<-naive_lm_result$coefficients[2,2]
    }
    
    Bias_beta_1<-mean(naive_matrix_1)-0.8
    SDE_1<-sd(naive_matrix_1)
    SE_1_mean<-mean(SD_1_matrix)
    SE_1_sd<-sd(SD_1_matrix)
    Bias_beta_0<-mean(naive_matrix_0)-1
    SDE_0<-sd(naive_matrix_0)
    SE_0_mean<-mean(SD_0_matrix)
    SE_0_sd<-sd(SD_0_matrix)
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=2)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
  #result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix)
}

Naive_result_Z<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set){
  sample_number_matrix<-c(100,500,2000)
  Bias_matrix<-matrix(nrow=3,ncol=3)
  SD_Mean_matrix<-matrix(nrow=3,ncol=3)
  SD_SE_matrix<-matrix(nrow=3,ncol=3)
  SDE_matrix<-matrix(nrow=3,ncol=3)
  for(q in 1:3){
    sample_number<-sample_number_matrix[q]
    naive_matrix_1<-matrix(nrow=50,ncol=1)
    naive_matrix_0<-matrix(nrow=50,ncol=1)
    naive_matrix_Z<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    SD_Z_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      Z=rnorm(sample_number,0,1)
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W <- x_real1*U+V
      
      y1 <- 1 + 0.8*x_real1 + 0.6*Z + rnorm(sample_number, sd = 0.04)
      naive_lm_result<-summary(lm(y1~W+Z))
      naive_matrix_1[index,1]<-naive_lm_result$coefficients[2,1]
      naive_matrix_0[index,1]<-naive_lm_result$coefficients[1,1]
      naive_matrix_Z[index,1]<-naive_lm_result$coefficients[3,1]
      SD_0_matrix[index,1]<-naive_lm_result$coefficients[1,2]
      SD_1_matrix[index,1]<-naive_lm_result$coefficients[2,2]
      SD_Z_matrix[index,1]<-naive_lm_result$coefficients[3,2]
    }
    
    Bias_beta_1<-mean(naive_matrix_1)-0.8
    SDE_1<-sd(naive_matrix_1)
    SE_1_mean<-mean(SD_1_matrix)
    SE_1_sd<-sd(SD_1_matrix)
    
    Bias_beta_0<-mean(naive_matrix_0)-1
    SDE_0<-sd(naive_matrix_0)
    SE_0_mean<-mean(SD_0_matrix)
    SE_0_sd<-sd(SD_0_matrix)
    
    Bias_beta_Z<-mean(naive_matrix_Z)-0.6
    SDE_Z<-sd(naive_matrix_Z)
    SE_Z_mean<-mean(SD_Z_matrix)
    SE_Z_sd<-sd(SD_Z_matrix)
    
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    Bias_matrix[q,3]<-Bias_beta_Z
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SDE_matrix[q,3]<-SDE_Z
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_Mean_matrix[q,3]<-SE_Z_mean
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    SD_SE_matrix[q,3]<-SE_Z_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=3)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
  #result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix)
}


Simultation_replicate_trueVariance_Z<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set,numrep){
  #number_repitition_matrix<-c(10,500,1000)
  sample_number_matrix<-c(100,500,2000)
  #sample_number_matrix<-c(50,50,50)
  Bias_matrix<-matrix(nrow=3,ncol=3)
  SD_Mean_matrix<-matrix(nrow=3,ncol=3)
  SD_SE_matrix<-matrix(nrow=3,ncol=3)
  SDE_matrix<-matrix(nrow=3,ncol=3)
  #numrep1=numrep
  for(q in 1:length(sample_number_matrix)){
    sample_number<-sample_number_matrix[q]
    simexresult_matrix_1<-matrix(nrow=50,ncol=1)
    simexresult_matrix_0<-matrix(nrow=50,ncol=1)
    simexresult_matrix_Z<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    SD_Z_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      XI<-seq(0,2,0.4)
      B=50
      #i=100
      Z<-rnorm(sample_number,0,1)
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      X <- matrix(rep(x_real1,numrep),nrow=length(x_real1),numrep)
      n <- dim(X)[1]
      rp <- dim(X)[2]
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W_rep <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
      y1 <- 1+0.8*x_real1 +0.6*Z+ rnorm(sample_number, sd = 0.04)
     # est<-estimate_triple(W_rep)
      #sig_mul_est<-est$sig_mul
      #sig_add_est<-est$sig_add
      #W_new<-rowMeans(W)
      #W_new<-rowMedians(W)
      #W_new<-W[,1]
      simexresult<-simex_with_replicates_Z(B,sig_mul_vec, sig_add_vec,W_rep,XI,y1,numrep,Z)
      simexresult_matrix_1[index,1]<-simexresult$simex_beta_1
      simexresult_matrix_0[index,1]<-simexresult$simex_beta_0
      simexresult_matrix_Z[index,1]<-simexresult$simex_beta_Z
      SD_0_matrix[index,1]<-simexresult$SD_0
      SD_1_matrix[index,1]<-simexresult$SD_1
      SD_Z_matrix[index,1]<-simexresult$SD_Z
    }
    
    Bias_beta_1<-mean(simexresult_matrix_1)-0.8
    SDE_1<-sd(simexresult_matrix_1,na.rm = TRUE)
    SE_1_mean<-mean(SD_1_matrix,na.rm = TRUE)
    SE_1_sd<-sd(SD_1_matrix,na.rm = TRUE)
    Bias_beta_0<-mean(simexresult_matrix_0,na.rm = TRUE)-1
    SDE_0<-sd(simexresult_matrix_0,na.rm = TRUE)
    SE_0_mean<-mean(SD_0_matrix,na.rm = TRUE)
    SE_0_sd<-sd(SD_0_matrix,na.rm = TRUE)
    
    Bias_beta_Z<-mean(simexresult_matrix_Z,na.rm = TRUE)-0.6
    SDE_Z<-sd(simexresult_matrix_Z,na.rm = TRUE)
    SE_Z_mean<-mean(SD_Z_matrix,na.rm = TRUE)
    SE_Z_sd<-sd(SD_Z_matrix,na.rm = TRUE)
    
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    Bias_matrix[q,3]<-Bias_beta_Z
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SDE_matrix[q,3]<-SDE_Z
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_Mean_matrix[q,3]<-SE_Z_mean
  
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    SD_SE_matrix[q,3]<-SE_Z_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=3)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
  # result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix)
}

Simultation_replicate_estimate_Z<-function(Distri_X,sig_mul_vec_set,sig_add_vec_set,numrep){
  #number_repitition_matrix<-c(10,500,1000)
  sample_number_matrix<-c(100,500,2000)
  #sample_number_matrix<-c(50,50,50)
  Bias_matrix<-matrix(nrow=3,ncol=3)
  SD_Mean_matrix<-matrix(nrow=3,ncol=3)
  SD_SE_matrix<-matrix(nrow=3,ncol=3)
  SDE_matrix<-matrix(nrow=3,ncol=3)
  #numrep1=numrep
  for(q in 1:length(sample_number_matrix)){
    sample_number<-sample_number_matrix[q]
    simexresult_matrix_1<-matrix(nrow=50,ncol=1)
    simexresult_matrix_0<-matrix(nrow=50,ncol=1)
    simexresult_matrix_Z<-matrix(nrow=50,ncol=1)
    SD_0_matrix<-matrix(nrow=50,ncol=1)
    SD_1_matrix<-matrix(nrow=50,ncol=1)
    SD_Z_matrix<-matrix(nrow=50,ncol=1)
    
    # each data set contain 100 observation
    
    for (index in 1:50){
      XI<-seq(0,2,0.4)
      B=50
      #i=100
      Z<-rnorm(sample_number,0,1)
      x_real1 <- Change_Distribution(Distri_X,sample_number)
      X <- matrix(rep(x_real1,numrep),nrow=length(x_real1),numrep)
      n <- dim(X)[1]
      rp <- dim(X)[2]
      sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
      sig_add_vec<-sig_add_vec_set*sd(x_real1)
      U<-rlnorm(sample_number,0,sig_mul_vec)
      V<-rnorm(sample_number,0,sig_add_vec)
      W_rep <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)
      y1 <- 1+0.8*x_real1 +0.6*Z+ rnorm(sample_number, sd = 0.04)
      est<-estimate_triple(W_rep)
      sig_mul_est<-est$sig_mul
      sig_add_est<-est$sig_add
      simexresult<-simex_with_replicates_Z(B,est$sig_mul, est$sig_add,W_rep,XI,y1,numrep,Z)
      simexresult_matrix_1[index,1]<-simexresult$simex_beta_1
      simexresult_matrix_0[index,1]<-simexresult$simex_beta_0
      simexresult_matrix_Z[index,1]<-simexresult$simex_beta_Z
      SD_0_matrix[index,1]<-simexresult$SD_0
      SD_1_matrix[index,1]<-simexresult$SD_1
      SD_Z_matrix[index,1]<-simexresult$SD_Z
    }
    
    Bias_beta_1<-mean(simexresult_matrix_1)-0.8
    SDE_1<-sd(simexresult_matrix_1,na.rm = TRUE)
    SE_1_mean<-mean(SD_1_matrix,na.rm = TRUE)
    SE_1_sd<-sd(SD_1_matrix,na.rm = TRUE)
    Bias_beta_0<-mean(simexresult_matrix_0,na.rm = TRUE)-1
    SDE_0<-sd(simexresult_matrix_0,na.rm = TRUE)
    SE_0_mean<-mean(SD_0_matrix,na.rm = TRUE)
    SE_0_sd<-sd(SD_0_matrix,na.rm = TRUE)
    
    Bias_beta_Z<-mean(simexresult_matrix_Z,na.rm = TRUE)-0.6
    SDE_Z<-sd(simexresult_matrix_Z,na.rm = TRUE)
    SE_Z_mean<-mean(SD_Z_matrix,na.rm = TRUE)
    SE_Z_sd<-sd(SD_Z_matrix,na.rm = TRUE)
    
    Bias_matrix[q,1]<-Bias_beta_0
    Bias_matrix[q,2]<-Bias_beta_1
    Bias_matrix[q,3]<-Bias_beta_Z
    SDE_matrix[q,1]<-SDE_0
    SDE_matrix[q,2]<-SDE_1
    SDE_matrix[q,3]<-SDE_Z
    SD_Mean_matrix[q,1]<-SE_0_mean
    SD_Mean_matrix[q,2]<-SE_1_mean
    SD_Mean_matrix[q,3]<-SE_Z_mean
    
    SD_SE_matrix[q,1]<-SE_0_sd
    SD_SE_matrix[q,2]<-SE_1_sd
    SD_SE_matrix[q,3]<-SE_Z_sd
    Bias_matrix<-round(Bias_matrix,3)
    SDE_matrix<-round(SDE_matrix,3)
    SD_Mean_matrix<-round(SD_Mean_matrix,3)
    SD_SE_matrix<-round(SD_SE_matrix,3)
  }
  Full_table<-matrix(nrow=12,ncol=3)
  Full_table[1,]<-Bias_matrix[1,]
  Full_table[2,]<-SDE_matrix[1,]
  Full_table[3,]<-SD_Mean_matrix[1,]
  Full_table[4,]<-SD_SE_matrix[1,]
  Full_table[5,]<-Bias_matrix[2,]
  Full_table[6,]<-SDE_matrix[2,]
  Full_table[7,]<-SD_Mean_matrix[2,]
  Full_table[8,]<-SD_SE_matrix[1,]
  Full_table[9,]<-Bias_matrix[3,]
  Full_table[10,]<-SDE_matrix[3,]
  Full_table[11,]<-SD_Mean_matrix[3,]
  Full_table[12,]<-SD_SE_matrix[3,]
  result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix,Full_table=Full_table)
  # result<-list(Bias_matrix=Bias_matrix,SD_Mean_matrix=SD_Mean_matrix,SD_SE_matrix=SD_SE_matrix,SDE_matrix=SDE_matrix)
}

lab <- data.frame(label=rep(c("Bias","S.D.","Mean($\\widehat{S.E.}$)","S.D.($\\widehat{S.E.}$)"),9))
size <- data.frame(size=rep(c("n=100","","","","n=500","","","","n=2000","","",""),3))
dist <- data.frame(dist=rep(c("LogNorm(0,0.5)",rep("",11),"Exp(0.5)",rep("",11),"Chisq(3)",rep("",11))))



sample_number=500
y1 <- 1+0.8*x_real1 +0.6*Z+ rnorm(sample_number, sd = 0.04)
x_real1 <- Change_Distribution(1,sample_number)
X <- matrix(rep(x_real1,3),nrow=length(x_real1),3)
n <- dim(X)[1]
rp <- dim(X)[2]
sig_mul_vec_set=0.1
sig_add_vec_set=0.3
sig_mul_vec<-sig_mul_vec_set*sd(x_real1)
sig_add_vec<-sig_add_vec_set*sd(x_real1)
U<-rlnorm(sample_number,0,sig_mul_vec)
V<-rnorm(sample_number,0,sig_add_vec)
W_rep <- X*matrix(rlnorm(n*rp,0,sig_mul_vec),n,rp)+matrix(rnorm(n*rp,0,sig_add_vec),n,rp)





s1=Simultation_replicate_trueVariance_(1,0.3,0.3,3)
s2=Simultation_replicate_trueVariance_(1,0.1,0.3,3)
s3=Simultation_replicate_trueVariance_(2,0.3,0.3,3)
s4=Simultation_replicate_trueVariance_(2,0.1,0.3,3)
s5=Simultation_replicate_trueVariance_(3,0.3,0.3,3)
s6=Simultation_replicate_trueVariance_(3,0.1,0.3,3)


c1=Simultation_with_estimate(1,0.3,0.3,3)
c2=Simultation_with_estimate(1,0.1,0.3,3)
c3=Simultation_with_estimate(2,0.3,0.3,3)
c4=Simultation_with_estimate(2,0.1,0.3,3)
c5=Simultation_with_estimate(3,0.3,0.3,3)
c6=Simultation_with_estimate(3,0.1,0.3,3)



Naive_result_1<-Naive_result(1,0.3,0.3)
Naive_result_2<-Naive_result(1,0.1,0.3)
Naive_result_3<-Naive_result(2,0.3,0.3)
Naive_result_4<-Naive_result(2,0.1,0.3)
Naive_result_5<-Naive_result(3,0.3,0.3)
Naive_result_6<-Naive_result(3,0.1,0.3)

#table for 0.3 0.3
s1_table<-as.table(s1$Full_table)
s3_table<-as.table(s3$Full_table)
s5_table<-as.table(s5$Full_table)
ss2=rbind(s1_table,s3_table,s5_table)
ss

Naive_result_1_table<-as.table(Naive_result_1$Full_table)
Naive_result_3_table<-as.table(Naive_result_3$Full_table)
Naive_result_5_table<-as.table(Naive_result_5$Full_table)
ss1=rbind(Naive_result_1_table,Naive_result_3_table,Naive_result_5_table)

#table for 
c1_table<-as.table(c1$Full_table)
c3_table<-as.table(c3$Full_table)
c5_table<-as.table(c5$Full_table)
ss3=rbind(c1_table,c3_table,c5_table)

ss_Full=cbind(ss1,ss2,ss3)
colnames(ss_Full) <- c("beta0_ls", "beta1_ls", "beta0_real","beta1_real","beta_0","beta_1")
ss_Full

ss_Full_latex <- xtable(cbind(dist,size,lab,ss_Full),digits = 3)
print(ss_Full_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

#table for 0.1 0.3
s2_table<-as.table(s2$Full_table)
s4_table<-as.table(s4$Full_table)
s6_table<-as.table(s6$Full_table)
sss2=rbind(s2_table,s4_table,s6_table)

Naive_result_2_table<-as.table(Naive_result_2$Full_table)
Naive_result_4_table<-as.table(Naive_result_4$Full_table)
Naive_result_6_table<-as.table(Naive_result_6$Full_table)
sss1=rbind(Naive_result_2_table,Naive_result_4_table,Naive_result_6_table)

#table for 
c2_table<-as.table(c2$Full_table)
c4_table<-as.table(c4$Full_table)
c6_table<-as.table(c6$Full_table)
sss3=rbind(c2_table,c4_table,c6_table)

sss_Full=cbind(sss1,sss2,sss3)
colnames(sss_Full) <- c("beta0_ls", "beta1_ls", "beta0_real","beta1_real","beta_0","beta_1")
#sss_Full

sss_Full_latex <- xtable(cbind(dist,size,lab,sss_Full),digits = 3)
#o<-sss_Full_latex[,-1]
print(sss_Full_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)




######### simulation with true variance and Z
#s=Simultation_replicate_trueVariance_Z(2,0.3,0.3,3)
#s2222=Simultation_replicate_estimate_Z(Distri_X,sig_mul_vec_set,sig_add_vec_set,numrep)

replicate_with_true_Z_1<-Simultation_replicate_trueVariance_Z(1,0.3,0.3,3)
replicate_with_true_Z_2<-Simultation_replicate_trueVariance_Z(1,0.1,0.3,3)
replicate_with_true_Z_3<-Simultation_replicate_trueVariance_Z(2,0.3,0.3,3)
replicate_with_true_Z_4<-Simultation_replicate_trueVariance_Z(2,0.1,0.3,3)
replicate_with_true_Z_5<-Simultation_replicate_trueVariance_Z(3,0.3,0.3,3)
replicate_with_true_Z_6<-Simultation_replicate_trueVariance_Z(3,0.1,0.3,3)


############## simulation with estimate variance and Z 
replicate_with_estimate_Z_1<-Simultation_replicate_estimate_Z(1,0.3,0.3,3)
replicate_with_estimate_Z_2<-Simultation_replicate_estimate_Z(1,0.1,0.3,3)
replicate_with_estimate_Z_3<-Simultation_replicate_estimate_Z(2,0.3,0.3,3)
replicate_with_estimate_Z_4<-Simultation_replicate_estimate_Z(2,0.1,0.3,3)
replicate_with_estimate_Z_5<-Simultation_replicate_estimate_Z(3,0.3,0.3,3)
replicate_with_estimate_Z_6<-Simultation_replicate_estimate_Z(3,0.1,0.3,3)


######################
Naive_result_1_Z<-Naive_result_Z(1,0.3,0.3)
Naive_result_2_Z<-Naive_result_Z(1,0.1,0.3)
Naive_result_3_Z<-Naive_result_Z(2,0.3,0.3)
Naive_result_4_Z<-Naive_result_Z(2,0.1,0.3)
Naive_result_5_Z<-Naive_result_Z(3,0.3,0.3)
Naive_result_6_Z<-Naive_result_Z(3,0.1,0.3)

############### table for 0.3 0.3

#table for 0.3 0.3
Z1_table<-as.table(replicate_with_true_Z_1$Full_table)
Z3_table<-as.table(replicate_with_true_Z_3$Full_table)
Z5_table<-as.table(replicate_with_true_Z_5$Full_table)
ZZ2=rbind(Z1_table,Z3_table,Z5_table)
ss

Naive_result_1_table<-as.table(Naive_result_1_Z$Full_table)
Naive_result_3_table<-as.table(Naive_result_3_Z$Full_table)
Naive_result_5_table<-as.table(Naive_result_5_Z$Full_table)
ZZ1=rbind(Naive_result_1_table,Naive_result_3_table,Naive_result_5_table)

#table for 0.30.3Z
Zc1_table<-as.table(replicate_with_estimate_Z_1$Full_table)
Zc3_table<-as.table(replicate_with_estimate_Z_3$Full_table)
Zc5_table<-as.table(replicate_with_estimate_Z_5$Full_table)
ZZ3=rbind(Zc1_table,Zc3_table,Zc5_table)

ZZ_Full=cbind(ZZ1,ZZ2,ZZ3)
colnames(ZZ_Full) <- c("beta0_ls", "beta1_ls","betaZ_ls", "beta0_real","beta1_real","betaZ_real","beta_0","beta_1","beta_Z")
ss_Full

ZZ_Full_latex <- xtable(cbind(dist,size,lab,ZZ_Full),digits = 3)
print(ZZ_Full_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)

#table for 0.10.3Z
Z2_table<-as.table(replicate_with_true_Z_2$Full_table)
Z4_table<-as.table(replicate_with_true_Z_4$Full_table)
Z6_table<-as.table(replicate_with_true_Z_6$Full_table)
ZZ_2=rbind(Z2_table,Z4_table,Z6_table)

Naive_result_2_table<-as.table(Naive_result_2_Z$Full_table)
Naive_result_4_table<-as.table(Naive_result_4_Z$Full_table)
Naive_result_6_table<-as.table(Naive_result_6_Z$Full_table)
ZZ_1=rbind(Naive_result_2_table,Naive_result_4_table,Naive_result_6_table)

#table for 0.10.3Z
Zc2_table<-as.table(replicate_with_estimate_Z_2$Full_table)
Zc4_table<-as.table(replicate_with_estimate_Z_4$Full_table)
Zc6_table<-as.table(replicate_with_estimate_Z_6$Full_table)
ZZ_3=rbind(Zc2_table,Zc4_table,Zc6_table)

ZZ01_Full=cbind(ZZ_1,ZZ_2,ZZ_3)
colnames(ZZ01_Full) <- c("beta0_ls", "beta1_ls","betaZ_ls", "beta0_real","beta1_real","betaZ_real","beta_0","beta_1","beta_Z")

ZZ01_Full_latex <- xtable(cbind(dist,size,lab,ZZ01_Full),digits = 3)
print(ZZ01_Full_latex, tabular.environment = "longtable", floating = FALSE,include.rownames = F)
