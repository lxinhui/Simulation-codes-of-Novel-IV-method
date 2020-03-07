setwd("D:/")
library(mvtnorm)
library(MASS)
options(scipen = 200)
set.seed(44)

discrimination_error<-function(NN,a_gx,a_xy,px1,py1){
  
  p1 <- 0.5
  N <- 10000#TIMES
  
  
  G <- rep(0,NN)
  U <- rep(0,NN)
  X <- rep(0,NN)
  Y <- rep(0,NN)
  X1 <- rep(0,NN)
  X2 <- rep(0,NN)
  Y1 <- rep(0,NN)
  Y2 <- rep(0,NN)
  
  bias12 <- rep(0,N)
  bias13 <- rep(0,N)
  effect <- rep(0,N)
  effect_est <- rep(0,N)
  effect_est2 <- rep(0,N)
  
  effect_est3 <- rep(0,N)
  effect_est4 <- rep(0,N)
  bias14 <- rep(0,N)
  bias15 <- rep(0,N)
  
  
  for(MM in 1:N){
    for(KK in 1:NN){
      G[KK]<-rbinom(1,1,0.5)
      U[KK]<-rbinom(1,1,0.5) 
      
      px <- exp(0.05+a_gx*G[KK]+0.5*U[KK])/(1+exp(0.05+a_gx*G[KK]+0.5*U[KK]))
      X[KK] <- rbinom(1,1,px)
      
      py <- exp(0.08+a_xy*X[KK]+0.8*U[KK])/(1+exp(0.08+a_xy*X[KK]+0.8*U[KK]))
      Y[KK] <- rbinom(1,1,py)
      
    }
    data <- data.frame(G,U,X,Y)
   
    ###########THE GENERATION OF X1,X2,Y1,Y2##########
    
    aa <- which(data$X==1)
    aa_sub <- sample(aa,round(length(aa)*px1))
    aa0 <- which(data$X==0)
    aa_sub0 <- sample(aa0,round(length(aa0)*px1)) 
    data$X1 <- data$X
    for(z in aa_sub){
      data[z,which(colnames(data)=="X1")] <- 0
    }
    for(z in aa_sub0){
      data[z,which(colnames(data)=="X1")] <- 1
    }
    
    aa <- which(data$X==1)
    aa_sub <- sample(aa,round(length(aa)*px1))
    aa0 <- which(data$X==0)
    aa_sub0 <- sample(aa0,round(length(aa0)*px1)) 
    data$X2 <- data$X
    for(z in aa_sub){
      data[z,which(colnames(data)=="X2")] <- 0
    }
    for(z in aa_sub0){
      data[z,which(colnames(data)=="X2")] <- 1
    }
    
    aa <- which(data$Y==1)
    aa_sub <- sample(aa,round(length(aa)*py1))
    aa0 <- which(data$Y==0)
    aa_sub0 <- sample(aa0,round(length(aa0)*py1)) 
    data$Y1 <- data$Y
    for(z in aa_sub){
      data[z,which(colnames(data)=="Y1")] <- 0
    }
    for(z in aa_sub0){
      data[z,which(colnames(data)=="Y1")] <- 1
    }
    
    aa <- which(data$Y==1)
    aa_sub <- sample(aa,round(length(aa)*py1))
    aa0 <- which(data$Y==0)
    aa_sub0 <- sample(aa0,round(length(aa0)*py1)) 
    data$Y2 <- data$Y
    for(z in aa_sub){
      data[z,which(colnames(data)=="Y2")] <- 0
    }
    for(z in aa_sub0){
      data[z,which(colnames(data)=="Y2")] <- 1
    }
    
    
  ##########GODE STANDARD##########
    A <- 0
    for(XX in c(0,1)){#s:=X,t:=U
      for(UU in c(0,1)){
        
        P_Y_XU <- exp(0.08+a_xy*XX+0.8*UU)/(1+exp(0.08+a_xy*XX+0.8*UU))
        
        if(XX==1){
          P_X_UZ <- exp(0.05+a_gx*1+0.5*UU)/(1+exp(0.05+a_gx*1+0.5*UU))
        }
        else{
          P_X_UZ <- 1-exp(0.05+a_gx*1+0.5*UU)/(1+exp(0.05+a_gx*1+0.5*UU))
        }
        if(UU==1){
          P_U <- 0.5
        }
        else{
          P_U <- 0.5
        }
        A <- A+P_Y_XU*P_X_UZ*P_U
      }
    }
    B <- 0
    for(XX in c(0,1)){#s:=X,t:=U
      for(UU in c(0,1)){
        
        P_Y_XU <- exp(0.08+a_xy*XX+0.8*UU)/(1+exp(0.08+a_xy*XX+0.8*UU))
        
        if(XX==1){
          P_X_UZ <- exp(0.05+a_gx*0+0.5*UU)/(1+exp(0.05+a_gx*0+0.5*UU))
        }
        else{
          P_X_UZ <- 1-exp(0.05+a_gx*0+0.5*UU)/(1+exp(0.05+a_gx*0+0.5*UU))
        }
        if(UU==1){
          P_U <- 0.5
        }
        else{
          P_U <- 0.5
        }
        B <- B+P_Y_XU*P_X_UZ*P_U
      }
    }
    C <- 0
    for(UU in c(0,1)){
      
      P_X_UZ <- exp(0.05+a_gx*1+0.5*UU)/(1+exp(0.05+a_gx*1+0.5*UU))
      
      if(UU==1){
        P_U <- 0.5
      }
      else{
        P_U <- 0.5
      }
      C <- C+P_X_UZ*P_U
    }
    D <- 0
    for(UU in c(0,1)){
      
      P_X_UZ <- exp(0.05+a_gx*0+0.5*UU)/(1+exp(0.05+a_gx*0+0.5*UU))
      
      if(UU==1){
        P_U <- 0.5
      }
      else{
        P_U <- 0.5
      }
      D <- D+P_X_UZ*P_U
    }
    effect[MM] <- (A-B)/(C-D)
    
    ########################MODEL1:NOVEL IV########################
    data2 <- data
    data2$a <- ifelse(data2$X1==1&data2$X2==1&data2$G==1,1,0)#number(X=1,G=0)
    data2$b <- ifelse(data2$X1==1&data2$X2==1&data2$G==0,1,0)#number(X=1,G=0)
    data2$c <- ifelse(data2$X1+data2$X2==1&data2$G==1,1,0)#number(X=1,G=0)
    data2$d <- ifelse(data2$X1+data2$X2==1&data2$G==0,1,0)#number(X=1,G=0)
    data2$e <- ifelse(data2$X1==0&data2$X2==0&data2$G==1,1,0)#number(X=1,G=0)
    data2$f <- ifelse(data2$X1==0&data2$X2==0&data2$G==0,1,0)#number(X=1,G=0)
    
    data2$g <- ifelse(data2$Y1==1&data2$Y2==1&data2$G==1,1,0)#number(X=1,G=0)
    data2$h <- ifelse(data2$Y1==1&data2$Y2==1&data2$G==0,1,0)#number(X=1,G=0)
    data2$i <- ifelse(data2$Y1+data2$Y2==1&data2$G==1,1,0)#number(X=1,G=0)
    data2$j <- ifelse(data2$Y1+data2$Y2==1&data2$G==0,1,0)#number(X=1,G=0)
    data2$k <- ifelse(data2$Y1==0&data2$Y2==0&data2$G==1,1,0)#number(X=1,G=0)
    data2$l <- ifelse(data2$Y1==0&data2$Y2==0&data2$G==0,1,0)#number(X=1,G=0)
    
    a <- sum(data2$a)
    b <- sum(data2$b)
    c <- sum(data2$c)
    d <- sum(data2$d)
    e <- sum(data2$e)
    f <- sum(data2$f)
    
    g <- sum(data2$g)
    h <- sum(data2$h)
    i <- sum(data2$i)
    j <- sum(data2$j)
    k <- sum(data2$k)
    l <- sum(data2$l)
    
    
    alpha <- 0.5+0.5*sqrt((NN-2*(c+d))/NN)   
    beta <- 0.5+0.5*sqrt((NN-2*(i+j))/NN) 

    C11 <- (a*(alpha^2)-e*((1-alpha)^2))/((2*alpha-1)*(a+e))
    D11 <- (b*(alpha^2)-f*((1-alpha)^2))/((2*alpha-1)*(b+f))
    A11 <- (g*(beta^2)-k*((1-beta)^2))/((2*beta-1)*(g+k))
    B11 <- (h*(beta^2)-l*((1-beta)^2))/((2*beta-1)*(h+l))
    
    effect_est2[MM] <- (A11-B11)/(C11-D11)
    bias13[MM] <- effect_est2[MM]-effect[MM]
    
    ########################MODEL2:IV########################
    data21 <- data[,c(1:4,5,7)]
    data22 <- data[,c(1:4,6,8)]
    colnames(data22) <- colnames(data21)
    data2 <- rbind(data21,data22)
    data2$x1g_11 <- ifelse(data2$X1==1&data2$G==1,1,0)#number(X1=1,G=1)
    data2$x1g_10 <- ifelse(data2$X1==1&data2$G==0,1,0)#number(X1=1,G=1)
    data2$p_G_1 <- ifelse(data2$G==1,1,0)#number(G=1)
    data2$p_G_0 <- ifelse(data2$G==0,1,0)#number(G=1)
    
    C1 <- (sum(data2$x1g_11)/NN)/(sum(data2$p_G_1)/NN)
    D1 <- (sum(data2$x1g_10)/NN)/(sum(data2$p_G_0)/NN)
    #
    data2$x1g_11 <- ifelse(data2$Y1==1&data2$G==1,1,0)#number(X1=1,G=1)
    data2$x1g_10 <- ifelse(data2$Y1==1&data2$G==0,1,0)#number(X1=1,G=1)
    data2$p_G_1 <- ifelse(data2$G==1,1,0)#number(G=1)
    data2$p_G_0 <- ifelse(data2$G==0,1,0)#number(G=1)
    
    A1 <- (sum(data2$x1g_11)/NN)/(sum(data2$p_G_1)/NN)
    B1 <- (sum(data2$x1g_10)/NN)/(sum(data2$p_G_0)/NN)
    
    
    effect_est3[MM] <- (A1-B1)/(C1-D1)
    bias14[MM] <- effect_est3[MM]-effect[MM]
    
    #########################MODEL3:CRUDE ASSOCIATION########################
    data2$y1x1_11 <- ifelse(data2$Y1==1&data2$X1==1,1,0)#number(X1=1,G=1)
    data2$y1x1_10 <- ifelse(data2$Y1==1&data2$X1==0,1,0)#number(X1=1,G=1)
    data2$p_x1_1 <- ifelse(data2$X1==1,1,0)#number(G=1)
    data2$p_x1_0 <- ifelse(data2$X1==0,1,0)#number(G=1)
    
    A1 <- (sum(data2$y1x1_11)/NN)/(sum(data2$p_x1_1)/NN)
    B1 <- (sum(data2$y1x1_10)/NN)/(sum(data2$p_x1_0)/NN)
    
    
    effect_est4[MM] <- A1-B1
    bias15[MM] <- effect_est4[MM]-effect[MM]
  }
  

  bias33 <- mean(bias13)
  bias44 <- mean(bias14)
  bias55 <- mean(bias15)
  

  se33<-sd(effect_est2)
  se44<-sd(effect_est3)
  se55<-sd(effect_est4)
  
  mse33<-bias33^2+se33^2
  mse44<-bias44^2+se44^2
  mse55<-bias55^2+se55^2
  
  aa <- list(bias33,bias44,bias55,se33,se44,se55,mse33,mse44,mse55)
  names(aa) <- c("bias22","bias33","bias44","se22","se33","se44","mse22","mse33","mse44")
  return(aa)
}


#----------------------------------------------SAMPLE SIZE-----------------------------------------------------------
N=seq(200,5000,200)
D_N<-data.frame(matrix(NA,0,9))
colnames(D_N) <- c("bias22","bias33","bias44","se22","se33","se44","mse22","mse33","mse44")

for(i in 1:length(N)){
  cat(i)
  cau.hat1=discrimination_error(N[i],1.6,1.4,0.1,0.12)
  D_N<-rbind(D_N,cau.hat1)
}
write.csv(D_N,"general-N1_use.csv",row.names=F)

#----------------------------------------------CAUSAL EFFECT OF G ON X----------------------------------------------------
a1=seq(0.3,2,0.1)#0.15,0.9,0.05
D_a1<-data.frame(matrix(NA,0,9))
colnames(D_a1) <- c("bias22","bias33","bias44","se22","se33","se44","mse22","mse33","mse44")

for(i in 1:length(a1)){
  cau.hat1=discrimination_error(5000,a1[i],1.4,0.1,0.12)
  D_a1<-rbind(D_a1,cau.hat1)
}
write.csv(D_a1,"general-a11_use.csv",row.names=F)
#----------------------------------------------CAUSAL EFFECT OF X ON Y-----------------------------------------------------------
bbb=seq(0.1,1.7,0.1)
D_b<-data.frame(matrix(NA,0,9))
colnames(D_b) <- c("bias22","bias33","bias44","se22","se33","se44","mse22","mse33","mse44")

for(i in 1:17){
  cau.hat1=discrimination_error(5000,1.6,bbb[i],0.1,0.12)
  D_b<-rbind(D_b,cau.hat1)
}
write.csv(D_b,"general-b11_use.csv",row.names=F)
#----------------------------------------------MISCLASSIFICATION OF X1-----------------------------------------------------------
a2=seq(0.02,0.26,0.02)

D_a2<-data.frame(matrix(NA,0,9))
colnames(D_a2) <- c("bias22","bias33","bias44","se22","se33","se44","mse22","mse33","mse44")

for(i in 1:13){
  cau.hat1=discrimination_error(5000,1.6,1.4,a2[i],0.12)
  D_a2<-rbind(D_a2,cau.hat1)
}
write.csv(D_a2,"general-x_use.csv",row.names=F)

#----------------------------------------------MISCLASSIFICATION OF Y1-----------------------------------------------------------
c1=seq(0.02,0.3,0.02)
D_c1<-data.frame(matrix(NA,0,9))
colnames(D_c1) <- c("bias22","bias33","bias44","se22","se33","se44","mse22","mse33","mse44")

for(i in 1:15){
  cat(i)
  cau.hat1=discrimination_error(5000,1.6,1.4,0.1,c1[i])
  D_c1<-rbind(D_c1,cau.hat1)
}
write.csv(D_c1,"general-y_use.csv",row.names=F)
