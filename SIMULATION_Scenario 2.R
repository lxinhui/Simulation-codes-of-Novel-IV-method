setwd("D:/")
library(mvtnorm)
library(MASS)
library(systemfit)
options(scipen = 200)
set.seed(44)

discrimination_error<-function(a1,b,NN,a2,a3,a6,a7,c1,c2,c3,c4,p,s){

  N<-10000#TIMES
  
  bias5<-rep(0,N)
  bias2<-rep(0,N)
  bias3<-rep(0,N)
  bias4<-rep(0,N)
  beta5<-rep(0,N)
  beta2<-rep(0,N)
  beta3<-rep(0,N)
  beta4<-rep(0,N)
  
  
  for(kk in 1:N){
 
    MM11<-rep(0,NN)
    MM12<-rep(0,NN)
    MM31<-rep(0,NN)
    MM32<-rep(0,NN)
    MM41<-rep(0,NN)
    MM42<-rep(0,NN)
    
    a4 <- 0.4
    a5 <- 0.4
    

    G<-rbinom(NN,2,p)
    U<-rnorm(NN,0,1) #U
    W <- rnorm(NN,0,1)
    
    e1<-rnorm(NN,0,sqrt(1-(2*p*(1-p))*a1^2-a4^2))
    e2<-rnorm(NN,0,sqrt(1-a5^2-b^2-2*b*a5*a4))
    M1<-a1*G+a4*U+e1 #X
    M2<-b*M1+a5*U+e2 #Y
    e3<-rnorm(NN,-c1,sqrt(1-a2^2-s^2)) 
    e4<-rnorm(NN,-c2,sqrt(1-a3^2-s^2))
    e5<-rnorm(NN,-c3,sqrt(1-a6^2-s^2))
    e6<-rnorm(NN,-c4,sqrt(1-a7^2-s^2))
    
    M11<-c1+a2*M1+s*W+e3#measured value X1
    M12<-c2+a3*M1+s*W+e4#measured value X2
    M21<-c3+a6*M2+s*W+e5#measured value Y1
    M22<-c4+a7*M2+s*W+e6#measured value Y2
    
    #MODEL1: GMM-3SLS
    MM31<-M11
    MM32<-M21
    MM33<-M12
    MM34<-M22
    data31<-as.data.frame(cbind(G,MM31,MM32,U))
    data32<-as.data.frame(cbind(G,MM33,MM34,U))
    colnames(data32) <- c("G","MM31","MM32","U")
    data3<-rbind(data31,data32)
    
    attach(data3)
    equate <- MM32~MM31
    fit3sls <- systemfit(equate, method = "3SLS", method3sls= "GMM",inst = ~G)
    detach(data3)
    
    beta5[kk]<-as.numeric(fit3sls$coefficients[2])
    bias5[kk]<-beta5[kk]-b
    
    #MODEL2: MOVEL IV
    data2<-as.data.frame(cbind(G,M1,M2,M11,M12,M21,M22,U))
    
    R_GM21<-cor(data2$M21,data2$G)
    R_GM22<-cor(data2$M22,data2$G)
    R_M11M12<-cor(data2$M11,data2$M12)
    R_M21M22<-cor(data2$M21,data2$M22)
    R_M11M21<-cor(data2$M11,data2$M21)
    R_M12M22<-cor(data2$M12,data2$M22)
    R_GM11<-cor(data2$M11,data2$G)
    R_GM12<-cor(data2$M12,data2$G)
    
    if(R_GM11/R_GM21>0){x<-1}
    if(R_GM11/R_GM21<0){x<-(-1)}
    
    ss <- ((R_M11M21*R_GM12*R_GM22)-(R_M12M22*R_GM11*R_GM21))/((R_GM12*R_GM22)-(R_GM11*R_GM21))
    beta2[kk]<-x*(sqrt(R_GM21*(R_M11M12-ss)*R_GM22))/(sqrt(R_GM11*(R_M21M22-ss)*R_GM12))
    bias2[kk]<-beta2[kk]-b
    
    #MODEL3: IV
    MM31<-M11
    MM32<-M21
    MM33<-M12
    MM34<-M22
    data31<-as.data.frame(cbind(G,MM31,MM32,U))
    data32<-as.data.frame(cbind(G,MM33,MM34,U))
    colnames(data32) <- c("G","MM31","MM32","U")
    data3<-rbind(data31,data32)
    
    model32<-lm(data3$MM32~data3$G,data=data3)
    A32<-summary(model32)$coef[2,1]
    
    model31<-lm(data3$MM31~data3$G,data=data3)
    A31<-summary(model31)$coef[2,1]
    
    beta3[kk]<-A32/A31
    bias3[kk]<-beta3[kk]-b
    
    #MODEL4: CRUDE ASSOCIATION
    MM41<-M11
    MM42<-M21
    MM43<-M12
    MM44<-M22
    data41<-as.data.frame(cbind(G,MM41,MM42,U))
    data42<-as.data.frame(cbind(G,MM43,MM44,U))
    colnames(data42) <- c("G","MM41","MM42","U")
    data4<-rbind(data41,data42)
    
    model41<-lm(MM42~MM41,data=data4)
    
    beta4[kk]<-summary(model41)$coef[2,1]
    bias4[kk]<-beta4[kk]-b
  }
  

  bias22<-mean(bias2)
  bias33<-mean(bias3)
  bias44<-mean(bias4)
  bias55<-mean(bias5)

  se22<-sd(beta2)
  se33<-sd(beta3)
  se44<-sd(beta4)
  se55<-sd(beta5)

  mse22<-bias22^2+se22^2
  mse33<-bias33^2+se33^2
  mse44<-bias44^2+se44^2
  mse55<-bias55^2+se55^2
  
  aa <- list(bias22,bias33,bias44,bias55,se22,se33,se44,se55,mse22,mse33,mse44,mse55)
  names(aa) <- c("bias22","bias33","bias44","bias55","se22","se33","se44","se55","mse22","mse33","mse44","mse55")
  return(aa)  
}
#----------------------------------------------CAUSAL EFFECT OF W ON X1-----------------------------------------------------------
s=seq(0.35,0.6,0.05)
D_N<-data.frame(matrix(NA,0,15))
colnames(D_N) <- c("bias11","bias22","bias33","bias44","bias55","se11","se22","se33","se44","se55","mse11","mse22","mse33","mse44","mse55")

for(i in 1:length(s)){
  cau.hat1=discrimination_error(0.8,0.6,5000,0.56,0.55,0.58,0.78,0.89,0.14,0.13,0.11,0.5,s[i])
  D_N<-rbind(D_N,cau.hat1)
}
write.csv(D_N,"general-s_use.csv",row.names=F)
