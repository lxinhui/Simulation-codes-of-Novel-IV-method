setwd("D:/")
library(mvtnorm)
library(MASS)
library(systemfit)
options(scipen = 200)
set.seed(44)

discrimination_error<-function(a1,b,NN,a2,a3,a6,a7,c1,c2,c3,c4,p){

  N<-10000#TIMES
  
  bias1<-rep(0,N)
  bias2<-rep(0,N)
  bias3<-rep(0,N)
  bias4<-rep(0,N)
  bias5<-rep(0,N)
  beta1<-rep(0,N)
  beta2<-rep(0,N)
  beta3<-rep(0,N)
  beta4<-rep(0,N)
  beta5<-rep(0,N)

  
  for(kk in 1:N){
    #DATA GENERATION 
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
    
    e1<-rnorm(NN,0,sqrt(1-(2*p*(1-p))*a1^2-a4^2))
    e2<-rnorm(NN,0,sqrt(1-a5^2-b^2-2*b*a5*a4))
    M1<-a1*G+a4*U+e1 #X
    M2<-b*M1+a5*U+e2 #Y
    e3<-rnorm(NN,-c1,sqrt(1-a2^2))
    e4<-rnorm(NN,-c2,sqrt(1-a3^2))
    e5<-rnorm(NN,-c3,sqrt(1-a6^2))
    e6<-rnorm(NN,-c4,sqrt(1-a7^2))
    
    M11<-c1+a2*M1+e3#measured value X1
    M12<-c2+a3*M1+e4#measured value X2
    M21<-c3+a6*M2+e5#measured value Y1
    M22<-c4+a7*M2+e6#measured value Y2
    
    #MODEL1:NOVEL IV
    data2<-as.data.frame(cbind(G,M1,M2,M11,M12,M21,M22,U))

    R_GM21<-cor(data2$M21,data2$G)
    R_GM22<-cor(data2$M22,data2$G)
    R_M11M12<-cor(data2$M11,data2$M12)
    R_M21M22<-cor(data2$M21,data2$M22)
    R_GM11<-cor(data2$M11,data2$G)
    R_GM12<-cor(data2$M12,data2$G)
    
    if(R_GM11/R_GM21>0){x<-1}
    if(R_GM11/R_GM21<0){x<-(-1)}
    
    beta2[kk]<-x*(R_GM21*sqrt(R_GM11*R_M11M12*R_GM22))/(R_GM11*sqrt(R_GM21*R_M21M22*R_GM12))
    bias2[kk]<-beta2[kk]-b
    
    #MODEL2: IV
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
    
    #MODEL3: CRUDE ASSOCIATION
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

    #MODEL4: GMM-3SLS
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
#----------------------------------------------SAMPLE SIZE-----------------------------------------------------------
N=seq(400,5000,100)
D_N<-data.frame(matrix(NA,0,15))
colnames(D_N) <- c("bias11","bias22","bias33","bias44","bias55","se11","se22","se33","se44","se55","mse11","mse22","mse33","mse44","mse55")

for(i in 1:length(N)){
  cat(i)
  cau.hat1=discrimination_error(0.8,0.6,N[i],0.56,0.55,0.58,0.78,0.89,0.14,0.13,0.11,0.5)
  D_N<-rbind(D_N,cau.hat1)
}
write.csv(D_N,"general-N1_use.csv",row.names=F)

#----------------------------------------------CAUSAL EFFECT G ON X----------------------------------------------------
a1=seq(0.25,0.9,0.05)#0.15,0.9,0.05
D_a1<-data.frame(matrix(NA,0,15))
colnames(D_a1) <- c("bias11","bias22","bias33","bias44","se11","se22","se33","se44","mse11","mse22","mse33","mse44")

for(i in 1:length(a1)){
  cau.hat1=discrimination_error(a1[i],0.6,5000,0.56,0.55,0.58,0.78,0.89,0.14,0.13,0.11,0.5)
  D_a1<-rbind(D_a1,cau.hat1)
}
write.csv(D_a1,"general-a11_use.csv",row.names=F)

#----------------------------------------------CAUSAL EFFECT X ON Y-----------------------------------------------------------
bbb=seq(0.2,0.75,0.05)
D_b<-data.frame(matrix(NA,0,15))
colnames(D_b) <- c("bias11","bias22","bias33","bias44","se11","se22","se33","se44","mse11","mse22","mse33","mse44")

for(i in 1:17){
  cau.hat1=discrimination_error(0.8,bbb[i],5000,0.56,0.55,0.58,0.78,0.89,0.14,0.13,0.11,0.5)
  D_b<-rbind(D_b,cau.hat1)
}
write.csv(D_b,"general-b11_use.csv",row.names=F)
#----------------------------------------------SLOPE OF X1-----------------------------------------------------------
a2=seq(0.2,0.8,0.05)

D_a2<-data.frame(matrix(NA,0,15))
colnames(D_a2) <- c("bias11","bias22","bias33","bias44","se11","se22","se33","se44","mse11","mse22","mse33","mse44")

for(i in 1:13){
  cau.hat1=discrimination_error(0.8,0.6,5000,a2[i],0.55,0.58,0.78,0.89,0.14,0.13,0.11,0.5)
  D_a2<-rbind(D_a2,cau.hat1)
}
write.csv(D_a2,"general-a111_use.csv",row.names=F)
#----------------------------------------------INTERCEPT OF X1-----------------------------------------------------------
c1=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
D_c1<-data.frame(matrix(NA,0,15))
colnames(D_c1) <- c("bias11","bias22","bias33","bias44","se11","se22","se33","se44","mse11","mse22","mse33","mse44")

for(i in 1:15){
  cau.hat1=discrimination_error(0.8,0.6,5000,0.56,0.55,0.58,0.78,c1[i],0.14,0.13,0.11,0.5)
  D_c1<-rbind(D_c1,cau.hat1)
}
write.csv(D_c1,"general-c11_use.csv",row.names=F)

#----------------------------------------------SLOPE OF Y1-----------------------------------------------------------
a6=seq(0.2,1,0.05)
D_a6<-data.frame(matrix(NA,0,15))
colnames(D_a6) <- c("bias11","bias22","bias33","bias44","bias55","se11","se22","se33","se44","se55","mse11","mse22","mse33","mse44","mse55")
for(i in 1:17){
  cau.hat1=discrimination_error(0.8,0.6,5000,0.56,0.55,a6[i],0.78,0.14,0.14,0.13,0.11,0.5)
  D_a6<-rbind(D_a6,cau.hat1)
}
write.csv(D_a6,"general-a21_use.csv",row.names=F)
#----------------------------------------------INTERCEPT OF Y1-------------------------------------------
c3=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)
D_c3<-data.frame(matrix(NA,0,15))
colnames(D_c3) <- c("bias11","bias22","bias33","bias44","bias55","se11","se22","se33","se44","se55","mse11","mse22","mse33","mse44","mse55")
for(i in 1:15){
  cau.hat1=discrimination_error(0.8,0.6,5000,0.56,0.55,0.89,0.78,0.14,0.14,c3[i],0.11,0.5)
  D_c3<-rbind(D_c3,cau.hat1)
}
write.csv(D_c3,"general-c21_use.csv",row.names=F)

#----------------------------------------------MAF-------------------------------------------
p=seq(0.05,0.5,0.05)
p1<-data.frame(matrix(NA,0,15))
colnames(p1) <- c("bias11","bias22","bias33","bias44","bias55","se11","se22","se33","se44","se55","mse11","mse22","mse33","mse44","mse55")
for(i in 1:11){
  cau.hat1=discrimination_error(0.8,0.6,5000,0.56,0.55,0.58,0.78,0.89,0.14,0.13,0.11,p[i])
  p1<-rbind(p1,cau.hat1)
}
write.csv(p1,"general-p_use.csv",row.names=F)

