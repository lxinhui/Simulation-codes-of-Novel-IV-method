setwd("D:/")
library(ggplot2)
library(readxl)
#----------------------------------------------vary across N-----------------------------------------------------------
N=seq(400,5000,100)
D_N <- read.csv("general-N1_use.csv")

tiff(file = "aaaageneral-N.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.4,0.3),xlab="Sample size",ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.22),xlab="Sample size",ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.13),xlab="Sample size",ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()

#----------------------------------------------vary across GX-----------------------------------------------------------
N=seq(0.25,0.9,0.05)#seq(50,5000,100)
D_N <- read.csv("general-a11_use.csv")

tiff(file = "aaaageneral-GX.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.4,0.3),xlab="The Causal Effect of G on X",ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.22),xlab="The Causal Effect of G on X",ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.13),xlab="The Causal Effect of G on X",ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()
#----------------------------------------------vary across XX1 slope-----------------------------------------------------------
N=seq(0.2,0.8,0.05)#seq(50,5000,100)
D_N <- read.csv("general-a111_use.csv")

tiff(file = "aaaageneral-xx1_slope����.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.4,1.6),xlab=expression(paste("The Slope of The Regression of ",X["1"]," on X")),ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.4),xlab=expression(paste("The Slope of The Regression of ",X["1"]," on X")),ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,1.8),xlab=expression(paste("The Slope of The Regression of ",X["1"]," on X")),ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()
#----------------------------------------------vary across XY-----------------------------------------------------------
N=seq(0.2,0.75,0.05)#seq(50,5000,100)
D_N <- read.csv("general-b11_use.csv")

tiff(file = "aaaageneral-XY.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.5,0.3),xlab="The Causal Effect of X on Y",ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.1),xlab="The Causal Effect of X on Y",ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.2),xlab="The Causal Effect of X on Y",ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()
#----------------------------------------------vary across XX1 intercept-----------------------------------------------------------
N=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)#seq(50,5000,100)
D_N <- read.csv("general-c11_use.csv")

tiff(file = "aaaageneral-xx1_int.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.4,0.35),xlab=expression(paste("The Intercept of The Regression of ",X["1"]," on X")),ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.07),xlab=expression(paste("The Intercept of The Regression of ",X["1"]," on X")),ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.14),xlab=expression(paste("The Intercept of The Regression of ",X["1"]," on X")),ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()

#----------------------------------------------vary across yy1 slope-----------------------------------------------------------
N=seq(0.2,1,0.05)#seq(50,5000,100)
D_N <- read.csv("general-a21_use.csv")

tiff(file = "aaaageneral-yy1_slope.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.4,0.6),xlab=expression(paste("The Slope of The Regression of ",Y["1"]," on Y")),ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.1),xlab=expression(paste("The Slope of The Regression of ",Y["1"]," on Y")),ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.25),xlab=expression(paste("The Slope of The Regression of ",Y["1"]," on Y")),ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()

#----------------------------------------------vary across YY1 intercept-----------------------------------------------------------
N=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5)#seq(50,5000,100)
D_N <- read.csv("general-c21_use.csv")

tiff(file = "aaaageneral-YY1_int.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.3,0.65),xlab=expression(paste("The Intercept of The Regression of ",Y["1"]," on Y")),ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.07),xlab=expression(paste("The Intercept of The Regression of ",Y["1"]," on Y")),ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.22),xlab=expression(paste("The Intercept of The Regression of ",Y["1"]," on Y")),ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()
#----------------------------------------------vary across p-----------------------------------------------------------
N=seq(0.05,0.5,0.05)#seq(50,5000,100)
D_N <- read.csv("general-p_use.csv")

tiff(file = "aaaageneral-p.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.4,0.4),xlab="Minor allele frequency of G",ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.17),xlab="Minor allele frequency of G",ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.14),xlab="Minor allele frequency of G",ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()

#----------------------------------------------vary across s-----------------------------------------------------------
N=seq(0.35,0.6,0.05)#seq(50,5000,100)
D_N <- read.csv("general-s_use.csv")

tiff(file = "aaaageneral-s.tiff", width = 3000, height = 1300, units = "px", res = 200)
split.screen(c(1,3))
screen(1)
plot(D_N$bias55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(-0.25,0.4),xlab="Correlation between W and measurement values",ylab = "Bias",cex=0.8)
lines(D_N$bias22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$bias33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$bias44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(2)
plot(D_N$se55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.07),xlab="Correlation between W and measurement values",ylab = "SE",cex=0.8)
lines(D_N$se22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$se33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$se44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

screen(3)
plot(D_N$mse55~N,type="b",col="darkgoldenrod1",pch=1,lwd=2,ylim = c(0,0.06),xlab="Correlation between W and measurement values",ylab = "MSE",cex=0.8)
lines(D_N$mse22~N,type="b", col="firebrick",pch=2,lwd=2,cex=0.8)
lines(D_N$mse33~N,type="b", col="royalblue3",pch=3,lwd=2,cex=0.8)
lines(D_N$mse44~N,type="b", col="chartreuse4",pch=4,lwd=2,cex=0.8)
abline(h=0,lty=3,lwd=2)
legend("topleft",c("Crude Association","IV","GMM-3SLS","Novel IV"),pch=c(4,3,1,2),col=c("chartreuse4","royalblue3","darkgoldenrod1","firebrick"))

dev.off()