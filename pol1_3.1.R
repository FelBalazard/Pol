#####Evaluation####
library(MASS)#multivariate normal

#Coefficients for linear model Y =b1*X+b2*T+b3*X*T+eps
b1=1
b2=1
b3=1.3
sig=1#variance of eps

#Sample size
Nsample=300

#Number of simulations
nsim=10000

#Storage of quantities of interest
Theta=numeric(nsim)#The theoretical improvement
Theta_hat=numeric(nsim)#The estimated improvement (mean)
Theta.low=numeric(nsim)#The lower bound on improvement
linreg=numeric(nsim)#p-value of interaction test in linear regression

#Sampling the covariate
set.seed(12)#The best seed
X=2*runif(Nsample,0,1)-1#The covariate is fixed.

#Treatment allocation
Tr=rep(c(0,1),Nsample/2)
for(sim in 1:nsim){
  set.seed(13+sim)
  eps=rnorm(Nsample,0,sig)#noise in linear model
  Y=b1*X+b2*Tr+b3*X*Tr+eps#generation of outcomes following the linear model
  reglin=lm(Y~X*Tr)
  beta_var=vcov(reglin)#covariance matrix for simulation
  
  #Theoretical improvement on distribution of X. Depends on the estimated policy.
  delta_theor=b2+b3*X
  Theta[sim]=-1/Nsample*sum(delta_theor[delta<0])
  
  #Estimation procedure 
  NMCMC=10000#number of sampling of beta*
  betaci=mvrnorm(NMCMC,reglin$coefficients[3:4],beta_var[3:4,3:4])#Sampling of beta*
  deltaci=betaci[,2]%*%t(X)+betaci[,1]#Delta*
  Thetaci=-1/Nsample*sapply(1:NMCMC,function(i) sum(deltaci[i,which(delta<0)]))#Theta*
  Theta_hat[sim]=mean(Thetaci)
  Theta.low[sim]=quantile(Thetaci,0.05)
  linreg[sim]=summary(lm(Y~X*Tr))$coefficient[4,4]
}

#####Plot####
library(tikzDevice)#Allows to have beautiful latex type setting in axis labels. <3
#This is to plot figure 1 that illustrates the problem with using policy 1.

#The following line is commented out to prevent unintentional overwriting of figure.
#tikz(file = "~/Documents/publications/personalized medecine/pol1.tex", width = 5, height = 8)
par(mfrow=c(2,1),mar=c(4,4.5,1,.5))
plot(Theta[c(1:1000)],Theta_hat[c(1:1000)],xlab="$\\Theta(\\hat{\\texttt{pol}}_1)$",
     ylab="$\\hat\\Theta(\\hat{\\texttt{pol}}_1)$",pch=3)
abline(a=0,b=1)
plot(Theta[c(1:1000)],Theta.low[c(1:1000)],xlab="$\\Theta(\\hat{\\texttt{pol}}_1)$",
     ylab="$\\hat{q}_{n,0.05}(\\hat{\\texttt{pol}}_1)$",pch=3)
#points(Theta[linreg<0.05],Theta.low[linreg<0.05],col=2,pch=3)
abline(h=0)
abline(a=0,b=1)
abline(v=0)
#The following line needs to be executed to finish the figure file.
#dev.off()

#####Numerical quantities####
#Power (size)
mean(linreg<0.05)
mean(Theta.low>0)

#Coverage
mean(Theta>Theta.low)#underestimation
mean(Theta<Theta.low)#overestimation
mean(Theta>Theta_hat)#underestimation
mean(Theta<Theta_hat)#overestimation

#False positives
mean(Theta[Theta.low>0]<0)

#Personalized proportion
mean(Pneg)
mean(Pneg[Theta.low>0])
