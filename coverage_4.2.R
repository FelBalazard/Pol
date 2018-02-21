#####Evaluation####
library(MASS)#multivariate normal

#Coefficients for linear model Y =b1*X+b2*T+b3*X*T+eps
b1=1

sig=1#variance of eps

#Sample size
Nsample=300

#Number of simulations
nsim=10000

#Thresholds 
thresh=0:5/10+qnorm(0.05)
npol=length(thresh)

#Storage of quantities of interest
Theta=matrix(numeric(npol*nsim),nsim,npol) # The theoretical improvement
Theta_hat=matrix(numeric(npol*nsim),nsim,npol) #The estimated improvement (median)
Theta.low=matrix(numeric(npol*nsim),nsim,npol) #The lower bound on improvement
Theta.low2=matrix(numeric(npol*nsim),nsim,npol) #The lower bound on improvement
Theta.up=matrix(numeric(npol*nsim),nsim,npol) #Another quantile
Pneg=matrix(numeric(npol*nsim),nsim,npol) #Proportion of individuals personalized
Cor_DT=matrix(numeric(npol*nsim),nsim,npol)#correlation between deltaci and thetaci
linreg=numeric(nsim) #p-value of interaction test in linear regression.
minz=numeric(nsim)#smallest standardized difference

set.seed(12)
X=2*runif(Nsample,0,1)-1#The covariate is fixed.

#Treatment allocation
Tr=rep(c(0,1),Nsample/2)
beta2=c(0.3)
beta3=c(0.5,0.8,1)
resultcovunderlow=numeric(length(beta3))
resultcoveqlow=numeric(length(beta3))
resultcovoverlow=numeric(length(beta3))
resultcovunderlow2=numeric(length(beta3))
resultcoveqlow2=numeric(length(beta3))
resultcovoverlow2=numeric(length(beta3))
resultcovundermedian=numeric(length(beta3))
resultcoveqmedian=numeric(length(beta3))
resultcovovermedian=numeric(length(beta3))
resultcovunderup=numeric(length(beta3))
resultcovequp=numeric(length(beta3))
resultcovoverup=numeric(length(beta3))

for(i in 1:length(beta2)){
  for(j in 1:length(beta3)){
b2=beta2[i]
b3=beta3[j]
for(sim in 1:nsim){
  set.seed(13+sim)
  eps=rnorm(Nsample,0,sig)#noise in linear model
  Y=b1*X+b2*Tr+b3*X*Tr+eps#generation of outcomes following the linear model
  reglin=lm(Y~X*Tr)
  beta_var=vcov(reglin)#covariance matrix for simulation
  
  #Policies definition (depending on the threshold on z_Delta)
  zscore=(reglin$coefficients[4]*X+reglin$coefficients[3])/sqrt(beta_var[3,3]+2*beta_var[4,3]*X+X^2*beta_var[4,4])
  policy=sapply(1:npol,function(i) zscore<thresh[i])
  
  #Theoretical improvement on distribution of X. Depends on the estimated policy.  
  delta_theor=b2+b3*X
  Theta[sim,]=-1/Nsample*delta_theor%*%policy
  
  #Estimation procedure 
  NMCMC=10000#number of sampling of beta*
  betaci=mvrnorm(NMCMC,reglin$coefficients[3:4],beta_var[3:4,3:4])#Sampling of beta*
  deltaci=(betaci[,2]%*%t(X)+betaci[,1])#Delta*
  Thetaci=-1/Nsample*t(sapply(1:NMCMC,function(i) deltaci[i,]%*%policy))#Theta*
  
  Theta_hat[sim,]=apply(Thetaci,2,mean)
  Theta.low[sim,]=sapply(1:npol,function(pol) quantile(Thetaci[,pol],0.05))
  Pneg[sim,]=apply(policy==T,2,mean)
  linreg[sim]=summary(lm(Y~X*Tr))$coefficient[4,4]
  minz[sim]=min(zscore)
}

#Max lower bound policy
policylow=sapply(1:nsim,function(i) min(which(Theta.low[i,]==max(Theta.low[i,]))))
resultcovunderlow[j]=mean(Theta[cbind(1:nsim,policylow)]>Theta.low[cbind(1:nsim,policylow)])
resultcovoverlow[j]=mean(Theta[cbind(1:nsim,policylow)]<Theta.low[cbind(1:nsim,policylow)])
resultcoveqlow[j]=mean(Theta[cbind(1:nsim,policylow)]==Theta.low[cbind(1:nsim,policylow)])
resultcovunderlow2[j]=mean(Theta[cbind(1:nsim,policylow)]>Theta.low2[cbind(1:nsim,policylow)])
resultcovoverlow2[j]=mean(Theta[cbind(1:nsim,policylow)]<Theta.low2[cbind(1:nsim,policylow)])
resultcoveqlow2[j]=mean(Theta[cbind(1:nsim,policylow)]==Theta.low2[cbind(1:nsim,policylow)])
resultcovundermedian[j]=mean(Theta[cbind(1:nsim,policylow)]>Theta_hat[cbind(1:nsim,policylow)])
resultcovovermedian[j]=mean(Theta[cbind(1:nsim,policylow)]<Theta_hat[cbind(1:nsim,policylow)])
resultcoveqmedian[j]=mean(Theta[cbind(1:nsim,policylow)]==Theta_hat[cbind(1:nsim,policylow)])
resultcovunderup[j]=mean(Theta[cbind(1:nsim,policylow)]>Theta.up[cbind(1:nsim,policylow)])
resultcovoverup[j]=mean(Theta[cbind(1:nsim,policylow)]<Theta.up[cbind(1:nsim,policylow)])
resultcovequp[j]=mean(Theta[cbind(1:nsim,policylow)]==Theta.up[cbind(1:nsim,policylow)])

  }
}


#####Plots####
#This is to plot figure 4
library(tikzDevice)

#tikz(file = "~/Documents/publications/personalized medecine/coverage.tex", width = 6, height = 6)
par(mfrow=c(2,2),mar=c(4,2,3,.5))
barplot(rbind(resultcovunderlow,resultcoveqlow,resultcovoverlow),
        col=c(4,3,2),xlab="$\\beta_3$",space=rep(0,3),main="$\\hat{q}_{n,0.05}$")
abline(h=0.95)
axis(1,at=1:3-1/2,labels=c(0.5,0.8,1))

barplot(rbind(resultcovunderlow2,resultcoveqlow2,resultcovoverlow2),
        col=c(4,3,2),xlab="$\\beta_3$",space=rep(0,3),main="$\\hat{q}_{n,0.25}$")
abline(h=0.75)
axis(1,at=1:3-1/2,labels=c(0.5,0.8,1))

barplot(rbind(resultcovundermedian,resultcoveqmedian,resultcovovermedian),
        col=c(4,3,2),xlab="$\\beta_3$",space=rep(0,3),main="$\\hat\\Theta=\\hat{q}_{n,0.5}$")
abline(h=0.5)
axis(1,at=1:3-1/2,labels=c(0.5,0.8,1))


barplot(rbind(resultcovunderup,resultcovequp,resultcovoverup),
        col=c(4,3,2),xlab="$\\beta_3$",space=rep(0,3),main="$\\hat{q}_{n,0.95}$")
abline(h=0.05)
axis(1,at=1:3-1/2,labels=c(0.5,0.8,1))
#dev.off()
