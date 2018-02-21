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
nsim=2

#Thresholds 
thresh=qnorm(c(0.5,0.4,0.3,0.25,0.2,0.15,0.1,0.07,0.05,0.04,0.025,0.02))
npol=length(thresh)

#Storage of quantities of interest
Theta=matrix(numeric(npol*nsim),nsim,npol) # The theoretical improvement
Theta_hat=matrix(numeric(npol*nsim),nsim,npol) #The estimated improvement (median)
Theta.low=matrix(numeric(npol*nsim),nsim,npol) #The lower bound on improvement
Pneg=matrix(numeric(npol*nsim),nsim,npol) #Proportion of individuals personalized
Cor_DT=matrix(numeric(npol*nsim),nsim,npol)#correlation between deltaci and thetaci
linreg=numeric(nsim) #p-value of interaction test in linear regression.
minz=numeric(nsim)#Smallest standardized individual treatment effect

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
  Cor_DT[sim,]=sapply(1:npol,function(pol) if(sum(policy[,pol])==0){1}else{
    cor(-deltaci[,which(X==sign(reglin$coefficients[4])*max(sign(reglin$coefficients[4])*X[policy[,pol]]))]
        ,Thetaci[,pol])})
}

#Max lower bound policy
policylow=sapply(1:nsim,function(i) min(which(Theta.low[i,]==max(Theta.low[i,]))))

#####Plots####
library(tikzDevice)#Allows to have beautiful latex type setting in axis labels. <3

#1st plot behavior of Theta_hat and Theta.low against Theta with varying threshold.
#The following line is commented out to prevent unintentional overwriting of figure.
#tikz(file = "~/Documents/publications/personalized medecine/mlb.tex", width = 5, height = 8)
par(mfrow=c(2,1),mar=c(4,4.5,1,.5))
nplotb=1
nplot=nplotb+99
plot(Theta[nplotb:nplot,1],Theta_hat[nplotb:nplot,1],xlim=range(Theta[nplotb:nplot,]),
     ylim=range(Theta_hat[nplotb:nplot,]),pch=3,xlab="$\\Theta(\\texttt{pol})$",ylab="$\\hat\\Theta(\\texttt{pol})$")
points(Theta[nplotb:nplot,npol],Theta_hat[nplotb:nplot,npol],pch=2)
abline(a=0,b=1)
for (pol in 2:npol){
  segments(Theta[nplotb:nplot,pol-1],Theta_hat[nplotb:nplot,pol-1],Theta[nplotb:nplot,pol],
           Theta_hat[nplotb:nplot,pol])
}
abline(v=0)
points(Theta[cbind(1:nsim,policylow)][nplotb:nplot],Theta_hat[cbind(1:nsim,policylow)][nplotb:nplot],pch=4,col=2)

#The following plot is Theta_hat against Theta.low.
# plot(Theta.low[nplotb:nplot,1],Theta_hat[nplotb:nplot,1],
#      xlim=range(Theta.low[nplotb:nplot,]),ylim=range(Theta_hat[nplotb:nplot,]),pch=3,
#      xlab="$\\hat\\Theta$(pol,0.05)",ylab="$\\hat\\Theta$(pol)")
# points(Theta.low[nplotb:nplot,npol],Theta_hat[nplotb:nplot,npol],pch=2)
# for (pol in 2:npol){
#   segments(Theta.low[nplotb:nplot,pol-1],Theta_hat[nplotb:nplot,pol-1],
#            Theta.low[nplotb:nplot,pol],Theta_hat[nplotb:nplot,pol])
# }
# abline(v=0)
# points(Theta.low[cbind(1:nsim,policylow)][nplotb:nplot],Theta_hat[cbind(1:nsim,policylow)][nplotb:nplot],pch=4,col=2)

plot(Theta[nplotb:nplot,1],Theta.low[nplotb:nplot,1],
     xlim=range(Theta[nplotb:nplot,]),ylim=range(Theta.low[nplotb:nplot,]),pch=3,
     xlab="$\\Theta(\\texttt{pol})$",ylab="$\\hat{q}_{n,0.05}(\\texttt{pol})$")
points(Theta[nplotb:nplot,npol],Theta.low[nplotb:nplot,npol],pch=2)
for (pol in 2:npol){
  segments(Theta[nplotb:nplot,pol-1],Theta.low[nplotb:nplot,pol-1],
           Theta[nplotb:nplot,pol],Theta.low[nplotb:nplot,pol])
}
abline(h=0)
abline(v=0)
abline(a=0,b=1)
points(Theta[cbind(1:nsim,policylow)][nplotb:nplot],Theta.low[cbind(1:nsim,policylow)][nplotb:nplot],pch=4,col=2)

#The following line needs to be executed to finish the figure file.
#dev.off()

#The following plot is used to ilustrate the procedure in subsection 4.3.
library(ellipse)
#tikz(file = "~/Documents/publications/personalized medecine/beta2_3.tex", width = 5, height = 5)
plot(reglin$coefficients[3],reglin$coefficients[4],pch=4,col=2,
     xlab="$\\beta^\\star_2$",ylab="$\\beta^\\star_3$",xlim=c(0,1.4),ylim=c(0,2.1))
for (lev in c(0.25,0.5,0.75,0.90)){
  lines(ellipse(beta_var[3:4,3:4],level=lev,centre=reglin$coefficients[3:4]))}
for(i in seq(1,Nsample)){if (zscore[i]<thresh[policylow[sim]]){
  curve(-x/X[i],add=TRUE,lty=3)
}else{curve(-x/X[i],add=TRUE)}
}
points(reglin$coefficients[3],reglin$coefficients[4],pch=4,col=2)
points(b2,b3,pch=3,col=2)
#dev.off()
  
#####Numerical quantities####

#Power (size)
mean(linreg<0.05)
mean(linreg<0.05 & minz<qnorm(0.05))
apply(Theta.low>0,2,mean)
mean(Theta.low[cbind(1:nsim,policylow)]>0)
mean(Theta.low[cbind(1:nsim,policylow)]>0 & linreg<0.05 & linreg<0.05)

#Coverage
apply(Theta>Theta.low,2,mean)#underestimation
apply(Theta<Theta.low,2,mean)#overestimation
mean(Theta[cbind(1:nsim,policylow)]>Theta.low[cbind(1:nsim,policylow)])#underestimation
mean(Theta[cbind(1:nsim,policylow)]<Theta.low[cbind(1:nsim,policylow)])#overestimation
apply(Theta>Theta.up,2,mean)#underestimation
apply(Theta<Theta.up,2,mean)#overestimation
mean(Theta[cbind(1:nsim,policylow)]>Theta.up[cbind(1:nsim,policylow)])#underestimation
mean(Theta[cbind(1:nsim,policylow)]<Theta.up[cbind(1:nsim,policylow)])#overestimation
apply(Theta>Theta_hat,2,mean)#underestimation
apply(Theta<Theta_hat,2,mean)#overestimation
mean(Theta[cbind(1:nsim,policylow)]>Theta_hat[cbind(1:nsim,policylow)])#underestimation
mean(Theta[cbind(1:nsim,policylow)]<Theta_hat[cbind(1:nsim,policylow)])#overestimation

#False positive
apply(Theta[linreg<0.05 & minz<qnorm(0.05),]<0,2,mean)
sapply(1:npol,function(pol) mean(Theta[Theta.low[,pol]>0,pol]<0))
mean(Theta[cbind(1:nsim,policylow)][Theta.low[cbind(1:nsim,policylow)]>0]<0)

Theta.low[which(Theta[cbind(1:nsim,policylow)]<0 & Theta.low[cbind(1:nsim,policylow)]>0),]
Theta[which(Theta[cbind(1:nsim,policylow)]<0 & Theta.low[cbind(1:nsim,policylow)]>0),]

#Personalized proportion
apply(Pneg,2,mean)
mean(Pneg[cbind(1:nsim,policylow)])
mean(Pneg[cbind(1:nsim,policylow)][Theta.low[cbind(1:nsim,policylow)]>0])

#Choice of threshold
hist(policylow,breaks=16)
hist(policylow[Theta[cbind(1:nsim,policylow)]>0],breaks=16)
table(policylow[Theta.low[cbind(1:nsim,policylow)]>0])/sum(Theta.low[cbind(1:nsim,policylow)]>0)
min(Cor_DT[,9])
