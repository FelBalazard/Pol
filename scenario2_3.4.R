#####Evaluation####
library(MASS)#multivariate normal

#Coefficients for linear model Y =b1*X+b2*T+b3*X*T+eps
b1=1
b2=0
b3=.4
sig=1#variance of eps

#Sample size
Nsample=500

#Number of simulations
nsim=100

#Thresholds 
thresh=0:20/10-0.4+qnorm(0.05)
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
  NMCMC=100000#number of sampling of beta*
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

#####Numerical quantities####

#Power (size)
mean(linreg<0.05)
mean(linreg<0.05 & minz<qnorm(0.05))
apply(Theta.low>0,2,mean)
mean(Theta.low[cbind(1:nsim,policylow)]>0)
mean(Theta.low[cbind(1:nsim,policylow)]>0 & linreg<0.05)

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
mean(Pneg[cbind(1:nsim,policylow)][Theta.low[cbind(1:nsim,policylow)]>0 & linreg<0.05])
mean(Pneg[cbind(1:nsim,policylow)][policylow>5 & linreg<0.05])
mean(Pneg[,5][policylow>5 & linreg<0.05])

#Choice of threshold
hist(policylow,breaks=16)
hist(policylow[Theta[cbind(1:nsim,policylow)]>0],breaks=16)
table(policylow[Theta.low[cbind(1:nsim,policylow)]>0])/sum(Theta.low[cbind(1:nsim,policylow)]>0)
mean(Cor_DT[,5][Theta.low[cbind(1:nsim,policylow)]>0 & linreg<0.05])
mean(Cor_DT[,5][policylow==5 & linreg<0.05])
mean(Cor_DT[,5][policylow>5 & linreg<0.05])

