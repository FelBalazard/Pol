#####Evaluation####
library(MASS)
b1=1
Nsample=300
beta2=c(0,0.1,0.2,0.3,0.4,0.5)

result=numeric(length(beta2))
resultz=numeric(length(beta2))
set.seed(12)
X=2*runif(Nsample,0,1)-1
Tr=rep(c(0,1),Nsample/2)
sig=1
nsim=4000
thresh=qnorm(c(0.5,0.4,0.3,0.25,0.2,0.15,0.1,0.07,0.05,0.04,0.025,0.02))
npol=length(thresh)

for (i in 1:length(beta2)){
b2=beta2[i]

Theta_hat=matrix(numeric(npol*nsim),nsim,npol)
Theta.low=matrix(numeric(npol*nsim),nsim,npol)
Theta.up=matrix(numeric(npol*nsim),nsim,npol)
Pneg=matrix(numeric(npol*nsim),nsim,npol)
sgnb2=numeric(nsim)
linreg=numeric(nsim)
minz=numeric(nsim)
maxz=numeric(nsim)

b3=-b2/min(X)
#X=rexp(Nsample,1)
#X=X-mean(X)

for(sim in 1:nsim){
  set.seed(13+sim)
  eps=rnorm(Nsample,0,sig)
  Y=b1*X+b2*Tr+b3*X*Tr+eps
  reglin=lm(Y~X*Tr)
  beta_var=vcov(reglin)
  #sgn=sign(reglin$coefficients[3])
  NMCMC=4000
  betaci=mvrnorm(NMCMC,reglin$coefficients[3:4],beta_var[3:4,3:4])
  deltaci=(betaci[,2]%*%t(X)+betaci[,1])#*sgn
  zscore=(reglin$coefficients[4]*X+reglin$coefficients[3])/sqrt(beta_var[3,3]+2*beta_var[4,3]*X+X^2*beta_var[4,4])#*sgn
  #hist(zscore,breaks=100)
  policy=sapply(1:npol,function(i) zscore<thresh[i])
  Thetaci=-1/Nsample*t(sapply(1:NMCMC,function(i) deltaci[i,]%*%policy))
  Theta_hat[sim,]=apply(Thetaci,2,mean)
  Theta.low[sim,]=sapply(1:npol,function(pol) quantile(Thetaci[,pol],0.05))
  Theta.up[sim,]=sapply(1:npol,function(pol) quantile(Thetaci[,pol],0.10))
  Pneg[sim,]=apply(policy==T,2,mean)
  linreg[sim]=summary(lm(Y~X*Tr))$coefficient[4,4]
  minz[sim]=min(zscore)
  sgnb2[sim]=reglin$coefficients[3]
}
policylow=sapply(1:nsim,function(i) min(which(Theta.low[i,]==max(Theta.low[i,]))))

result[i]=mean(Theta.low[cbind(1:nsim,policylow)]>0 & linreg<0.05)
resultz[i]=mean(minz<qnorm(0.05) & linreg<0.05)
}

#####Plots####
#Plot null hypothesis
library(tikzDevice)
#tikz(file = "~/Documents/publications/personalized medecine/null2.tex", width = 4, height = 5)
plot(beta2,result,type='l',xlab="$\\beta_2$",ylab="Type I error",ylim=c(0,0.055))
abline(h=0.05,lty=2)
#dev.off()
  