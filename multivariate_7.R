#IST####
load("~/Documents/Médecine personnalisée/IST/IST.RData")
library(randomForest)
library(missMDA)
library(MASS)
library(e1071)#sigmoid function

#Data preparation.
data=ist[,c(2:21,23)]#variable selection as described in the appendix.
data$HOURLOCAL[data$HOURLOCAL==99]=NA#missing values are coded by 99 for this variable.
set.seed(12)#Best seed ever.
fulldata=imputeFAMD(data)#Missing data imputation for covariates.
rct=fulldata$completeObs
rct$group=ist$aspirin
rct=rct[!is.na(ist$deathdep),]#Exclusion for missing outcome.
Nsample=dim(rct)[1]
outcome=1-ist$deathdep[!is.na(ist$deathdep)]#Reversal of outcome.

#Prediction model. keep.inbag allows to know which observation are used for each tree.
a=randomForest(rct,factor(outcome),keep.inbag = T,ntree = 1500)

#Feature extraction
rct_twins=rct
rct_twins$group=1-rct_twins$group# Same coavriates but opposite treatment.
b=predict(a,rct_twins,predict.all=T,type="vote",norm.votes = T)
vote=matrix(as.numeric(b$individual),Nsample,500)#Matrix of votes tree times participants.
OOB=(a$inbag==0)#Out-Of-Bag observations
pred=apply(OOB*vote,1,sum)/apply(OOB,1,sum)#Out-of-bag predictions for twins
delta=(a$votes[,2]-pred)*(rct$group==1)+(pred-a$votes[,2])*(rct$group==0)
rct_twins$group=0#To have the prognostic term.
pronostic=predict(a,rct_twins,predict.all=T,type="vote",norm.votes = T)
vote_prono=as.numeric(pronostic$individual)
Z1=apply(OOB*vote_prono,1,sum)/apply(OOB,1,sum)#Out-Of-Bag prognostic prediction
Z3=(delta-mean(delta))#Out-Of-Bag treatment interaction
Z3group=Z3*rct$group#Interaction with treatment

#The unidimensional procedure
reglin=glm(outcome~Z1+rct$group+Z3group,family="binomial")
beta_var=vcov(reglin)
NMCMC=5000
betaci=mvrnorm(NMCMC,reglin$coefficients,beta_var)
deltaci=sigmoid(betaci[,2]%*%t(Z1)+betaci[,1]+betaci[,4]%*%t(Z3)+betaci[,3])-sigmoid(betaci[,2]%*%t(Z1)+betaci[,1])
probdelta=apply(deltaci>0,2,mean)
thresh=c(0.04,0.05,0.06,0.065,0.07,0.075,0.08)
npol=length(thresh)
policy=sapply(1:npol,function(i) probdelta<thresh[i])
Thetaci=-1/Nsample*t(sapply(1:NMCMC,function(i) deltaci[i,]%*%policy))
Theta_hat=apply(Thetaci,2,mean)
Theta.low=sapply(1:npol,function(pol) quantile(Thetaci[,pol],0.05))
Theta.up=sapply(1:npol,function(pol) quantile(Thetaci[,pol],0.95))
Cor_DT=sapply(1:npol,function(pol) if(sum(policy[,pol])==0){1}else{
  cor(-deltaci[,min(which(Z3==sign(reglin$coefficients[4])*max(sign(reglin$coefficients[4])*Z3[policy[,pol]])))]
      ,Thetaci[,pol])})
Pneg=apply(policy==T,2,mean)

#Results
plot(Theta.low,Theta_hat,pch=3)
policylow=min(which(Theta.low==max(Theta.low)))
Cor_DT
Pneg[policylow]
Theta.low[policylow]
Theta_hat[policylow]
mean(Thetaci[,policylow]<0)
policylow

