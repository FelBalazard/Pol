#This script is used to plot figure 2. It needs the environment to contain a reglin from pol1_3.1.R.
x=(0:100)/200-1.1
z=(reglin$coefficients[3]+x*reglin$coefficients[4])/sqrt(beta_var[3,3]+2*beta_var[3,4]*x+x^2*beta_var[4,4])
library(tikzDevice)#Allows to have beautiful latex type setting in axis labels. <3
#tikz(file = "~/Documents/publications/personalized medecine/zdelta.tex", width = 5, height = 5)
plot(x,z,pch=".",ylab="$z_\\Delta$",xlab="x",type="l",frame.plot = F,axes=F)
abline(a=0,b=0)
xalpha=x[max(which(z<qnorm(0.1)))]
segments(xalpha,z[max(which(z<qnorm(0.1)))],xalpha,0,lty=2)
axis(3,pos=0,at=c(-1,xalpha,-reglin$coefficients[3]/reglin$coefficients[4]), labels=c("$x_0$","$X_\\alpha$","$-\\frac{\\hat\\beta_2}{\\hat\\beta_3}$"))
segments(xalpha,z[max(which(z<qnorm(0.1)))],-100,z[max(which(z<qnorm(0.1)))],lty=2)
axis(2, at=c(z[max(which(z<qnorm(0.1)))],-2,-1,0,1), labels=c("$q_\\alpha$","-2","-1","0","1"))

#The following line needs to be executed to finish the figure file.
#dev.off()
