


###############
# HYPOTHESIS TESTING
#############


## Setup data
datasort=as.numeric(wza1[[1]][,4])
plot(datasort,dnorm(datasort))

### Histogram
h=hist(datasort,breaks=50)
xfit=datasort
yfit=dnorm(datasort)
yfit <- yfit*diff(h$mids[1:2])*length(datasort)
lines(xfit, yfit, col="blue", lwd=1) 

## Kernel density plot
d <- density(datasort,adjust=1.2) # returns the density data
plot(d, main='Kernel density of NormDiff scores') # plots the results 
polygon(d, col="#BB2222CC", border="#222244") 


## Q-Q Plot wholte genome
datasort=as.numeric(wza1[[1]][,4])
clone=datasort
qqnorm(clone)
qqline(clone,col=2)

## Q-Q Plot peaks

clone1=sapply(wz1,mean)
clone2=sapply(wz2,mean)

qqnorm(clone)
points(qqnorm(clone2,plot=FALSE),col="blue")
points(qqnorm(clone1,plot=FALSE),col="red")
qqline(rnorm(5000),col=2)


qqline(rnorm(5000),col=2)
par(new=TRUE)
clone=sapply(wz1,mean)
qqnorm(clone)