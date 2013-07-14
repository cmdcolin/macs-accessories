library(DIME)
#slide1<-slideMean(x=normdiff1,windowsize=100,slide=100)
#slide2<-slideMean(x=normdiff2,windowsize=100,slide=100)
nddime<-slideMean(x=normdiff1-normdiff2,windowsize=100,slide=100)



plot((1:length(nddime)),nddime,col=bestClass+1,pch='.',cex=bestClass*2+1, xlab="Position",ylab="Difference")
legend("bottomright",legend=c("Conserved","Differential"),fill=c("red","black"))

DIME.plot.fit(nddime,retdime)

plot((1:length(nddime)),nddime,col=bestClass+1,pch='.',cex=bestClass*2+1, xlab="Position",ylab="Difference")
legend("bottomright",legend=c("Conserved","Differential"),fill=c("red","black"))
title('DIME GNG Model K=6')



system.time(retdime2<-DIME(nddime,gng.K=2))

tr<-retdime2
bestClass = tr$best$class 
plot((1:length(nddime)),nddime,col=bestClass+1,pch='.',cex=bestClass*2+1, xlab="Position",ylab="Difference")
legend("bottomright",legend=c("Conserved","Differential"),fill=c("black","red"))
title('DIME GNG Model K=2')





system.time(retdime6<-DIME(nddime,gng.K=6))

tr<-retdime6
bestClass = tr$best$class 
plot((1:length(nddime)),nddime,col=bestClass+1,pch='.',cex=bestClass*2+1, xlab="Position",ylab="Difference")
legend("bottomright",legend=c("Conserved","Differential"),fill=c("black","red"))
title('DIME GNG Model K=6')


#######
# Not working
testNDs<-log(wiggles$treat1/wiggles$treat2)

system.time(retlog<-DIME(slideMean(testNDs,100,100)))
DIME.plot.fit(retlog)