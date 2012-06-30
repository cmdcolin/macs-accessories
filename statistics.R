source('readplot.R')
source('knitr.R')
debug=TRUE
wig1=WiggleClass('S96')
wig2=WiggleClass('HS959')
wig1$loadWiggles() 
wig2$loadWiggles()
#wig1$loadWiggles(globalenv()) 
#wig2$loadWiggles(globalenv())
###
wig1$estimateScalingFactor()
wig1$estimateVarianceAll()
wig2$estimateScalingFactor()
wig2$estimateVarianceAll()


############
# Z score peaks
wz1=wig1$Z(wig1$peaks)
wz2=wig2$Z(wig1$peaks)
wz4=wig2$Z(wig2$peaks)
wz3=wig1$Z(wig2$peaks)

r1=plotMaxAvgZscore('Max Avg  S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2,'orange', 'blue')
r2=plotMaxAvgZscore('Max Avg HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'orange','darkviolet')



r3=plotAvgZscore('Mean S96 peak NormDiff score vs HS959 synteny w=100', wig1, wig2, wz1, wz2, 'lightblue', 'orange')
r4=plotAvgZscore('Mean HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'pink','orange')

plotSortedMaxAvgZscoreX('Sorted HS959 Avg Normdiff in S96 peak regions',wig1,wig2,r3,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Avg Normdiff in HS959 peak regions',wig2,wig1,r4,'green','blue')
plotSortedMaxAvgZscoreX('Sorted HS959 Max Avg Normdiff in S96 peak regions',wig1,wig2,r1,'red','blue')
plotSortedMaxAvgZscoreX('Sorted S96 Max Avg Normdiff in HS959 peak regions',wig2,wig1,r2,'green','blue') 





plotSortedMaxAvgZscore('Sorted HS959 Max Avg Normdiff connected with S96 peak regions',wig1,wig2,r1,'#00334407','#55221133') 
plotSortedMaxAvgZscore('Sorted S96 Max Avg Normdiff connected with HS959 peak regions',wig2,wig1,r2,'#00661107','#55221133')
plotSortedMaxAvgZscore('Sorted HS959 Max Avg Normdiff connected with S96 peak regions',wig1,wig2,r1,'#00334407','#55221133','#003344','#552211') 
plotSortedMaxAvgZscore('Sorted S96 Max Avg Normdiff connected with HS959 peak regions',wig2,wig1,r2,'#00661107','#55221133','#006611','#552211')







unique=uniqueBed(wig1$peaks,wig2$peaks)
shared=intersectBed(wig1$peaks,wig2$peaks)
uniqueid=match(unique$V4,wig1$peaks$V4)
sharedid=match(shared$V4,wig1$peaks$V4)  
maxw1<-wig1$getMaxAvgZscore(wz1)
maxw2<-wig2$getMaxAvgZscore(wz2)



xx=getBayesian(wig1,wig2,wz1,wz2,maxw1,maxw2,uniqueid)
yy=getBayesian(wig1,wig2,wz2,wz1,maxw1,maxw2,sharedid)
par(mfrow=c(1,2))
barplot(sort(xx),xlab='conditional probability', ylim=c(0,1))
title('Conditional probability of unique peaks in HS959')
barplot(sort(yy),xlab='conditional probability', ylim=c(0,1))
title('Conditional probability for shared peaks in HS959')



plotMaxAvgZscoreColor('Colors1',wig1,wig2,wz1,wz2)
plotMaxAvgZscoreColor('Colors2',wig2,wig1,wz4,wz3)




datasort=as.numeric(wza1[[1]][,4])
plot(datasort,dnorm(datasort))

h=hist(datasort,breaks=50)
xfit=datasort
yfit=dnorm(datasort)
yfit <- yfit*diff(h$mids[1:2])*length(datasort)
lines(xfit, yfit, col="blue", lwd=1) 


d <- density(datasort,adjust=1.5) # returns the density data
plot(d, main='Kernel density of NormDiff scores') # plots the results 
polygon(d, col="#BB2222CC", border="#222244") 

clone=datasort
qqnorm(clone)
qqline(clone,col=2)
qqplot(clone, datasort)



x <- rnorm(100, mean=5, sd=2)

# sort in ascending order
x.sorted <- sort(x)


#plot(p1,type='l',col="#aa0000bb",ylim=c(0,20),xlim=c(0,100000))
#lines(p2,col="#101010cc")


#wza1=wig1$Zall()
#wza2=wig2$Zall()
#wzamax1=wig1$getMaxAvgZscoreAll(wza1)
#wzamax2=wig2$getMaxAvgZscoreAll(wza2)
#c1=plotZall(wzamax1,wig1,wig2)
#c2=plotZall(wzamax2,wig2,wig1)
#plotZall(wzamax1,wig1,wig2,c1)
#plotZall(wzamax2,wig2,wig1,c2)

#plotZall2(c1,wig1,wig2)
#plotZall3(c1,wig1,wig2)

#plotMaxAvgReads('S96 avg peak reads vs HS959 syntenic', wig1, wig2, 'green', 'blue') 
#plotMaxAvgReads('HS959 avg peak reads vs S96 syntenic', wig2, wig1, 'green', 'red')
