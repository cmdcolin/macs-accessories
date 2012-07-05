source('readplot.R')
source('knitr.R')
debug=TRUE
wig1=WiggleClass('S96')
wig2=WiggleClass('S96rep2u')
#wig1$loadWiggles() 
#wig2$loadWiggles()
wig1$loadWiggles(globalenv()) 
wig2$loadWiggles(globalenv())
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

r1=plotMaxAvgZscore('Max Avg  S96rep1 peak NormDiff score vs S96rep2 synteny w=100', wig1, wig2, wz1, wz2,'orange', 'blue')
r2=plotMaxAvgZscore('Max Avg S96rep2 peak NormDiff score vs S96rep1 synteny w=100', wig2, wig1, wz4, wz3, 'orange','darkviolet')

## Test  warning signals
#withCallingHandlers(plotMaxAvgZscore('Max Avg HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4, wz3, 'orange','darkviolet') , warning=function(c) recover())

#r2=plotMaxAvgZscore('Max Avg HS959 peak NormDiff score vs S96 synteny w=100', wig2, wig1, wz4[[5]], wz3[[5]], 'orange','darkviolet')
#i=0
#apply(cbind(wz3,wz4),1,function(x)
#  {
#  i=i+1;
#  cat('step', i);
##  plotMaxAvgZscore(paste('str',i), wig2, wig1, x[2], x[1], 'orange','darkviolet')
#  })
#######



##
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





plotBarCharts(wig1,wig2,wz1,wz2)

unique=uniqueBed(wig1$peaks,wig2$peaks)
shared=intersectBed(wig1$peaks,wig2$peaks)
uniqueid=match(unique$V4,wig1$peaks$V4)
sharedid=match(shared$V4,wig1$peaks$V4)  
maxw1<-wig1$getMaxAvgZscore(wz1)
maxw2<-wig2$getMaxAvgZscore(wz2)



xx=getBayesian(wig1,wig2,wz1,wz2,maxw1,maxw2,uniqueid)/4
yy=getBayesian(wig1,wig2,wz2,wz1,maxw1,maxw2,sharedid)/4
par(mfrow=c(1,2))
barplot(sort(xx),xlab='conditional probability', ylim=c(0,1))
title('Conditional probability of unique peaks in HS959')
barplot(sort(yy),xlab='conditional probability', ylim=c(0,1))
title('Conditional probability for shared peaks in HS959')

plotBarChartY('Title',wig1,wig2,wz1,wz2)




### Color plots
plotMaxAvgZscoreColorY('S96rep1 peaks vs S96rep2 synteny',wig1,wig2,wz1,wz2)
plotMaxAvgZscoreColorY('S96rep2 peaks vs S96rep1 synteny',wig2,wig1,wz4,wz3)


### Zoom on unique plots
b1=plotMaxAvgZscoreColorUnique('S96 peaks vs S96rep2 synteny (Unique only)',wig1,wig2,wz1,wz2)
b2=plotMaxAvgZscoreColorUnique('HS959 peaksvs S96 synteny (Unique only)',wig2,wig1,wz4,wz3)


ret=plotMaxAvgZscoreColorUnique('S96 peaks vs HS959 synteny (Unique only)',wig1,wig2,wz1,wz2)
print(ret)
b1=as.numeric(ret[,1])


bedselect=wig2$peaks[b1,]
selection1=wig1$Z(bedselect)
selection2=wig2$Z(bedselect)
reads1=wig1$getChipReads(bedselect)
reads2=wig2$getChipReads(bedselect)

plotOverlaps(b1,wig1,wig2,1,250)
plotOverlaps(b1,wig1,wig2,2,250)
plotOverlaps(b1,wig1,wig2,3,250)
plotOverlaps(b1,wig1,wig2,4,250)
plotOverlaps(b1,wig1,wig2,5,250)

## Plot - axis 1
xx=seq(bedselect[1,2]-10,bedselect[1,3]+10,by=wig2$spacing)
plot(xx,unlist(selection1[1]),type='l',xlab='Genome position',ylab='NormDiff score')
lines(xx,unlist(selection2[1]))
polygon(c(xx,rev(xx)),c(unlist(selection1[1]),numeric(length(unlist(selection1[1])))),col='#0000bb44')
polygon(c(xx,rev(xx)),c(unlist(selection2[1]),numeric(length(unlist(selection2[1])))),col='#bb000044')
## Plot - axis 2
par(new=TRUE)
xx=seq(bedselect[1,2],bedselect[1,3],by=wig2$spacing)
plot(xx,unlist(reads1[1]),type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
lines(xx,unlist(reads2[1]),col='red')
polygon(c(xx,rev(xx)),c(unlist(reads1[1]),numeric(length(unlist(reads1[1])))),col='#0000bb22')
polygon(c(xx,rev(xx)),c(unlist(reads2[1]),numeric(length(unlist(reads2[1])))),col='#bb000088')
axis(4)
mtext('ChIP-seq reads',side=4,line=3)


## Setup data
datasort1=as.numeric(wza1[[1]][,4])
datasort2=as.numeric(wza2[[1]][,4])
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


## Q-Q Plot
clone=datasort
qqnorm(clone)
qqline(clone,col=2)

# densityplot
#library(lattice)
#d<-densityplot(~datasort,data=faithful,kernel="rect")
#d<densityplot(~datasort,kernel="gaussian")
#plot(d, main='Kernel density of NormDiff scores') # plots the results


#### Draw rainbow
plot(c(100, 250), c(300, 450), type = "n",main="r")
i <- 4*(0:10)
## draw rectangles with bottom left (120, 300)+i  and top right (177, 380)+i
rect(120, 300+i, 177, 380+i, col=rainbow(11, start=.0,end=.7))



#####################3
# SCRAP CODE

#plot(p1,type='l',col="#aa0000bb",ylim=c(0,20),xlim=c(0,100000))
#lines(p2,col="#101010cc")


wza1=wig1$Zall()
wza2=wig2$Zall()
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
