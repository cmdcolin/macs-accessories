


###############
# HYPOTHESIS TESTING
#############


## Setup data
datasort=as.numeric(wza1[[1]][,4])  



# Basic density plot
plot(datasort,dnorm(datasort))


## Kernel density plot
d <- density(datasort,adjust=1.2) # returns the density data
plot(d, main='Kernel density of NormDiff scores') # plots the results 
polygon(d, col="#BB2222CC", border="#222244") 


## Q-Q Plot whole genome
clone=datasort
qqnorm(clone)
qqline(clone,col=2)




plot(ren,dnorm(ren))


## Q-Q Plot peak means
clone1=sapply(wz1,mean)
clone2=sapply(wz2,mean)
qqnorm(clone)
points(qqnorm(clone2,plot=FALSE),col="blue")
points(qqnorm(clone1,plot=FALSE),col="red")
qqline(rnorm(5000),col=2)



## Q-Q plot peak scores
select=(wig1$peaks[,1]=="chr01.fsa")
datasel1=unlist(wz1[select])
datasel2=unlist(wz2[select])
clone=datasort
qqnorm(clone)
points(qqnorm(datasel1,plot=FALSE),col="blue")
points(qqnorm(datasel2,plot=FALSE),col="red")
qqline(rnorm(5000),col=2)
legend('topleft',legend=c('Whole genome','S96 peaks','HS959 peaks'),fill=c('black','blue','red'))



## Q-Q plot MEAN peak scores
dataselmean=wzamax1[,4]
dataselmean1=maxw1
dataselmean2=maxw2
clone=dataselmean
qqnorm(clone,ylim=c(-3,8))
points(qqnorm(dataselmean1,plot=FALSE),col="blue")
points(qqnorm(dataselmean2,plot=FALSE),col="red")
qqline(clone,col=2)
legend('topleft',legend=c('Whole genome','S96 peaks','HS959 syntenic'),fill=c('black','blue','red'))


## Kernel density plot MEANS
d <- density(dataselmean,adjust=1.4) # returns the density data
plot(d, main='Kernel density of MEAN NormDiff scores') # plots the results 
polygon(d, col="#2222BBCC", border="#222244")
d1 <- density(dataselmean1,adjust=1.4) # returns the density data
d2 <- density(dataselmean2,adjust=1.4) # returns the density data
lines(d1)
lines(d2)
polygon(d1, col="#BB222288", border="#222244")
polygon(d2, col="#22BB2288", border="#222244")
legend('topright',legend=c('Whole genome','S96 peaks','HS959 syntenic'),fill=c('blue','green','red'))









# Basic density plot
# note: needs setup
plot(datasort,dnorm(datasort))
points(datasel1,dnorm(datasel1),col="blue")
points(datasel2,dnorm(datasel2),col="red")


plotMAPlot2('MA plotof NormDiff scores',wig1,wig2,wz1,wz2,wz3,wz4,0)


## Kernel density plot
d <- density(datasort,adjust=1.4) # returns the density data
plot(d, main='Kernel density of NormDiff scores') # plots the results 
polygon(d, col="#2222BBCC", border="#222244")
d1 <- density(datasel1,adjust=1.4) # returns the density data
d2 <- density(datasel2,adjust=1.4) # returns the density data
lines(d1)
lines(d2)
polygon(d1, col="#BB222288", border="#222244")
polygon(d2, col="#22BB2288", border="#222244")
legend('topright',legend=c('Whole genome','S96 peaks','HS959 peaks'),fill=c('blue','green','red'))




# Change colors
ret=plotMaxAvgZscore('Max Avg NormDiff S96 peaks vs HS959 synteny',wig1,wig2,wz1,wz2,'blue','red')
ret=plotMaxAvgZscore('Max Avg NormDiff HS959 peaks vs S96 synteny',wig2,wig1,wz4,wz3,'blue','brown')
ret=plotMaxAvgZscore('S96rep1 peaks vs S96rep2 synteny',wig3,wig4,wz5,wz6,'green','orange')
ret=plotMaxAvgZscore('S96rep2 peaks vs S96rep1 synteny',wig4,wig3,wz7,wz8,'green','orange',my.pch=0)

ret=plotSortedMaxAvgZscoreX('Sorted HS959 NormDiff Scores adjascent to S96 peaks',wig1,wig2,wz1,wz2,'blue','red')

ret=plotSortedMaxAvgZscoreX('Sorted S96 NormDiff Scores adjascent to HS959 peaks',wig2,wig1,wz4,wz3,'blue','brown')


# MACS pvalues
ret=plotMaxAvgZscoreU('S96 peaks vs HS959 synteny',wig1,wig2,wz1,wz2,0.5)

# NormDiff pvalues
ret=plotMaxAvgZscoreW('S96 peaks vs HS959 synteny',wig1,wig2,wz1,wz2,0.05)

# MA plot
ret=plotMaxAvgZscoreR('S96 peaks vs HS959 synteny',wig1,wig2,wz1,wz2,wz3,wz4,0.05)

# Get peaks cutoff
plotZscoreCutoff('S96 peaks vs HS959 synteny (New peaks)',wig2,wig1,wz4,wz3,0.01)
plotZscoreCutoff('HS959 peaks vs S96 synteny (New peaks)',wig1,wig2,wz1,wz2,0.01)
plotZscoreCutoff('HS959 peaks vs S96 synteny (New peaks)',wig1,wig2,wz1,wz2,0.05)
plotZscoreCutoff('S96 peaks vs HS959 synteny (New peaks)',wig2,wig1,wz4,wz3,0.05)

ret=plotZscoreCutoff('S96 peaks vs HS959 synteny (New peaks)',wig3,wig4,wz5,wz6,0.01)
ret=plotZscoreCutoff('HS959 peaks vs S96 synteny (New peaks)',wig4,wig3,wz8,wz7,0.01)


plotZscoreCutoffShared('test',wig2,wig1,wz4,wz3,0.05)
# Two color plot
ret=plotZscoreCutoff2('HS959 peaks conserved in S96',wig1,wig2,wz1,wz2,0.05,0.01,'blue','red')
ret=plotZscoreCutoff2('S96 peaks conserved in HS959',wig2,wig1,wz4,wz3,0.05,0.01,'blue','red')
ret=plotZscoreCutoff2('Max Avg NormDiff S96rep1 peaks vs S96rep2 synteny',wig3,wig4,wz5,wz6,0.05,0.01,'blue','red')
ret=plotZscoreCutoff2('Max Avg NormDiff S96rep2 peaks vs S96rep1 synteny',wig4,wig3,wz8,wz7,0.05,0.01,'blue','brown')

indexes=which(ret==TRUE)
wig3$peaks[ret,]




#Mod fix
ret=plotZscoreCutoffShared('S96 peaks vs HS959 synteny (New peaks)',wig1,wig2,wz1,wz2,0.01)

# select peaks cutoff
wig1$peaks[ret,]





vall1=estimateVarianceAllMod(wig1,wig2)
Zmod(wig1$peaks, wig1,wig2,Zaddxi,vall1)


#### Scrap code

# Histogram
#h=hist(datasort,breaks=50)
#xfit=datasort
#yfit=dnorm(datasort)
#yfit <- yfit*diff(h$mids[1:2])*length(datasort)
#lines(xfit, yfit, col="blue", lwd=1)

#qqline(rnorm(5000),col=2)
#par(new=TRUE)
#clone=sapply(wz1,mean)
#qqnorm(clone)
