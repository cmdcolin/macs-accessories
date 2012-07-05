###############
# Bayes and correlation scores
##########
xx=getBayesian(wig1,wig2,wz1,wz2,maxw1,maxw2,uniqueid)/4
yy=getBayesian(wig1,wig2,wz2,wz1,maxw1,maxw2,sharedid)/4
par(mfrow=c(1,2))
barplot(sort(xx),xlab='conditional probability', ylim=c(0,1))
title('Conditional probability of unique peaks in HS959')
barplot(sort(yy),xlab='conditional probability', ylim=c(0,1))
title('Conditional probability for shared peaks in HS959')
plotBarChartY('Title',wig1,wig2,wz1,wz2)


### Color plots
plotMaxAvgZscoreColorY('S96 peaks vs HS959 synteny',wig1,wig2,wz1,wz2)
plotMaxAvgZscoreColorY('HS959 peaksvs S96 synteny',wig2,wig1,wz4,wz3)
### Zoom on unique plots
b1=plotMaxAvgZscoreColorUnique('S96 peaks vs HS959 synteny (Unique only)',wig1,wig2,wz1,wz2)
b2=plotMaxAvgZscoreColorUnique('HS959 peaksvs S96 synteny (Unique only)',wig2,wig1,wz4,wz3)

#### Plot Overlaps
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
