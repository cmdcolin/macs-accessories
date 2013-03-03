
bed1<-read.table('s96rep1-high_peaks.bed')
bed2<-read.table('s96rep2-high_peaks.bed')
names(bed1)<-c('chr','start','end','name','num')
names(bed2)<-c('chr','start','end','name','num')
bed1$start=as.numeric(bed1$start)
bed1$end=as.numeric(bed1$end)



#############


match=getMatchList(chrnames, macswiggle[[1]]$treat,macswiggle[[2]]$treat)
retable=getScores(match, macswiggle[[1]]$treat,macswiggle[[2]]$treat)
plot(retable$col1,retable$col2,pch='.')
out1<-getPeakScores(bed1,retable)
out2<-getPeakScores(bed2,retable)
points(out1$col1,out1$col2,pch='.',col=2)
points(out2$col1,out2$col2,pch='.',col=3)


###################

cmatch1=getMatchList(chrnames,macswiggle[[1]]$treat,macswiggle[[1]]$control)
cmatch2=getMatchList(chrnames,macswiggle[[2]]$treat,macswiggle[[2]]$control)
retable1=getScores(cmatch1, macswiggle[[1]]$treat,macswiggle[[1]]$control)
retable2=getScores(cmatch2, macswiggle[[2]]$treat,macswiggle[[2]]$control)

#out<-getMatchListMod(chrnames,retable1,retable2)
#outg<-getScoresMod(retable1,retable2)
#plot(outg$col1,outg$col3,pch='.')




rettab<-getscoresmod(chrnames,macswiggle[[1]]$treat,macswiggle[[1]]$control,macswiggle[[2]]$treat,macswiggle[[2]]$control)



#no background subtraction
plot(rettab$treat1,rettab$treat2,pch='.',xlab='S96 rep1 treated',ylab='S96 rep2 treated')
title('S96 replicates ChIP-seq scores (no control)')





# background subtraction
plot(rettab$treat1-rettab$control1,rettab$treat2-rettab$control2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (background subtract)')


#background scaling
m1=median(rettab$control1/rettab$treat1)
m2=median(rettab$control2/rettab$treat2)
l1<-rettab$treat1-rettab$control1/m1
l2<-rettab$treat2-rettab$control2/m2
plot(l1,l2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (background scaling)')


#normalization
v1<-sqrt(mean(rettab$treat1)+mean(rettab$control1)/m1^2)
v2<-sqrt(mean(rettab$treat2)+mean(rettab$control2)/m2^2)
plot(l1/v1,l2/v2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=all)')
smoothScatter(l1/v1,l2/v2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=all)')


treatmeans<-slideMean(rettab$treat1)
controlmeans<-slideMean(rettab$control1)
vlist1<-sqrt(treatmeans+controlmeans/m1^2)
treatmeans<-slideMean(rettab$treat2)
controlmeans<-slideMean(rettab$control2)
vlist2<-sqrt(treatmeans+controlmeans/m2^2)


#normdiff local
vlist1f<-sapply(vlist1,function(v,vall){max(v,vall)},v1)
vlist2f<-sapply(vlist2,function(v,vall){max(v,vall)},v2)

plot(l1/vlist1,l2/vlist2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=local)')

smoothScatter(l1/vlist1,l2/vlist2,pch='.',xlab='S96 rep1',ylab='S96 rep2',ylim=c(-5,35),xlim=c(-5,35))
title('S96 replicates ChIP-seq scores (normdiff w=local)')





#no background subtraction-MA plot
A=(rettab$treat1+rettab$treat2)*1/2
M=(rettab$treat1-rettab$treat2)
plot(A,M,pch='.',xlab='M',ylab='A')
title('S96 MA plot for replicates (raw scores)')




#no background subtraction-MA plot
l1<-rettab$treat1-rettab$control1
l2<-rettab$treat2-rettab$control2
A=(l1+l2)*1/2
M=(l1-l2)
plot(A,M,pch='.',xlab='M',ylab='A')
title('S96 MA plot for replicates (background subtraction)')



#scaling-MA plot
m1=median(rettab$control1/rettab$treat1)
m2=median(rettab$control2/rettab$treat2)
l1<-rettab$treat1-rettab$control1/m1
l2<-rettab$treat2-rettab$control2/m2
A=(l1+l2)*1/2
M=(l1-l2)
plot(A,M,pch='.',xlab='M',ylab='A')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (background scaling) lm')


#normalization MA-plot
v1<-sqrt(mean(rettab$treat1)+mean(rettab$control1)/m1^2)
v2<-sqrt(mean(rettab$treat2)+mean(rettab$control2)/m2^2)
M<-(l1/v1-l2/v2)
A<-(l1/v1+l2/v2)*1/2
plot(A,M,pch='.',xlab='S96 rep1',ylab='S96 rep2')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=all)')





#normalization MA-plot
M<-(l1/vlist1-l2/vlist2)
A<-(l1/vlist1+l2/v2)*1/2
plot(A,M,pch='.',xlab='S96 rep1',ylab='S96 rep2')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=local)')


M<-(l1/vlist1-l2/vlist2)
A<-(l1/vlist1+l2/v2)*1/2
plot(A,M,pch='.',xlab='S96 rep1',ylab='S96 rep2')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=local)')