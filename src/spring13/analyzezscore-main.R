
#out<-getMatchListMod(chrnames,retable1,retable2)
#outg<-getScoresMod(retable1,retable2)
#plot(outg$col1,outg$col3,pch='.')




wiggles<-getscoresmod(chrnames,macswiggle[[1]]$treat,macswiggle[[1]]$control,macswiggle[[2]]$treat,macswiggle[[2]]$control)


load('wiggledata.rda')
#wiggles=wiggles



m1=median(wiggles$control1/wiggles$treat1)
m2=median(wiggles$control2/wiggles$treat2)

background1<-wiggles$treat1-wiggles$control1
background2<-wiggles$treat2-wiggles$control2

scale1<-wiggles$treat1-wiggles$control1/m1
scale2<-wiggles$treat2-wiggles$control2/m2


treatmeans<-slideMean(wiggles$treat1)
controlmeans<-slideMean(wiggles$control1)
vlist1<-sqrt(treatmeans+controlmeans/m1^2)
treatmeans<-slideMean(wiggles$treat2)
controlmeans<-slideMean(wiggles$control2)
vlist2<-sqrt(treatmeans+controlmeans/m2^2)


v1<-sqrt(mean(wiggles$treat1)+mean(wiggles$control1)/m1^2)
v2<-sqrt(mean(wiggles$treat2)+mean(wiggles$control2)/m2^2)
vlist1f<-sapply(vlist1,function(v,vall){max(v,vall)},v1)
vlist2f<-sapply(vlist2,function(v,vall){max(v,vall)},v2)

normdiff1<-scale1/vlist1f
normdiff2<-scale2/vlist2f

#no background subtraction
plot(wiggles$treat1,wiggles$treat2,pch='.',xlab='S96 rep1 treated',ylab='S96 rep2 treated')
title('S96 replicates ChIP-seq scores (no control)')





# background subtraction
plot(background1,background2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (background subtract)')


#background scaling
plot(scale1,scale2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (background scaling)')


#normalization
plot(l1/v1,l2/v2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=all)')
smoothScatter(l1/v1,l2/v2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=all)')




#normdiff local

plot(normdif1,normdiff2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=local)')

smoothScatter(normdiff1,normdiff2,pch='.',xlab='S96 rep1',ylab='S96 rep2',ylim=c(-5,35),xlim=c(-5,35))
title('S96 replicates ChIP-seq scores (normdiff w=local)')





#no background subtraction-MA plot
A=(wiggles$treat1+wiggles$treat2)*1/2
M=(wiggles$treat1-wiggles$treat2)
plot(A,M,pch='.',xlab='M',ylab='A')
title('S96 MA plot for replicates (raw scores)')




#no background subtraction-MA plot
A=(background1+background2)*1/2
M=(background1-background2)
plot(A,M,pch='.',xlab='M',ylab='A')
title('S96 MA plot for replicates (background subtraction)')



#scaling-MA plot
m1=median(wiggles$control1/wiggles$treat1)
m2=median(wiggles$control2/wiggles$treat2)
l1<-wiggles$treat1-wiggles$control1/m1
l2<-wiggles$treat2-wiggles$control2/m2
A=(l1+l2)*1/2
M=(l1-l2)
plot(A,M,pch='.',xlab='A',ylab='M')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (background scaling) lm')


#normalization MA-plot
v1<-sqrt(mean(wiggles$treat1)+mean(wiggles$control1)/m1^2)
v2<-sqrt(mean(wiggles$treat2)+mean(wiggles$control2)/m2^2)
M<-(l1/v1-l2/v2)
A<-(l1/v1+l2/v2)*1/2
plot(A,M,pch='.',xlab='A',ylab='M')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=all)')





#normalization MA-plot
M<-(l1/vlist1-l2/vlist2)
A<-(l1/vlist1+l2/v2)*1/2
plot(A,M,pch='.',xlab='A',ylab='M')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=local)')
