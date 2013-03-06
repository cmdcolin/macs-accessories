
#out<-getMatchListMod(chrnames,retable1,retable2)
#outg<-getScoresMod(retable1,retable2)
#plot(outg$col1,outg$col3,pch='.')




wiggles<-getscoresmod(chrnames,macswiggle[[3]]$treat,macswiggle[[3]]$control,macswiggle[[4]]$treat,macswiggle[[4]]$control)

wiggles<-getscoresmodappend(wiggles,chrnames,macswiggle[[1]]$treat,macswiggle[[1]]$control,3)

wiggles<-getscoresmodappend(wiggles,chrnames,macswiggle[[2]]$treat,macswiggle[[2]]$control,4)

head(wiggles)
load('wiggledata.rda')
#wiggles=wiggles

control1<-wiggles$control1
treat1<-wiggles$treat1
control2<-wiggles$control2
treat2<-wiggles$treat2


control1<-wiggles$control3
treat1<-wiggles$treat3
control2<-wiggles$control4
treat2<-wiggles$treat4


m1=median(control1/treat1)
m2=median(control2/treat2)

background1<-treat1-control1
background2<-treat2-control2

scale1<-treat1-control1/m1
scale2<-treat2-control2/m2


treatmeans1<-slideMean(treat1)
controlmeans1<-slideMean(control1)
treatmeans2<-slideMean(treat2)
controlmeans2<-slideMean(control2)
vlist1<-sqrt(treatmeans1+controlmeans1/m1^2)
vlist2<-sqrt(treatmeans2+controlmeans2/m2^2)


v1<-sqrt(mean(treat1)+mean(control1)/m1^2)
v2<-sqrt(mean(treat2)+mean(control2)/m2^2)
vlist1f<-sapply(vlist1,function(v,vall){max(v,vall)},v1)
vlist2f<-sapply(vlist2,function(v,vall){max(v,vall)},v2)

gnormdiff1<-scale1/v1
gnormdiff2<-scale2/v2

normdiff1<-scale1/vlist1f
normdiff2<-scale2/vlist2f



#no background subtraction
plot(treat1,treat2,pch='.',xlab='S96 rep1 treated',ylab='S96 rep2 treated')
title('S96 replicates ChIP-seq scores (no control)')





# background subtraction
plot(background1,background2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (background subtract)')


#background scaling
plot(scale1,scale2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (background scaling)')


#normalization
plot(gnormdiff1,gnormdiff2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=all)')




#normdiff local

plot(normdiff1,normdiff2,pch='.',xlab='S96 rep1',ylab='S96 rep2')
title('S96 replicates ChIP-seq scores (normdiff w=local)')






#no background subtraction-MA plot
A=(treat1+treat2)*1/2
M=(treat1-treat2)
plot(A,M,pch='.',xlab='M',ylab='A')
title('S96 MA plot for replicates (raw scores)')




#no background subtraction-MA plot
A=(background1+background2)*1/2
M=(background1-background2)
plot(A,M,pch='.',xlab='M',ylab='A')
title('S96 MA plot for replicates (background subtraction)')



#scaling-MA plot
A=(scale1+scale2)*1/2
M=(scale1-scale2)
plot(A,M,pch='.',xlab='A',ylab='M')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (background scaling) lm')


#normalization MA-plot
M<-(gnormdiff1-gnormdiff2)
A<-(gnormdiff1+gnormdiff2)*1/2
plot(A,M,pch='.',xlab='A',ylab='M')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=all)')





#normalization MA-plot
M<-(normdiff1-normdiff2)
A<-(normdiff1+normdiff2)*1/2
plot(A,M,pch='.',xlab='A',ylab='M')
lm1<-lm(M~A)
lines(A,lm1$fitted)
title('S96 MA plot for replicates (normdiff w=local)')
