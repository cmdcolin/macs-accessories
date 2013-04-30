require(RColorBrewer)

x1<-macswiggle[[1]]$treat$chr01
x2<-macswiggle[[1]]$control$chr01
x3<-macswiggle[[4]]$treat$chr01
c1<-'#000000'
c2<-'#CA5444CC'
c3<-'#5444CACC'
plot(x2[,1],x2[,2],type='l',ylim=c(0,75),col=c1,ylab="Read count",xlab="Genome position (chr01)")
lines(x1[,1],x1[,2],col=c2)
lines(x3[,1],x3[,2],col=c3)
title('Distribution of ChIP-seq and input control read counts (chr01)')
legend('topright',legend=c('Control','ChIP-seq'),fill=c(c1,c2))



x1<-macswiggle[[4]]$treat$chr01
x2<-macswiggle[[4]]$control$chr01
c1<-'#000000'
c2<-'#CA5454CC'
plot(x2[,1],x2[,2],type='l',ylim=c(0,75),col=c1,ylab="Read count",xlab="Genome position (chr01)")
lines(x1[,1],x1[,2],col=c2)
title('Distribution of ChIP-seq and input control read counts (chr01)')
legend('topright',legend=c('Control','ChIP-seq'),fill=c(c1,c2))