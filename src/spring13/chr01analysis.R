

r1<-macswiggle[[1]]$treat$chr01
r2<-macswiggle[[2]]$treat$chr01
c1<-macswiggle[[1]]$control$chr01
c2<-macswiggle[[2]]$control$chr01

match=findInterval(r1$V1,r2$V1)
plot(r1$V2,r2$V2,pch='.')

# Background subtraction
cmatch=findInterval(r1$V1,c1$V1)
cmatch2=findInterval(r2$V1,c2$V1)

p1<-r1$V2-c1[cmatch,2]
p2<-r2$V2-c2[cmatch2,2]
plot(p1,p2[match],pch='.')
