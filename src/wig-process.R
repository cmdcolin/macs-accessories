wig1=wig[1]
wig2=wig[2]

#plot(wig1$treat$chr01)
#plot(wig2)


x=findInterval(wig1[,1],wig[,2])
plot(wig1[,2],wig2[x,2])


s96rep1<-read.table('S96rep1_normdiff_old.txt')
s96rep2<-read.table('S96rep2_normdiff_old.txt',skip=1)
names(s96rep1)<-c("blank","chr","start","end","score")
names(s96rep2)<-c("blank","chr","start","end","score")
s96rep1.chr01=s96rep1[s96rep1$chr=="chr01",]
s96rep2.chr01=s96rep1[s96rep2$chr=="chr01",]
plot(s96rep1.chr01$start,s96rep1.chr01$score,type="l")






match=findInterval(s96rep1.chr01,s96rep2.chr01[,3])
smoothScatter(x1c1[,5],x2c1[match,5],xlab="S96 rep1",ylab="S96 rep2")
title('Comparison of NormDiff scores across replicates')


p23<-function(chr) {
  
  x1c1=x1[x1[,2]==chr,]
  x2c1=x2[x2[,2]==chr,] 
  b1=c1[c1[,1]==chr,]
  b2=c2[c2[,1]==chr,]
  b3=c3[c3[,1]==chr,]
  
  
  doplot<-function(b,c) {
    apply(b,1,function(x) {
      set.seed(1234)
      s=as.numeric(x[2])-50
      e=as.numeric(x[3])+50
      m1=x1c1[,3]>s & x1c1[,3]<e
      m2=findInterval(x1c1[m1,3],x2c1[,3])
      l1=mean(x1c1[m1,5],na.rm=TRUE)
      l2=mean(x2c1[m2,5],na.rm=TRUE)
      points(l1,l2,col=c,pch=19,cex=0.5)
    })
  }
  
  
  doplot(b1,'red')
  doplot(b2,'green')
  doplot(b3,'blue')
}


sapply(as.character(unique(x1[,2])),p23)


x2c1=x2[x2[,2]=="chr01",]
plot(x2c1[,3],x2c1[,5],type="l")

c1=read.table('newdata/S96rep1_peaks.bed')
c2=read.table('newdata/S96rep2_peaks.bed')
c3=read.table('newdata/int1.bed')
b1=c1[c1[,1]=='chr01',]
b2=c2[c2[,1]=='chr01',]
b3=c3[c3[,1]=='chr01',]


doplot<-function(b,c) {
  apply(b,1,function(x) {
    set.seed(1234)
    s=as.numeric(x[2])-100
    e=as.numeric(x[3])+100
    printf("%d\t%d\n",s,e)
    m1=x1c1[,3]>s & x1c1[,3]<e
    m2=findInterval(x1c1[m1,3],x2c1[,3])
    l1=mean(x1c1[m1,5],na.rm=TRUE)
    l2=mean(x2c1[m2,5],na.rm=TRUE)
    printf("%f\t%f\n",l1,l2)
    points(l1,l2,col=c,pch=19,cex=0.5)
  })
}


doplot(b1,'red')
doplot(b2,'green')
doplot(b3,'blue')
