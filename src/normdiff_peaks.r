nfiles=list.files(pattern="*_normdiff.txt")
bfiles=list.files(pattern="*_peaks.bed")
ret=lapply(nfiles,function(x) read.table(x))
ret2=lapply(bfiles,function(x) read.table(x))
hrep1=ret[[1]]
hrep1_peaks=ret2[[1]]
chr1=hrep1[hrep1$V1=='chr01',]
chr1_peaks=hrep1_peaks[hrep1_peaks$V1=='chr01',]
plot(chr1$V2,chr1$V4,type='l')
ret=apply(chr1_peaks,1,function(x){
  rect(x[2],-100,x[3],100,col=rgb(1,0,0,0.5))
})



require(R.utils)
qqnd<-function(normdiff) {
  clt=sapply(1:5000,function(x){mean(sample(normdiff,50))})
  nd.mu=mean(clt)
  nd.sd=sd(clt)
  qqnorm(clt)
  qqline(rnorm(1000,mean=nd.mu,sd=nd.sd))
}
qqnd(normdiff=hrep1$V4)

densitynd<-function(normdiff) {
  clt=sapply(1:10000,function(x){mean(sample(normdiff,50))})
  nd.mu=mean(clt)
  nd.sd=sd(clt)
  plot(function(x){dnorm(x,mean=nd.mu,sd=nd.sd)},-3,3)
  lines(density(clt))
  printf("%f\t%f",nd.mu,nd.sd)
  #lines(density(rnorm(1000,nd.mu,nd.sd)))
}
densitynd(normdiff=hrep1$V4)



#################
#get all chr names
#loop from seq 0 to maxchr length by window size
#match intervals in samples for samples in location
#average normdiff scores
#plot pair on plot
##
nd1=ret[[1]]
nd2=ret[[2]]
str=c('chr','start','end','score')
names(nd1)<-str
names(nd2)<-str
chrs=levels(nd1$chr)
results=lapply(chrs,function(name){
  sel1=nd1[nd1$chr==name,]
  sel2=nd2[nd2$chr==name,]
  x=sapply(seq(1,max(sel1$start),by=100),function(i) {
    x1=sel1$start[i:(i+100)]
    x2=findInterval(x1,sel2$start,all.inside=TRUE)
    m1=mean(sel1[x1,]$score)
    m2=mean(sel2[x2,]$score)
    cbind(m1,m2)
  }) 
})
