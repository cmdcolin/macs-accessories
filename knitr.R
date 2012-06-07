debug=TRUE
loadMacsEnv<-function(name1,name2) {
  local({
    wig1=WiggleClass(name1)
    wig2=WiggleClass(name2)
    wig1$loadWiggles() 
    wig2$loadWiggles()
    ###
    wig1$peaks=read.table(paste(name1,'/',name1,'_peaks.bed',sep=''))
    wig1$shared=read.table(paste(name1,'/',name1,'_overlap.bed',sep=''))
    wig1$unique=read.table(paste(name1,'/',name1,'_unique.bed',sep=''))
    wig2$peaks=read.table(paste(name2,'/',name2,'_peaks.bed',sep=''))
    wig2$shared=read.table(paste(name2,'/',name2,'_overlap.bed',sep=''))
    wig2$unique=read.table(paste(name2,'/',name2,'_unique.bed',sep=''))
    # Match indexes 
    environment()
  })
}

#setwd('macs1.4.1')
#e1=loadMacsEnv('S96','HS959')
#setwd('..')
#setwd('macs-bad-run')
#e3=loadMacsEnv('S96','HS959')
#setwd('..')

setwd('macs1.4.2')
e2=loadMacsEnv('S96','HS959')
setwd('..')
attach(e2)

plotTotalReads<-function(t, w1, w2,c1,c2) {
  r1=w1$getTotalReads(w1$peaks) 
  r2=w2$getTotalReads(w1$peaks)
  plot(r1,r2,xlab=paste(w1$name,'reads'),ylab=paste(w2$name, 'reads'),pch='*') 
  id1=match(wig1$shared$V4,wig1$peaks$V4) 
  id2=match(wig1$unique$V4,wig1$peaks$V4)
  points(r1[id1],r2[id1],pch=1,col=c1) 
  points(r1[id2],r2[id2],pch=1,col=c2) 
  title(t)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
}
plotTotalReads('Total S96 peak reads vs HS959 synteny MACS 1.4.1 bw=25',wig1,wig2,'pink','green')
plotTotalReads('Total HS959 peak reads vs S96 synteny MACS 1.4.2 bw=25',wig2,wig1,'lightblue','green')


plotMaxAvgReads<-function(t, w1, w2, bed,c1,c2) {
  ###########
  # Use Max avg reads over windows
  rma1=w1$getMaxAvgReads(w1$peaks,100)
  rma2=w2$getMaxAvgReads(w1$peaks,100)
  ###################### 
  plot(rma1,rma2,xlab=paste('Max Avg', w1$name,'reads'),ylab=paste('Max Avg',w2$name,'reads'),pch='*')
  id1=match(wig1$shared$V4,wig1$peaks$V4) 
  id2=match(wig1$unique$V4,wig1$peaks$V4)
  points(rma1[id1],rma2[id1],pch=1,col=c1)
  points(rma1[id2],rma2[id2],pch=1,col=c2)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}

#####
# e2 is the latest and greatest run....
plotMaxAvgReads('Max Average reads S96 peaks  w=100', wig1,wig2,wig1peak,'orange','lightgreen')
plotMaxAvgReads('Max Average reads HS959 peaks  w=100', wig2,wig1,wig2peak,'yellow','lightgreen')

