loadMacsEnv<-function(name1,name2) {
  local({
    wig1=WiggleClass(name1)
    wig2=WiggleClass(name2)
    #wig1$loadWiggles() 
    #wig2$loadWiggles()
    ###
    wig1$peaks=read.table(paste(name1,'/',name1,'_peaks.bed',sep=''))
    wig2$peaks=read.table(paste(name2,'/',name2,'_peaks.bed',sep=''))
    wig1$shared=intersectBed(wig1,wig2)
    wig1$unique=uniqueBed(wig1,wig2)
    wig2$shared=intersectBed(wig2,wig1)
    wig2$unique=uniqueBed(wig2,wig1)
    # Return environment
    environment()
  })
  
}

plotTotalReads<-function(t, w1, w2,c1,c2) {
  r1=w1$getTotalReads(w1$peaks) 
  r2=w2$getTotalReads(w1$peaks)
  plot(r1,r2,xlab=paste(w1$name,'reads'),ylab=paste(w2$name, 'reads'),pch='*') 
  id1=match(w1$shared$V4,w1$peaks$V4) 
  id2=match(w1$unique$V4,w1$peaks$V4)
  points(r1[id1],r2[id1],pch=1,col=c1) 
  points(r1[id2],r2[id2],pch=1,col=c2) 
  title(t)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
}




plotMaxAvgReads<-function(t, w1, w2,c1,c2) {
  ###########
  # Use Max avg reads over windows
  rma1=w1$getMaxAvgReads(w1$peaks,100)
  rma2=w2$getMaxAvgReads(w1$peaks,100)
  ###################### 
  plot(rma1,rma2,xlab=paste('Max Avg', w1$name,'reads'),ylab=paste('Max Avg',w2$name,'reads'),pch='*')
  id1=match(w1$shared$V4,w1$peaks$V4) 
  id2=match(w1$unique$V4,w1$peaks$V4)
  points(rma1[id1],rma2[id1],pch=1,col=c1)
  points(rma1[id2],rma2[id2],pch=1,col=c2)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}



plotAvgZscore<-function(t, w1, w2, wz1, wz2, c1,c2) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  Z1<-sapply(wz1,mean)
  Z2<-sapply(wz2,mean)
  plot(Z1,Z2,xlab=paste('Max Avg', w1$name,'Zscore'),ylab=paste('Max Avg',w2$name,'Zscore'),pch='*')
  id1=match(w1$shared$V4,w1$peaks$V4) 
  id2=match(w1$unique$V4,w1$peaks$V4)
  points(Z1[id1],Z2[id1],pch=1,col=c1)
  points(Z1[id2],Z2[id2],pch=1,col=c2)
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
  ret=list()
  ret[['1shared']]=Z1[id1]
  ret[['1unique']]=Z1[id2]
  ret[['2shared']]=Z2[id1]
  ret[['2unique']]=Z2[id2]
  ret
}

plotMaxAvgZscore<-function(t, w1, w2, wz1,wz2, c1,c2) {
  ##########
  # Get Z scoresmn/.,mnb,.,..,m.,
  maxw1<-w1$getMaxAvgZscore(wz1)
  maxw2<-w2$getMaxAvgZscore(wz2)
  id1=match(w1$shared$V4,w1$peaks$V4) 
  id2=match(w1$unique$V4,w1$peaks$V4)
  plot(maxw1,maxw2,pch='*',xlab=paste(w1$name,'peak'),ylab=paste(w2$name,'syntenic'))
  points(maxw1[id1],maxw2[id1],col=c1)
  points(maxw1[id2],maxw2[id2],col=c2)
  #legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
  ret=list()
  ret[['1shared']]=maxw1[id1]
  ret[['1unique']]=maxw1[id2]
  ret[['2shared']]=maxw2[id1]
  ret[['2unique']]=maxw2[id2]
  ret
}

plotSortedMaxAvgZscore<-function(t, w1, w2, r, c1,c2) {
  
  wz1sort=sort(r[['2shared']])
  wz2sort=sort(r[['2shared']])
  idx1temp=match(r[['1']],wz1sort)
  idx2temp=match(r[['2']],wz2sort)
  id1=match(w1$shared$V4,w1$peaks$V4) 
  id2=match(w1$unique$V4,w1$peaks$V4)
  idx1=match(id1,idx1temp)
  idx2=match(id2,idx2temp)
  xs=1:length(wz2sort)
  xs2=1:length(wz1sort)
  plot(xs,wz2sort,pch='.',xlab='Rank',ylab='Avg NormDiff')
  points(xs[idx1],wz2sort[idx1],col=c1,pch='.')
  points(xs[idx2],wz2sort[idx2],col=c2,pch='.')
  polygon(c(xs2,rev(xs)),c(wz1sort,rev(wz2sort)),col='lightyellow',border=FALSE)
  points(xs2[idx1],wz1sort[idx1],col=c1,pch='.')
  points(xs2[idx2],wz1sort[idx2],col=c2,pch='.')
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}

plotSortedMaxAvgZscoreX<-function(t, w1, w2, r, c1,c2) {
  
  wz1sort=sort(r[['2shared']])
  wz2sort=sort(r[['2unique']])
  xs1=1:length(wz1sort)
  xs2=1:length(wz2sort)
  plot(xs1,wz1sort,pch='.',xlab='Rank',ylab='Avg NormDiff')
  points(xs1,wz1sort,col=c1,pch='.')
  points(xs2,wz2sort,col=c2,pch='.')
  #polygon(c(xs1,rev(xs)),c(wz1sort,rev(wz2sort)),col='lightyellow',border=FALSE)
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}
plotSortedMaxAvgZscoreOld<-function(t, w1, w2, r, c1,c2) {
  
  wz1sort=sort(r[['1']])
  wz2sort=sort(r[['2']])
  id1=match(w1$shared$V4,w1$peaks$V4) 
  id2=match(w1$unique$V4,w1$peaks$V4)
  xs=1:length(wz2sort)
  xs2=1:length(wz1sort)
  plot(xs,wz2sort,pch='.',xlab='Rank',ylab='Avg NormDiff')
  points(xs[id1],wz2sort[id1],col=c1,pch='.')
  points(xs[id2],wz2sort[id2],col=c2,pch='.')
  polygon(c(xs2,rev(xs)),c(wz1sort,rev(wz2sort)),col='lightyellow',border=FALSE)
  points(xs2[id1],wz1sort[id1],col=c1,pch='.')
  points(xs2[id2],wz1sort[id2],col=c2,pch='.')
  legend('bottomright', legend=c('shared', 'unique'), fill=c(c1, c2))
  title(t)
}


#setwd('macs1.4.2')
#e2=loadMacsEnv('S96','HS959')
#setwd('..')
#attach(e2)
#plotTotalReads('Total S96 peak reads vs HS959 synteny MACS 1.4.1 bw=25',wig1,wig2,'pink','green')
#plotTotalReads('Total HS959 peak reads vs S96 synteny MACS 1.4.2 bw=25',wig2,wig1,'lightblue','green')



#####
# e2 is the latest and greatest run....
#plotMaxAvgReads('Max Average reads S96 peaks  w=100', wig1,wig2,'orange','lightgreen')
#plotMaxAvgReads('Max Average reads HS959 peaks  w=100', wig2,wig1,'yellow','lightgreen')
