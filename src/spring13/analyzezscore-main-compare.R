
wiggles<-getscoresmod(chrnames,macswiggle[[3]]$treat,macswiggle[[3]]$control,macswiggle[[1]]$treat,macswiggle[[1]]$control)



bed1<-read.table('s96rep1-high_peaks.bed')
bed2<-read.table('hs959rep1-new_peaks.bed')
bed3<-read.table('s96-intersect.bed')
names(bed1)<-c('chr','start','end','name','num')
names(bed2)<-c('chr','start','end','name','num')
names(bed3)<-c('chr','start','end','name','num')
bed1$start=as.numeric(bed1$start)
bed1$end=as.numeric(bed1$end)
bed2$start=as.numeric(bed2$start)
bed2$end=as.numeric(bed2$end)
bed3$start=as.numeric(bed3$start)
bed3$end=as.numeric(bed3$end)


getpeaknormdiffmax<-function(bed,scores) {
  chrmatch=""
  chrsub=data.frame()
  ret<-apply(bed,1,function(row) {
    s=as.numeric(row[2])
    e=as.numeric(row[3])
    printf("Processing peak %s (%d,%d)\n",row[4],s,e)
    chrselect<-strsplit(row[1],'.fsa')[[1]]
    if(chrmatch!=chrselect) {
      printf("Here %s %s\n",chrmatch,chrselect)
      chrmatch<<-chrselect
      chrsub<<-scores[scores$chr==chrselect,]
    }
    
    
    chrsub2<-chrsub[chrsub$pos>as.numeric(row[2])&chrsub$pos<as.numeric(row[3]),]
    x<-(apply(chrsub2[,3:6],2,function(x){
      #print(dim(chrsub2))
      if(length(x)<10) mean(x)
      else slideMean(x,10,10)
    }))
    y<-data.frame(x[,1],x[,2],x[,3],x[,4])
    #chrsub2
  })
  
  # from R inferno, Burns (2011)
  do.call('rbind', ret) 
}

ret<-getpeaknormdiffmax(bed1,wiggles)
ret2<-getpeaknormdiffmax(bed2,wiggles)
ret3<-getpeaknormdiffmax(bed3,wiggles)


plot(wiggles$treat1,wiggles$treat2,pch='.',cex=2,xlab='S96rep1',ylab='HS959rep1')
title('Raw read counts of S96 vs HS959')
points(ret[,1],ret[,3],col=2,pch='.',cex=3)
points(ret2[,1],ret2[,3],col=3,pch='.',cex=3)
points(ret3[,1],ret3[,3],col=rgb(0,0,1,0.5),pch='.',cex=3)
legend("topleft",legend=c('Overlap','S96 unique','HS959 unique'),fill=c('blue','red','green'))

####mod


plot(wiggles$treat1,wiggles$treat2,pch='.',cex=2,xlab='S96rep1',ylab='HS959rep1')
title('Raw read counts of S96 vs HS959')
points(ret[,1+2],ret[,3+2],col=2,pch='.',cex=3)
points(ret2[,1+2],ret2[,3+2],col=3,pch='.',cex=3)
points(ret3[,1+2],ret3[,3+2],col=4,pch='.',cex=3)
legend("topleft",legend=c('Overlap','S96 unique','HS959 unique'),fill=c('blue','red','green'))








ret<-getpeaknormdiffmax(bed1,wiggles)
ret2<-getpeaknormdiffmax(bed2,wiggles)
ret3<-getpeaknormdiffmax(bed3,wiggles)


plot(wiggles$treat1-wiggles$control1/m1,wiggles$treat2-wiggles$control2/m2,pch='.',cex=2,xlab='S96rep1',ylab='HS959rep1')
title('Background scaled reads S96 vs HS959')
points(ret[,1+]-ret[,2]/m1,ret[,3]-ret[,4]/m2,col=2,pch='.',cex=3)
points(ret2[,1]-ret2[,2]/m1,ret2[,3]-ret2[,4]/m2,col=3,pch='.',cex=3)
points(ret3[,1]-ret3[,2]/m1,ret3[,3]-ret3[,4]/m2,col=4,pch='.',cex=3)
legend("topleft",legend=c('Overlap','S96 unique','HS959 unique'),fill=c('blue','red','green'))



plot(wiggles$treat1-wiggles$control1/m1,wiggles$treat2-wiggles$control2/m2,pch='.',cex=2,xlab='S96rep1',ylab='HS959rep1')
title('Background scaled reads S96 vs HS959')
points(ret[,1+2]-ret[,2+2]/m1,ret[,3+2]-ret[,4+2]/m2,col=2,pch='.',cex=3)
points(ret2[,1+2]-ret2[,2+2]/m1,ret2[,3+2]-ret2[,4+2]/m2,col=3,pch='.',cex=3)
points(ret3[,1+2]-ret3[,2+2]/m1,ret3[,3+2]-ret3[,4+2]/m2,col=4,pch='.',cex=3)
legend("topleft",legend=c('Overlap','S96 unique','HS959 unique'),fill=c('blue','red','green'))