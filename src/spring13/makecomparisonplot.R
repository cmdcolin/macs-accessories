#makeComparisonPlot
# main + helper functions

makeComparisonPlot<-function(bed1,bed2,table,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale=FALSE,MA=FALSE) {
  
  bedOverlap<-intersectBed(bed1,bed2)
  bedUnique1<-uniqueBed(bed1,bed2)
  bedUnique2<-uniqueBed(bed2,bed1)
  makeComparisonPlotHelp(bedOverlap,bedUnique1,bedUnique2,table,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale,MA)
}


makeComparisonPlotHelp<-function(bedOverlap,bedUnique1,bedUnique2,table,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale=FALSE,MA=FALSE){
  
  
  wtoverlap<-getPeakScores(bedOverlap,table)
  wtrep1<-getPeakScores(bedUnique1,table)
  wtrep2<-getPeakScores(bedUnique2,table)
  p<-table[sample(1:nrow(table),nrow(table)/10),c(c1,c2)]
  p1<-wtoverlap[,c(c1,c2)]
  p2<-wtrep1[,c(c1,c2)]
  p3<-wtrep2[,c(c1,c2)]
  if(logscale==TRUE) {
    p<-log2(p)
    p1<-log2(p1)
    p2<-log2(p2)
    p3<-log2(p3)
  }
  if(MA==TRUE) {
    m<-p[,1]-p[,2]
    a<-(p[,1]+p[,2])/2
    p[,1]<-a
    p[,2]<-m
    m<-p1[,1]-p1[,2]
    a<-(p1[,1]+p1[,2])/2
    p1[,1]<-a
    p1[,2]<-m
    
    m<-p2[,1]-p2[,2]
    a<-(p2[,1]+p2[,2])/2
    p2[,1]<-a
    p2[,2]<-m
    
    m<-p3[,1]-p3[,2]
    a<-(p3[,1]+p3[,2])/2
    p3[,1]<-a
    p3[,2]<-m
  }
  
  
  plot(p1,pch=19,cex=0.9,col="#aaaaaa22",xlab=xlab,ylab=ylab)
  #,log=if(log)"xy"else""
  points(p1,pch=19,cex=0.9,col=paste0(fillx[1],"77"))
  points(p3,pch=19,cex=0.9,col=paste0(fillx[3],"77"))
  points(p2,pch=19,cex=0.9,col=paste0(fillx[2],"77"))
  legend("bottomright",fill=fillx,legend=legendx)
  title(titlex)
  if(logscale==TRUE&&MA==FALSE){
    
    lines(seq(sqrt(2),9+sqrt(2),length=9)),lwd=3,col=fillx[3])
    lines(seq(-sqrt(2),9-sqrt(2),length(9)),lwd=3,col=fillx[3])
  }
  else if(logscale==TRUE&&MA==TRUE) {
    
    lines(0:9,rep(1,10),lwd=2,col=fillx[3])
    lines(0:9,rep(-1,10),lwd=2,col=fillx[3])
  }
}


makeComparisonPlotHelp2<-function(uniqueTab,overlapTab,table,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale=FALSE,MA=FALSE){
  
  wtrep1<-uniqueTab
  wtoverlap<-getPeakScores(overlapTab,table)
  p<-table[sample(1:nrow(table),nrow(table)/2),c(c1,c2)]
  p1<-wtrep1[,c(c1,c2)]
  p2<-wtoverlap[,c(c1,c2)]
  if(logscale==TRUE) {
    p<-log2(p)
    p1<-log2(p1)
  }
  if(MA==TRUE) {
    m<-p[,1]-p[,2]
    a<-(p[,1]+p[,2])/2
    p[,1]<-a
    p[,2]<-m
    m<-p1[,1]-p1[,2]
    a<-(p1[,1]+p1[,2])/2
    p1[,1]<-a
    p1[,2]<-m
  }
  
  
  plot(p,pch=19,cex=0.9,col=fillx[3],xlab=xlab,ylab=ylab)
  #,log=if(log)"xy"else""
  points(p2,pch=19,cex=0.9,col=paste0(fillx[1],"77"))
  points(p1,pch=19,cex=0.9,col=paste0(fillx[2],"77"))
  legend("bottomright",fill=fillx,legend=legendx)
  title(titlex)
}

