

getNormDiff<-function(treat,control) {
  
  treatSmooth<-slideMean(treat) #default params
  controlSmooth<-slideMean(control)
  
  scalingFactor=median(control/treat,na.rm=TRUE)
  globalVariance=sqrt(mean(treat)+mean(control)/scalingFactor^2)
  varlist<-sqrt(treatSmooth+controlSmooth/scalingFactor^2)
  
  
  #normdiff local
  normVar<-sapply(varlist,function(x){
    max(x,globalVariance)
  })
  
  (treat-control/scalingFactor)/normVar
}

getPeakScores<-function(bed,scores,do.rbind=TRUE) {
  chrsplit<-split(scores,factor(scores$chr))
  chrselect<-""
  ret<-apply(bed,1,function(row) {
    start=as.numeric(row['start'])
    end=as.numeric(row['end'])
    
    chrselect<-strsplit(row[['chr']],"\\.")[[1]][1]
    if(debug) {
      printf("Processing %s (%d,%d)\n",chrselect,start,end)
    }
    chrsub<-chrsplit[[chrselect]]
    chrsub[chrsub$pos>start&chrsub$pos<end,]
  })
  
  # from R inferno, Burns (2011)
  if(do.rbind)
    do.call('rbind', ret) 
  else
    ret
}

slideMean<-function(x,windowsize=100,slide=1){
  idx1<-seq(1,length(x),by=slide);
  idx1+windowsize->idx2;
  idx2[idx2>(length(x)+1)]<-length(x)+1;
  c(0,cumsum(x))->cx;
  return((cx[idx2]-cx[idx1])/windowsize);
}

