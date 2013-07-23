

getNormDiff<-function(treat,control,slideWindow=1,slideWindow2=10) {
  
  treatSmooth<-slideMean(treat,slideWindow) #default params
  controlSmooth<-slideMean(control,slideWindow)
  treatSmooth2<-slideMean(treat,slideWindow2) #default params
  controlSmooth2<-slideMean(control,slideWindow2)
  
  
  
  scalingFactor=median(control/treat,na.rm=TRUE)
  globalVariance=sqrt(mean(treat)+mean(control)/scalingFactor^2)
  varlist<-sqrt(treatSmooth+controlSmooth/scalingFactor^2)
  varlist2<-sqrt(treatSmooth2+controlSmooth2/scalingFactor^2)
  
  
  #normdiff local
  normVar<-unlist(apply(cbind(varlist,varlist2),1,function(x){
    max(x[1],x[2],globalVariance)
  }))
  
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

