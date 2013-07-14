

doheatmap<-function(table,granularity=1,Rowv=NA,Colv=NA,scale="none",dist=cor) {
  
  if(granularity!=1) {
    sret<-apply(table[,1:ncol(table)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
  }
  else {
    sret<-table
  }
  
  heatmap.2(as.matrix(sret),col=redgreen(75), scale=scale,key=TRUE,
            density.info='none',trace='none',Rowv=Rowv,Colv=Colv)
}

doTreeView<-function(table,file=tempfile(),granularity=1) {
  
  if(granularity!=1) {
    lemm<-seq(1,nrow(table),by=granularity)
    coord<-paste0(table$chr[lemm],'P',table$pos[lemm])
    
    temp<-apply(table[,c(-1,-2)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
    sret<-cbind(coord,temp)
  }
  else {
    sret<-table
  }
  write.table(sret,file=file,row.names=FALSE,quote=FALSE,sep='\t')
}


getNormDiff<-function(treat,control) {
  
  treatSmooth<-slideMean(treat) #default params
  controlSmooth<-slideMean(control)
  
  
  scalingFactor=median(control/treat)
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





# Accessory function for plotting large heatmaps
resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
