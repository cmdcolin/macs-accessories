library(plyr)
library(gplots)
library(R.utils)


### Get wig scores for treat only
getjoinscores<-function(chrnames,t1,t2,currpos) {
  
  ret<-lapply(chrnames, function(chr) {
    if(debug)
      printf("processing %s\tat V%d\n", chr,currpos);
    str<-paste0("V",currpos)
    colnames(t2[[chr]])[2]<-str
    ret<-join(t1[[chr]], t2[[chr]],by="V1")
    ret[complete.cases(ret),]
  })
  names(ret)<-chrnames
  #print(str(ret))
  # from R inferno, Burns (2011)
  #do.call('rbind', ret) 
  ret
}



flatten<-function(chrnames,ret) {
  RTE<-lapply(chrnames, function(chr) {
    if(debug)
      printf("processing %s\n", chr);
    data.frame(chr=chr,ret[[chr]])
  })
  
  do.call('rbind', RTE) 
}

joinWiggleFiles<-function(chrnames,macswig) {
  ret<-getjoinscores(chrnames,macswig[[1]]$control,macswig[[1]]$treat,3)
  currpos<-4
  for(i in 2:length(macswiggle)) {
    ret<-getjoinscores(chrnames,ret,macswig[[i]]$control,currpos)
    ret<-getjoinscores(chrnames,ret,macswig[[i]]$treat,currpos+1)
    currpos<-currpos+2
  }
  flatten(chrnames,ret)
}






doheatmap<-function(table,granularity=1) {
  
  if(granularity!=1) {
    sret<-apply(table[,1:ncol(table)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
  }
  else {
    sret<-table
  }
  
  heatmap.2(as.matrix(sret),col=redgreen(75), scale="none",key=TRUE, 
            density.info='none',trace='none',Rowv=FALSE)
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


getnormdiff<-function(treat,control) {
  
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



ddply(bed,c('chr','start','end'),getPeakNormDiff2,wiggleTable)

getPeakNormDiff2<-function(chr,start,end) {
  chrselect<-strsplit(row[1],'.fsa')[[1]]
  if(chrmatch!=chrselect) {
    chrmatch<<-chrselect
    chrsub<<-scores[scores$chr==chrselect,]
  }
  chrsub[chrsub$pos>row[2]&chrsub$pos<row[3],]
}
getPeakNormDiff<-function(bed,scores) {
  chrmatch="NA"
  chrsub=data.frame()
  ret<-apply(bed,1,function(row) {
    printf("Processing peak %s (%d,%d)\n",row[4],row[2],row[3])
    chrselect<-strsplit(row[1],'.fsa')[[1]]
    if(chrmatch!=chrselect) {
      chrmatch<<-chrselect
      chrsub<<-scores[scores$chr==chrselect,]
    }
    chrsub[chrsub$pos>row[2]&chrsub$pos<row[3],]
  })
  
  # row bind dataframes. from R inferno, Burns (2011)
  do.call('rbind', ret) 
}



slideMean<-function(x,windowsize=100,slide=1){
  idx1<-seq(1,length(x),by=slide);
  idx1+windowsize->idx2;
  idx2[idx2>(length(x)+1)]<-length(x)+1;
  c(0,cumsum(x))->cx;
  return((cx[idx2]-cx[idx1])/windowsize);
}



# Set wiggle table names as abbreviated sample names
prettyNames<-function(wiggleTable) {
  
  fixNames<-sapply(names(macswiggle),function(x) strsplit(x,'-new')[[1]])
  prettyNames<-paste0(sort(rep(fixNames,2)),'-',rep(c('c','t'),nsamples))
  names(wiggleTable)<-c('chr','pos',prettyNames)
}

## Set wiggle table names as chr, pos, V1-Vn
plainNames<-function(wiggleTable) {
  
  fixNames<-c('chr','pos',paste0(rep('V',nsamples),1:(nsamples*2)))
  names(wiggleTable)<-fixNameszx
  
}


# Accessory function for plotting large heatmaps
resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}


