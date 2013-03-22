library(plyr)
library(gplots)
library(R.utils)




# Get chromosomes list
chrnames<-names(macswiggle[[1]]$treat)


nsamples<-length(macswiggle)

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


