

#################
#! Constructor
WiggleClass<-function(name, environ=environment()) {
  nc=list(
    name=name,
    environ=environ,
    scaling=1,
    variance=1,
    spacing=10,
    controlpath=paste(name,'/',name,'_MACS_wiggle/control/',sep=''),
    treatpath=paste(name,'/',name,'_MACS_wiggle/treat/',sep=''),
    controlname=paste(name,'_control_afterfiting_',sep=''),
    treatname=paste(name,'_treat_afterfiting_',sep='')
    )
  
  ########################
  # Read wiggle files from path into memory and assign filesnames
  nc$loadWiggles=function() {
    loadWiggle<-function(wigpath,environ) {
      files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
      for (i in files) {
        file<-paste(wigpath,i,sep='')
        x<-read.table(file, skip=2)
        assign(i,x,inherits=TRUE,envir=environ)
      }
    }
    loadWiggle(nc$treatpath,nc$environ)
    loadWiggle(nc$controlpath,nc$environ)
  }
  #######
  # Get avg reads
  nc$getTotalReads = function(bedfile) {
    getTotalReads<-function(x,filepath){
      chr=x[1]
      start=x[2]
      end=x[3]
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile)
      b=findInterval(start,wig$V1)
      e=findInterval(end,wig$V1)
      sum(wig$V2[b:e])
    }
    # Apply to chip
    apply(bedfile,1,getTotalReads,filepath=nc$treatname)
  } 
  
  # Get avg reads
  nc$getAvgReads = function(bedfile) {
    getAvgReads<-function(x, filepath)
    {  
      chr=x[1]
      start=x[2]
      end=x[3]
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile)
      b=findInterval(start,wig$V1)
      e=findInterval(end,wig$V1)
      sum(wig$V2[b:e])/(e-b);
    }
    #Apply to treat data
    apply(bedfile,1,getAvgReads,filepath=nc$treatname)
  }
  
  nc$getMaxAvgReads<-function(bedfile, window,inc) {
    # Get max average reads over window size
    getMaxAvgReads<-function(x, filepath, window)
    {
      maxreads=array()
      chr=x[1];
      start=as.integer(x[2]);
      end=as.integer(x[3]);
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile);
      bstart=start-window
      bend=end
      maxreads=sapply(seq(bstart,end,by=window),function(p){
        b=findInterval(p,wig$V1)
        e=findInterval(p+window,wig$V1)
        sum(wig$V2[b:e])/(e-b)
      })
      ret=max(maxreads)
      cat('Found max ',ret,' in ',chr,'\n');
      ret
    }
    # Apply to treated data
    apply(bedfile,1,getMaxAvgReads,filepath=nc$treatname,window)
  }

  
  
  
  ####
  # Use all chromosomes for scaling factor
  nc$estimateScalingFactor <- function() {
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    ratio_data=apply(cbind(files1,files2),1,function(x){
      treat=get(x[1])
      control=get(x[2])
      corr=match(treat[,1],control[,1])
      tsig=treat[,2]
      csig=control[,2]
      sapply(corr,function(p){ csig[p]/tsig[p]})
    })
    nc$scaling=median(unlist(ratio_data),na.rm=TRUE)
  }
  
  # Sqrt(Aw+Bw/c), w=all
  nc$estimateVarianceAll<-function() {
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    getSignal<-function(file){get(file)$V2}
    chip_signal=unlist(lapply(files1,getSignal))
    control_signal=unlist(lapply(files2,getSignal))
    #Average signal
    average_chip=mean(chip_signal)
    average_control=mean(control_signal)
    nc$variance=sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  
  ####
  nc$estimateVarianceWindow<-function(xpos, treat,control,ws,corr1,corr2) {
    # find start end positions
    # check for inconsistencies in data
    #
    # checkspacing is correct b window
    #select signals
    #cat(xpos,'\t(', treat$V1[beginning],',',treat$V1[ending],')',nc$matches, '-', nc$unmatches,'\t',
    #    nc$matches1,' ', nc$unmatches1,'(',beginning, ' ', ending, ') (', nc$shifter, ' ', nc$shifter2, '\n')
    b1=xpos[1]-ws
    e1=b1+2*ws
    b2=xpos[2]-ws
    e2=b2+2*ws
    sel1=findInterval(b1:e1,corr1)
    sel2=findInterval(b2:e2,corr2)
    chip_signal=treat$V2[sel1]
    control_signal=control$V2[sel2]
    average_chip=mean(chip_signal,na.rm=TRUE)
    average_control=mean(control_signal,na.rm=TRUE)
    sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  nc$Zxi<-function(x,treat,control,window,corr1,corr2) {
    pos1=x[1]
    pos2=x[2]
    (treat$V2[pos1]-control$V2[pos2]/nc$scaling)/
      max(nc$estimateVarianceWindow(x,treat,control,window[1],corr1,corr2),
          nc$estimateVarianceWindow(x,treat,control,window[2],corr1,corr2),
          nc$variance)
  }
  
  
  
  # Calculate Z scores over all wiggle files
  nc$Z<-function(bedfile, window=c(10,100)) {
    # Get max average reads over window size
    getZscore<-function(x,f1,f2,window){
      chr=x[1];
      start=as.integer(x[2]);
      end=as.integer(x[3]);
      tf=paste(f1,chr,'.wig.gz',sep='')
      cf=paste(f2,chr,'.wig.gz',sep='')
      treat=get(tf)
      control=get(cf)
      
      corr1=findInterval(start:end,treat$V1)
      corr2=findInterval(start:end,control$V1)
      cat(chr,'-\t(',start,',',end, ')\n')
      app=cbind(corr1,corr2)
      app=app[-head(app,max(window)),]
      app=app[-tail(app,max(window)),]
      res<-apply(app,1,nc$Zxi,treat,control,window,corr1,corr2)
      cbind(app,res)
    }
    
    apply(bedfile,1,getZscore,f1=nc$treatname,f2=nc$controlname,window)
  }
  
  nc$getMaxAvgZscore<-function(Zscore,ws=10) {
    acc1<-function(j,ws,vpos,vsig) {
      # binary search
      b=findInterval(j,vpos)
      e=findInterval(j+ws,vpos)
      mean(vsig[b:e]);
    }
    acc2<-function(x,window,acc1){
      vpos=x[,1]
      vsig=x[,3]
      bstart=head(vpos,1)
      bend=tail(vpos,1)
      if(bstart==bend) { 0}
      else{
      cat(bstart,' ',bend, ' ', ws,'\n')
      windows=seq(bstart,bend-ws/nc$spacing,by=ws)
      reads=sapply(windows, acc1,ws,vpos,vsig)
      max(reads)
      }
    }
    sapply(Zscore,acc2,ws,acc1)
  }
  
  
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}




###########


######Sketches
#
# Get peak index from bed
#  (uneeded) match(bed1,bed2)
###########
#getPeakIndex<-function(bedfile)
#{
#  index=array()
#  for(i in 1:length(bedfile$V1)) {
#    index[i]=as.integer(unlist(strsplit(as.character(bedfile$V4[i]),'_'))[3])
#  }
#  index
#}
#########
#
#nc$overlapBed(b) {
#  s=paste('intersectBed -a ',nc$name,'/',nc$name,'_peaks.bed')
#  s=paste(s, '-b ', b$name, '/', b$name, '_peaks.bed -wa > ')
#  s=paste(nc$name, '/', nc$name,'_overlap.bed')
#  system(s)
#}
#nc$uniqueBed(b) {
#  s=paste('subtractBed -a ',nc$name,'/',nc$name,'_peaks.bed')
#  s=paste(s, '-b ', nc$name, '/', nc$name, '_overlap.bed > ')
#  s=paste(a$name, '/', a$name,'_unique.bed')
#  system(s)
#}
#plotnew(x1,x2,b1,b2,b3,name1,name2,func) {
#  id1=match(b2$V4,b1$V4)
#  id2=match(b3$V4,b1$V4)
#  r1=hs959$func(b1)
#  r2=s96$func(b1)
#  
#}

