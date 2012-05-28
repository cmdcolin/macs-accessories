########################
# Read wiggle files from path into memory and assign filesnames
loadWiggle<-function(wigpath) {
  files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
  for (i in files) {
    file<-paste(wigpath,i,sep='')
    x<-read.table(file, skip=2)
    assign(i,x,inherits=TRUE)
  }
}


#################
#! Constructor
WiggleClass<-function(name) {
  nc=list(
    name=name,
    scaling=1,
    variance=1,
    controlpath=paste(name,'/',name,'_MACS_wiggle/control/',sep=''),
    treatpath=paste(name,'/',name,'_MACS_wiggle/treat/',sep=''),
    controlname=paste(name,'_control_afterfiting_',sep=''),
    treatname=paste(name,'_treat_afterfiting_',sep='')
    )
  
  nc$loadWiggles=function() {
    loadWiggle(nc$treatpath)
    loadWiggle(nc$controlpath)
  }
  #######
  # Get avg reads
  nc$getTotalReads = function(bedfile) {
    gtr<-function(x,filepath){
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
    apply(bedfile,1,gtr,filepath=nc$treatname)
  } 
  
  # Get avg reads
  nc$getAvgReads = function(bedfile) {
    gar<-function(x, filepath)
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
    apply(bedfile,1,gar,filepath=nc$treatname)
  }
  nc$getMaxAvgReads<-function(bedfile, wsize) {
    # Get max average reads over window size
    getMaxAvgReads<-function(x, filepath, wsize)
    {
      maxreads=array()
      chr=x[1];
      start=as.integer(x[2]);
      end=as.integer(x[3]);
      wigfile=paste(filepath,chr,'.wig.gz',sep='')
      wig=get(wigfile);
      bstart=start-wsize/10
      bend=end+wsize/10
      maxreads=sapply(seq(bstart,end,by=wsize),function(p){
        b=findInterval(p,wig$V1)
        e=findInterval(p+wsize,wig$V1)
        sum(wig$V2[b:e])/(e-b)
      })
      ret=max(maxreads)
      cat('Found max ',ret,' in ',chr,'\n');
      ret
    }
    # Apply to treated data
    apply(bedfile,1,getMaxAvgReads,filepath=nc$treatname,wsize=wsize)
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
  # Sqrt(Aw+Bw/c), w=1
  nc$estimateVarianceWindow<-function(treat,control,pos,wsize) {
    b=findInterval(pos-wsize,treat$V1)
    e=findInterval(pos+wsize,treat$V1)
    corr=match(treat$V1[b:e],control$V1)
    chip_signal=treat$V2[b:e]
    control_signal=control$V2[corr]
    average_chip=mean(chip_signal,na.rm=TRUE)
    average_control=mean(control_signal,na.rm=TRUE)
    sqrt(average_chip+average_control/nc$scaling^2)
  }
  
  nc$Zxi<-function(treat,control,pos,wsize) {
    (treat$V2[pos]-control$V2[pos]/nc$scaling)/
      max(nc$estimateVarianceWindow(treat,control,pos,wsize[1]),
          nc$estimateVarianceWindow(treat,control,pos,wsize[2]),
          nc$variance)
  }
  
  
  
  
  ############
  # Calculate Z scores over all wiggle files
  nc$Z<-function(wsize=c(10,100)) {

    files1=list.files(nc$treatpath,"*.fsa.wig.gz")
    files2=list.files(nc$controlpath,"*.fsa.wig.gz")
    
    apply(cbind(files1,files2),1,function(x){
      treat<-get(x[1])
      control<-get(x[2])
      cat(x[1],'\n')
      
      start=findInterval(head(treat$V1,1)+max(wsize),treat$V1)
      end=findInterval(tail(treat$V1,1)-max(wsize),treat$V1)
      V1<-treat$V1[start:end]
      V2<-sapply(start:end,function(x){
        nc$Zxi(treat,control,x,wsize)
      })
      cat(x[1],' finished\n')
      as.table(cbind(V1,V2))
    })
  }
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}





###########
# Get average Z
getAvgZ<-function(bedfile,Zscore) {
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    name=paste(Zscore[['name']],'_treat_afterfiting_',chr,'.wig.gz',sep='')
    name=paste('Z', name)
    
    Zchr=Zscore[[name]]
    Zchrpos=as.array(Zchr[,1])
    Zchrpeak=as.array(Zchr[,2])
    
    peak=Zchrpeak[Zchrpos>start & Zchrpos<end];
    reads[i]=mean(peak,na.rm=TRUE)
  }
  reads
}



###########
# Get average Z
getMaxAvgZ<-function(bedfile,Zscore,wsize) {
  reads=array();
  chrold=''
  Zchr=NULL
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    # Avoid reloading env variables
    if(as.character(chrold)!=as.character(chr)){
      name=paste(Zscore[['name']],'_treat_afterfiting_',chr,'.wig.gz',sep='')
      name=paste('Z', name)
      
      Zchr=Zscore[[name]]
      print(name)
      chrold=chr
    }    
    
    bstart=start-wsize/10
    bend=end+wsize/10
    maxreads=array()
    # Avoid extract$column in loop
    wigpos=Zchr[,1]
    wigpeak=Zchr[,2]
    
    for(j in seq(bstart,end,by=100)) {
      # binary search
      b=findInterval(j,wigpos)
      e=findInterval(j+wsize,wigpos)
      # Peak window
      peak=wigpeak[b:e];
      maxreads[j]=mean(peak,na.rm=TRUE);
    }
    reads[i]=max(maxreads,na.rm=TRUE);
    print(cat('Found max reads ', reads[i], ' at ', b, ' ', e,'. Used ', (end-bstart)/100, ' windows'))
    
  }
  reads
}



###########
# Get average Z
getAvgZ_WholeGenome<-function(Zscore,wsize) {
  reads=array();
  for(i in Zscore) {
    if(is.null(dim(i))) {next;}
    wigpos=i[,1]
    wigpeak=i[,2]
    bstart=wsize
    bend=length(wigpos)-wsize
    maxreads=array()
    for(j in seq.int(bstart,bend,by=wsize)) {
      # binary search
      b=findInterval(j,wigpos)
      e=findInterval(j+wsize,wigpos)
      
      # Peak window
      peak=wigpeak[b:e];
      
      maxreads[j/wsize]=mean(peak,na.rm=TRUE);
    }
    reads=c(reads,maxreads)
  }
  reads
}




############
# Old get avg Z
getAvgNormDiff<-function(bedfile, path1, path2, scaling_factor, variance) {
  
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    
    # Reads from S96
    wigfile1=paste(path1,chr,'.wig.gz',sep='');
    wig1=get(wigfile1);
    treat_peak=wig1[wig1$V1>start & wig1$V1<end,];
    wigfile2=paste(path2,chr,'.wig.gz',sep='');
    wig2=get(wigfile2);
    control_peak=wig2[wig2$V1>start & wig2$V1<end,];
  
    treat_peak=as.array(treat_peak$V2)
    control_peak=as.array(control_peak$V2)
    Z=array()
    
    for(j in 1:(end-start)) {
      Z[j]= (treat_peak[j]-control_peak[j]/scaling_factor)/sqrt(variance)
    }
    reads[i]=mean(Z,na.rm=TRUE);
  }
  reads
}






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

