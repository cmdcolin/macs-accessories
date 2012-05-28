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
    gamr<-function(x, filepath, wsize)
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
    apply(bedfile,1,gamr,filepath=nc$treatname,wsize=wsize)
  }
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
  
  
  
  ####
  # Use all chromosomes for scaling factor
  nc$estimateScalingFactor <- function() {
    files1=list.files(nc$treatpath,pattern="*.fsa.wig.gz")
    files2=list.files(nc$controlpath,pattern="*.fsa.wig.gz")
    ratio_data=apply(cbind(files1,files2),1,function(x){
      cat(x[1], ' ',x[2],'\n')
      treat=get(x[1])
      control=get(x[2])
      corr=match(treat[,1],control[,1])
      tsig=treat[,2]
      csig=control[,2]
      sapply(corr,function(p){ csig[p]/tsig[p]})
    })
    #Fix this crap
    median(sapply(ratio_data,function(x){median(x,na.rm=TRUE)}))
  }
  
  
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}











####
# Use all chromosomes for scaling factor
estimate_scaling_factor <- function(path1,path2) {
  files1=list.files(treatpath,pattern="*.fsa.wig.gz")
  files2=list.files(controlpath,pattern="*.fsa.wig.gz")
  ratio_data=array()
  for (i in 1:length(files1)) {
    treat=get(files1[i])
    control=get(files2[i])
    corr=match(control[,1],treat[,1])
    ratios=sapply(corr,function(x){control[,2]/treatsig[x,2]})
    ratio_data=c(ratio_data,ratios)
  }
  median(ratio_data,na.rm=TRUE)
}


# Sqrt(Aw+Bw/c), w=all
estimate_variance_all<-function(path1,path2,scaling_factor) {
  files1=list.files(path=path1,pattern="*.fsa.wig.gz")
  files2=list.files(path=path2,pattern="*.fsa.wig.gz")
  ratio_data=array()
  chip_signal=array()
  control_signal=array()
  for (i in 1:length(files1)) {
    treat<-get(files1[i])
    control<-get(files2[i])
    chip_signal=c(chip_signal,treat[,2])
    control_signal=c(control_signal,control[,2])
  }
  average_chip=mean(chip_signal,na.rm=TRUE)
  average_control=mean(control_signal,na.rm=TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}


# Sqrt(Aw+Bw/c), w=1
estimate_variance_one<-function(treat,control,scaling_factor,pos) {
  chip_signal=as.array(t(treat[(pos-1):(pos+1)]))
  control_signal=as.array(t(control[(pos-1):(pos+1)]))
  average_chip=mean(chip_signal,na.rm=TRUE)
  average_control=mean(control_signal,na.rm=TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}

# Sqrt(Aw+Bw/c), w=10
estimate_variance_ten<-function(treat,control,scaling_factor,pos) {
  chip_signal=as.array(t(treat[(pos-10):(pos+10)]))
  control_signal=as.array(t(control[(pos-10):(pos+10)]))
  average_chip=mean(chip_signal,na.rm=TRUE)
  average_control=mean(control_signal,na.rm=TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}




Zxi<-function(treat,control,scaling_factor,pos,varianceall) {
  (treat[pos]-control[pos]/scaling_factor)/
      max(estimate_variance_one(treat,control,scaling_factor,pos),
          estimate_variance_ten(treat,control,scaling_factor,pos),
          varianceall)
}



############
# Calculate Z scores over all wiggle files
Z<-function(treat,control,scaling_fact,variance_all,name) {

  files1=list.files(path=treat,pattern="*.fsa.wig.gz")
  files2=list.files(path=control,pattern="*.fsa.wig.gz")
  Z=list()
  Z[['name']]=name
  for (i in 1:length(files1)) {
    treat<-get(files1[i])
    control<-get(files2[i])
    print(files1[i])
    start=10
    end=length(treat$V2)-10
    treatdata=treat$V2
    controldata=control$V2
    Zn=function(x){Zxi(treatdata,controldata,scaling_fact,x,variance_all)}
    V1<-treat$V1[start:end]
    V2<-sapply(start:end,Zn)
    Ztable<-as.table(cbind(V1,V2))
    N<-paste('Z',files1[i])
    Z[[N]]=Ztable
  }
  Z
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






######
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

