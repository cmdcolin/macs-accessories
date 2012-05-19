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
    name=name
    )
  nc$getReads = function() {
    loadWiggle(paste(nc$name,'/',nc$name,'_MACS_wiggle/control/',sep=''))
    loadWiggle(paste(nc$name,'/',nc$name,'_MACS_wiggle/treat/',sep=''))
  }
  nc<-list2env(nc)
  class(nc)<-"WiggleClass"
  return(nc)
}




############
# Get total reads in peaks from bedfile
getReads<-function(bedfile,wigpath) {
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    
    wigfile=paste(wigpath,chr,'.wig.gz',sep='')
    wig=get(wigfile);
    peak2=wig[wig$V1>start & wig$V1<end,];
    totalreads2=sum(peak2$V2);
    reads[i]=totalreads2;
  }
  reads
}



#######
# Get avg reads given bedfile and totalreads
getAvgReads<-function(bedfile, totalreads)
{
  newreads=array()
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    newreads[i]=totalreads[i]/(end-start);
  }
  newreads
}


#############
# Get max average reads over window size
getMaxAvgReads<-function(bedfile, wigpath, wsize)
{
  newreads=array()
  chrold=''
  wig=NULL
  for(i in 1:length(bedfile$V1))
  {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    if(as.character(chrold)!=as.character(chr)){
      wigfile=paste(wigpath,chr,'.wig.gz',sep='')
      wig=get(wigfile);
      maxreads=array()
      chrold=chr
      print(wigfile)
    }
    
    
    bstart=start-wsize/10
    bend=end+wsize/10
    b=0
    e=0
    maxreads=array()
    wigpos=wig$V1
    wigpeak=wig$V2
    for(j in seq(bstart,end,by=100)) {
      b=findInterval(j,wigpos)
      e=findInterval(j+wsize,wigpos)
      peak=wigpeak[b:e];
      maxreads[j]=mean(peak,na.rm=TRUE);
    }
    newreads[i]=max(maxreads,na.rm=TRUE);
    print(cat('Found max reads ', newreads[i], ' at ', b, ' ', e,'. Used ', (end-bstart)/100, ' windows'))
  }
  newreads
}







####
# Use all chromosomes for scaling factor
estimate_scaling_factor <- function(path1,path2) {
  files1=list.files(path=path1,pattern="*.fsa.wig.gz")
  files2=list.files(path=path2,pattern="*.fsa.wig.gz")
  ratio_data=array()
  for (i in 1:length(files1)) {
    treat=get(files1[i])
    control=get(files2[i])
    ratio_est=as.array(t(control$V2/treat$V2))
    ratio_data=c(ratio_data,t(ratio_est))
  }
  median(ratio_data,TRUE)
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
    chip_signal=c(chip_signal,as.array(t(treat$V2)))
    control_signal=c(control_signal,as.array(t(control$V2)))
  }
  average_chip=mean(chip_signal,na.rm=TRUE)
  average_control=mean(control_signal,na.rm=TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}


# Sqrt(Aw+Bw/c), w=1
estimate_variance_one<-function(treat,control,scaling_factor,pos) {
  chip_signal=as.array(t(treat$V2[(pos-1):(pos+1)]))
  control_signal=as.array(t(control$V2[(pos-1):(pos+1)]))
  average_chip=mean(chip_signal,na.rm=TRUE)
  average_control=mean(control_signal,na.rm=TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}

# Sqrt(Aw+Bw/c), w=10
estimate_variance_ten<-function(treat,control,scaling_factor,pos) {
  chip_signal=as.array(t(treat$V2[(pos-10):(pos+10)]))
  control_signal=as.array(t(control$V2[(pos-10):(pos+10)]))
  average_chip=mean(chip_signal,na.rm=TRUE)
  average_control=mean(control_signal,na.rm=TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}




Zxi<-function(treat,control,scaling_factor,pos,varianceall) {
  (treat$V2[pos]-control$V2[pos]/scaling_factor)/
      max(estimate_variance_one(treat,control,scaling_factor,pos),
          estimate_variance_ten(treat,control,scaling_factor,pos),
          varianceall)
}



############
# Calculate Z scores over all wiggle files
Z<-function(treat,control,scaling_fact,variance_all) {

  files1=list.files(path=treat,pattern="*.fsa.wig.gz")
  files2=list.files(path=control,pattern="*.fsa.wig.gz")
  Z=list()
  for (i in 1:length(files1)) {
    treat<-get(files1[i])
    control<-get(files2[i])
    start=10
    end=length(treat$V2)-10
    Zn=function(x){Zxi(treat,control,scaling_fact,x,variance_all)}
    Zscore<-sapply(start:end,Zn)
    Ztable<-as.table(cbind(treat$V1,Zscore))
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
    print(typeof(Zchr))
    peak=Zchr[Zchr$V1>start & Zchr$V1<end,];
    reads[i]=mean(peak,na.rm=TRUE)
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
getPeakIndex<-function(bedfile)
{
  index=array()
  for(i in 1:length(bedfile$V1)) {
    index[i]=as.integer(unlist(strsplit(as.character(bedfile$V4[i]),'_'))[3])
  }
  index
}

