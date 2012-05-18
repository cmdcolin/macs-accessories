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
    
    #reads for HS959
    wigfile=paste(wigpath,chr,'.wig.gz',sep='')
    wig=get(wigfile);
    peak2=wig[wig$V1>start & wig$V1<end,];
    totalreads2=apply(t(peak2$V2),1,sum);
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
  average_chip=mean(chip_signal,0,TRUE)
  average_control=mean(control_signal,0,TRUE)
  sqrt(average_chip+average_control/scaling_factor^2)
}




getAvgNormDiff<-function(bedfile, path1, path2, scaling_factor, variance) {
  
  # Get total reads in S96 from bedfile
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
  
    Z=array()
    treat_peak=as.array(treat_peak$V2)
    control_peak=as.array(control_peak$V2)
    
    for(j in 1:(end-start)) {
      Z[j]= (treat_peak[j]-control_peak[j]/scaling_factor)/sqrt(variance)
    }
    reads[i]=mean(Z,na.rm=TRUE);
  }
  reads
}




######
# Get peak indexes from bedfile column
getPeakIndex<-function(bedfile)
{
  index=array()
  for(i in 1:length(bedfile$V1)) {
    index[i]=as.integer(substring(bedfile$V4[i],11,length(bedfile$V4[i])))
    print(substring(bedfile$V4[i],11,length(bedfile$V4[i])))
  }
  index
}

