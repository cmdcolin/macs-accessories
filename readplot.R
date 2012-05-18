########################
# Read wiggle files into memory with assign
loadWiggle<-function(wigpath) {
  files=list.files(path=wigpath,pattern="*.fsa.wig.gz")
  for (i in files) {
    file<-paste(wigpath,i,sep='')
    x<-read.table(file, skip=2)
    assign(i,x,inherits=TRUE)
  }
}
loadWiggle('S96/S96_MACS_wiggle/treat/')
loadWiggle('S96/S96_MACS_wiggle/control/')
loadWiggle('HS959/HS959_MACS_wiggle/treat/')
loadWiggle('HS959/HS959_MACS_wiggle/control/')

###################
# Use MACS to find S96 peaks
# macs14 --bw=25 --mfold=4,30 -g 1.2e7 -w -t chip -c control



#Functions
# Get total reads in S96 from bedfile


############
# Get total reads in HS959 from bedfile
getReads<-function(bedfile,wigpath) {
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    
    #reads for HS959
    wigfile=paste(wigpath,chr,'.wig.gz',sep='');
    wig=get(wigfile);
    peak2=wig[wig$V1>start & wig$V1<end,];
    totalreads2=apply(t(peak2$V2),1,sum);
    reads[i]=totalreads2;
  }
  reads
}


s96bed=read.table('S96/S96_peaks.bed')
read1=getReads(s96bed,'S96_treat_afterfiting_')
read2=getReads(s96bed,'HS959_treat_afterfiting_')



#####
# Overlap of peaks
# intersectBed -a S96_peaks.bed -b HS959_peaks.bed -wa
s96overlap=read.table('S96/S96_overlap.bed')
read3=getReads(s96overlap,'S96_treat_afterfiting_')
read4=getReads(s96overlap,'HS959_treat_afterfiting_')

##
# Unique S96 peaks
# subtractBed -a S96_peaks.bed -b S96_overlap.bed
s96unique=read.table('S96/S96_unique.bed')
read5=getReads(s96unique,'S96_treat_afterfiting_')
read6=getReads(s96unique,'HS959_treat_afterfiting_')



####
# Total read plot
plot(read1,read2,ylab='HS959 reads',xlab='S96 reads',xlim=c(0,50000),ylim=c(0,30000),pch='*')
points(read3,read4,pch=1,col='red')
points(read5,read6,pch=1,col='green')
title('Total reads S96 peaks')





#############################################
# HS959 peaks
hs959bed=read.table('HS959/HS959_peaks.bed')
read7=getReads(hs959bed,'HS959_treat_afterfiting_')
read8=getReads(hs959bed,'S96_treat_afterfiting_')

#####
# Overlap of peaks
# intersectBed -a HS959/HS959_peaks.bed -b S96/S96_peaks.bed -wa
HS959overlap=read.table('HS959/HS959_overlap.bed')
read9=getReads(HS959overlap,'HS959_treat_afterfiting_')
read10=getReads(HS959overlap,'S96_treat_afterfiting_')

##
# Unique HS959 
# subtractBed -a HS959/HS959_peaks.bed -b HS959/HS959_overlap.bed
hs959unique=read.table('HS959/HS959_unique.bed')
read11=getReads(hs959unique,'HS959_treat_afterfiting_')
read12=getReads(hs959unique,'S96_treat_afterfiting_')

###
plot(read7,read8,ylim=c(0,50000),xlim=c(0,30000),ylab='S96 reads',xlab='HS959 reads',pch='*')
points(read9,read10,pch=1,col='blue')
points(read11,read12,pch=1,col='green')
title('Total reads HS959 peaks')




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

avgread7=getAvgReads(hs959bed,read7)
avgread8=getAvgReads(hs959bed,read8)
avgread9=getAvgReads(HS959overlap,read9)
avgread10=getAvgReads(HS959overlap,read10)
avgread11=getAvgReads(hs959unique,read11)
avgread12=getAvgReads(hs959unique,read12)
plot(avgread7,avgread8,ylab='S96 reads',xlab='HS959 reads',pch='*')
points(avgread9,avgread10,pch=1,col='blue')
points(avgread11,avgread12,pch=1,col='green')
title('Average reads HS959 peaks')




avgread1=getAvgReads(s96bed,read1)
avgread2=getAvgReads(s96bed,read2)
avgread3=getAvgReads(s96overlap,read3)
avgread4=getAvgReads(s96overlap,read4)
avgread5=getAvgReads(s96unique,read5)
avgread6=getAvgReads(s96unique,read6)
plot(avgread1,avgread2,ylab='HS959 reads',xlab='S96 reads',pch='*')
points(avgread3,avgread4,pch=1,col='red')
points(avgread5,avgread6,pch=1,col='green')
title('Average reads S96 peaks')





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




################
s96scale=estimate_scaling_factor('S96/S96_MACS_wiggle/treat/','S96/S96_MACS_wiggle/control')
hs959scale=estimate_scaling_factor('HS959/HS959_MACS_wiggle/treat/','HS959/HS959_MACS_wiggle/control')
s96var=estimate_variance_all('S96/S96_MACS_wiggle/treat/','S96/S96_MACS_wiggle/control',s96scale)
hs959var=estimate_variance_all('HS959/HS959_MACS_wiggle/treat/','HS959/HS959_MACS_wiggle/control',hs959scale)

















AvgNormDiff<-function(bedfile, path1, path2, scaling_factor, variance) {
  
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
    #reads[i]=reads[i]/(end-start)
  }
  reads
}
s96nd<-AvgNormDiff(s96bed, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959nd<-AvgNormDiff(s96bed, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
s96ndoverlap<-AvgNormDiff(s96overlap, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959ndoverlap<-AvgNormDiff(s96overlap, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
s96ndunique<-AvgNormDiff(s96unique, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959ndunique<-AvgNormDiff(s96unique, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)

plot(hs959nd,s96nd,pch='*',xlab='HS959',ylab='S96')
points(hs959ndoverlap,s96ndoverlap,col='red')
points(hs959ndunique,s96ndunique,col='green')
title('S96 peaks NormDiff')


s96nd2<-AvgNormDiff(hs959bed, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959nd2<-AvgNormDiff(hs959bed, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
s96ndoverlap2<-AvgNormDiff(HS959overlap, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959ndoverlap2<-AvgNormDiff(HS959overlap, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
s96ndunique2<-AvgNormDiff(hs959unique, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959ndunique2<-AvgNormDiff(hs959unique, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)




plot(s96nd2,hs959nd2,pch='*',xlab='S96',ylab='HS959')
points(s96ndoverlap2,hs959ndoverlap2,col='yellow')
points(s96ndunique2,hs959ndunique2,col='blue')
title('HS959 peaks NormDiff')




s96overlapnd<-AvgNormDiff(s96overlap, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
s96uniquend<-AvgNormDiff(s96unique, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
s96sort=sort(s96nd)
hs959=sort(hs959nd)
plot(1:length(s96nd),s96nd,pch='*',xlab='Rank',ylab='Avg NormDiff')

title('S96 Normdiff Sort')
plot(1:length(s96sort),s96sort,pch='*',xlab='Rank',ylab='Avg NormDiff')
title('HS959 NormDiff')

getPeakIndex(s96overlap)


getPeakIndex<-function(bedfile)
{
  index=array()
  for(i in 1:length(bedfile$V1)) {
    index[i]=as.integer(substring(bedfile$V4[i],11,length(bedfile$V4[i])))
    print(substring(bedfile$V4[i],11,length(bedfile$V4[i])))
  }
  index
}

