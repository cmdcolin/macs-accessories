########################
# Read wiggle files into memory with assign
files=list.files(path='S96/S96_MACS_wiggle/treat/',pattern="*.fsa.wig.gz")
for (i in files) {
  file<-paste('S96/S96_MACS_wiggle/treat/',i,sep='')
  x<-read.table(file, skip=2)
  assign(i,x)
}

files=list.files(path='HS959/HS959_MACS_wiggle/treat/',pattern="*.fsa.wig.gz")
for (i in files) {
  file<-paste('HS959/HS959_MACS_wiggle/treat/',i,sep='')
  x<-read.table(file, skip=2)
  assign(i,x)
}
files=list.files(path='S96/S96_MACS_wiggle/control/',pattern="*.fsa.wig.gz")
for (i in files) {
  file<-paste('S96/S96_MACS_wiggle/control/',i,sep='')
  x<-read.table(file, skip=2)
  assign(i,x)
}

files=list.files(path='HS959/HS959_MACS_wiggle/control/',pattern="*.fsa.wig.gz")
for (i in files) {
  file<-paste('HS959/HS959_MACS_wiggle/control/',i,sep='')
  x<-read.table(file, skip=2)
  assign(i,x)
}

###################
# Use MACS to find S96 peaks
# macs14 --bw=25 --mfold=4,30 -g 1.2e7 -w -t chip -c control
s96bed=read.table('S96/S96_peaks.bed')
read1=getS96Reads(s96bed)
read2=getHS959Reads(s96bed)



#####
# Overlap of peaks
# intersectBed -a S96_peaks.bed -b HS959_peaks.bed -wa
s96overlap=read.table('S96/S96_overlap.bed')
read3=getS96Reads(s96overlap)
read4=getHS959Reads(s96overlap)

##
# Unique S96 peaks
# subtractBed -a S96_peaks.bed -b S96_overlap.bed
s96unique=read.table('S96/S96_unique.bed')
read5=getS96Reads(s96unique)
read6=getHS959Reads(s96unique)



####
# Total read plot
plot(read1,read2,ylab='HS959 reads',xlab='S96 reads',xlim=c(0,50000),ylim=c(0,30000),pch='*')
points(read3,read4,pch=1,col='red')
points(read5,read6,pch=1,col='green')
title('Total reads S96 peaks')





#############################################
# HS959 peaks
hs959bed=read.table('HS959/HS959_peaks.bed')
read7=getHS959Reads(hs959bed)
read8=getS96Reads(hs959bed)

#####
# Overlap of peaks
# intersectBed -a HS959/HS959_peaks.bed -b S96/S96_peaks.bed -wa
HS959overlap=read.table('HS959/HS959_overlap.bed')
read9=getHS959Reads(HS959overlap)
read10=getS96Reads(HS959overlap)

##
# Unique HS959 
# subtractBed -a HS959/HS959_peaks.bed -b HS959/HS959_overlap.bed
hs959unique=read.table('HS959/HS959_unique.bed')
read11=getHS959Reads(hs959unique)
read12=getS96Reads(hs959unique)

###


plot(read7,read8,ylim=c(0,50000),xlim=c(0,30000),ylab='S96 reads',xlab='HS959 reads',pch='*')
points(read9,read10,pch=1,col='blue')
points(read11,read12,pch=1,col='green')
title('Total reads HS959 peaks')






##########################3
#Avg read plot1
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


###
#Avg read plot2
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
    treat<-get(files1[i])
    control<-get(files2[i])
    ratio_est=as.array(t(control$V2/treat$V2))
    ratio_data=c(ratio_data,t(ratio_est))
  }
  median(ratio_data,TRUE)
}


################
s96scale=estimate_scaling_factor('S96/S96_MACS_wiggle/treat/','S96/S96_MACS_wiggle/control')
hs959scale=estimate_scaling_factor('HS959/HS959_MACS_wiggle/treat/','HS959/HS959_MACS_wiggle/control')












#Functions




# Get total reads in S96 from bedfile
getS96Reads<-function(bedfile) {
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    
    # Reads from S96
    wigfile=paste('S96_treat_afterfiting_',chr,'.wig.gz',sep='');
    wig=get(wigfile);
    peak=wig[wig$V1>start & wig$V1<end,];
    totalreads=apply(t(peak$V2),1,sum);
    reads[i]=totalreads;
  }
  reads
}

############
# Get total reads in HS959 from bedfile
getHS959Reads<-function(bedfile) {
  reads=array();
  for(i in 1:length(bedfile$V1)) {
    chr=bedfile$V1[i];
    start=bedfile$V2[i];
    end=bedfile$V3[i];
    
    #reads for HS959
    wigfile=paste('HS959_treat_afterfiting_',chr,'.wig.gz',sep='');
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

