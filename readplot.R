##
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
  print (i)
}

####
# Use MACS to find S96 peaks
# macs14 --bw=25 --mfold=4,30 -g 1.2e7 -w -t chip -c control
s96bed=read.table('S96/S96_peaks.bed')
read1=getS96Reads(s96bed)
read2=getHS959Reads(s96bed)
plot(read1,read2,ylim=c(0,30000),xlim=c(0,50000),ylab='HS959 reads',xlab='S96 reads',pch='*')
title('S96 peaks')



#####
# Overlap of peaks in both HS959 and S96
# intersectBed -a S96_peaks.bed -b HS959_peaks.bed -wa >> S96_overlap.bed
s96overlap=read.table('S96_overlap.bed')
read3=getS96Reads(s96overlap)
read4=getHS959Reads(s96overlap)
points(read3,read4,pch=1,col='red')

##
# Unique S96 peaks
# subtractBed -a S96_peaks.bed -b S96_overlap.bed >> S96_unique.bed
s96unique=read.table('S96_unique.bed')
read5=getS96Reads(s96unique)
read6=getHS959Reads(s96unique)
points(read5,read6,pch=1,col='green')



############
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
