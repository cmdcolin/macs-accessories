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
read=array()
read2=array()
for(i in 1:length(s96bed$V1))
{
  chr=s96bed$V1[i];
  start=s96bed$V2[i];
  end=s96bed$V3[i];
                
  # Reads from S96
  wigfile=paste ('S96_treat_afterfiting_' , chr, '.wig.gz', sep='')
  wig=get(wigfile)
  peak=wig[wig$V1>start & wig$V1<end,];
  totalreads=apply(t(peak$V2),1,sum) ;
  read[i]=totalreads
                
  # reads for HS959
  wigfile=paste('HS959_treat_afterfiting_',chr, '.wig.gz',sep='')
  wig=get(wigfile)
  peak2=wig[wig$V1>start & wig$V1<end,];
  totalreads2=apply(t(peak2$V2),1,sum);
  read2[i]=totalreads2
} 


plot(read,read2,ylim=c(0,30000),xlim=c(0,50000),ylab='HS959 reads',xlab='S96 reads',pch='*')
title('S96 peaks')




# Overlap of peaks in both HS959 and S96
# intersectBed -a S96_peaks.bed -b HS959_peaks.bed -wa >> S96_overlap.bed
s96bed=read.table('S96_overlap.bed')
read3=array()
read4=array()
for(i in 1:length(s96bed$V1))
{
  chr=s96bed$V1[i];
  start=s96bed$V2[i];
  end=s96bed$V3[i];
  
  # Reads from S96
  wigfile=paste('S96_treat_afterfiting_',chr,'.wig.gz',sep='')
  wig=get(wigfile)
  peak=wig[wig$V1>start & wig$V1<end,];
  totalreads=apply(t(peak$V2),1,sum)
  read3[i]=totalreads
  
  # reads for HS959
  wigfile=paste('HS959_treat_afterfiting_',chr,'.wig.gz',sep='')
  wig=get(wigfile)
  peak2=wig[wig$V1>start & wig$V1<end,];
  totalreads2=apply(t(peak2$V2),1,sum);
  read4[i]=totalreads2
}
points(read3,read4,pch=1,col='red')

##
# Unique S96 peaks
# subtractBed -a S96_peaks.bed -b S96_overlap.bed >> S96_unique.bed
s96bed=read.table('S96_unique.bed')
read5=array()
read6=array()
for(i in 1:length(s96bed$V1)) {
  chr=s96bed$V1[i];
  start=s96bed$V2[i];
  end=s96bed$V3[i];
  
  # Reads from S96
  wigfile=paste('S96_treat_afterfiting_',chr,'.wig.gz',sep='')
  wig=get(wigfile)
  peak=wig[wig$V1>start & wig$V1<end,];
  totalreads=apply(t(peak$V2),1,sum);
  read5[i]=totalreads
  
  #reads for HS959
  wigfile=paste('HS959_treat_afterfiting_',chr,'.wig.gz',sep='')
  wig=get(wigfile)
  peak2=wig[wig$V1>start & wig$V1<end,];
  totalreads2=apply(t(peak2$V2),1,sum);
  read6[i]=totalreads2
}
points(read5,read6,pch=1,col='green')


