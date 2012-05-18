###################
# Use MACS to generate wiggle files and find peaks
#
# >macs14 --bw=25 --mfold=4,30 -g 1.2e7 -w -t chip -c control

s96treat<-'S96/S96_MACS_wiggle/treat/'
s96control<-'S96/S96_MACS_wiggle/control/'
hs959treat<-'HS959/HS959_MACS_wiggle/treat/'
hs959control<-'HS959/HS959_MACS_wiggle/control/'
s96tafter<-'S96_treat_afterfiting_'

loadWiggle(s96treat)
loadWiggle(s96control)
loadWiggle(hs959treat)
loadWiggle(hs959control)




s96bed=read.table('S96/S96_peaks.bed')
read1=getReads(s96bed, 'S96_treat_afterfiting_')
read2=getReads(s96bed, 'HS959_treat_afterfiting_')



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






#####

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




################
# Get NormDiff scaling factor, variance
s96scale=estimate_scaling_factor(s96treat,s96control)
hs959scale=estimate_scaling_factor(hs959treat,hs959control)
s96var=estimate_variance_all(s96treat,s96control,s96scale)
hs959var=estimate_variance_all(hs959treat,hs959control,hs959scale)

Zs96=Z(s96treat,s96control,s96scale,s96var)
Zhs959=Z(hs959treat,hs959control,hs959scale,hs959var)



##########
# Get S96 normdiff
s96ndpeak<-getAvgNormDiff(s96bed, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
s96ndoverlap<-getAvgNormDiff(s96overlap, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
s96ndunique<-getAvgNormDiff(s96unique, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
######
# Get HS959 normdiff
hs959nd<-getAvgNormDiff(s96bed, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
hs959ndoverlap<-getAvgNormDiff(s96overlap, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
hs959ndunique<-getAvgNormDiff(s96unique, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)




################
# Plot s96 vs hs959 normdiff scores
plot(hs959nd,s96ndpeak,pch='*',xlab='HS959',ylab='S96')
points(hs959ndoverlap,s96ndoverlap,col='red')
points(hs959ndunique,s96ndunique,col='green')
title('S96 peaks vs HS959 syntenic NormDiff scores')


s96nd2<-getAvgNormDiff(hs959bed, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959nd2<-getAvgNormDiff(hs959bed, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
s96ndoverlap2<-getAvgNormDiff(HS959overlap, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959ndoverlap2<-getAvgNormDiff(HS959overlap, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
s96ndunique2<-getAvgNormDiff(hs959unique, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
hs959ndunique2<-getAvgNormDiff(hs959unique, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)




plot(s96nd2,hs959nd2,pch='*',xlab='S96',ylab='HS959')
points(s96ndoverlap2,hs959ndoverlap2,col='yellow')
points(s96ndunique2,hs959ndunique2,col='blue')
title('HS959 peaks NormDiff')





##########3
s96sort=sort(s96nd)
plot(1:length(s96sort),s96sort,pch='*',xlab='Rank',ylab='Avg NormDiff')
title('S96 Normdiff Sort')


################
getPeakIndex(s96overlap)
