
s96=WiggleClass('S96')
hs959=WiggleClass('HS959')
s96$loadWiggles()
hs959$loadWiggles()

#######
# Read tables-see notes about command lines
s96bed=read.table('S96/S96_peaks.bed')
s96overlap=read.table('S96/S96_overlap.bed')
s96unique=read.table('S96/S96_unique.bed')
########
# Match unique and overlap indexes
id1=match(s96overlap$V4,s96bed$V4)
id2=match(s96unique$V4,s96bed$V4)

######
# Get total reads
r1=s96$getTotalReads(s96bed)
r2=hs959$getTotalReads(s96bed)


####
# S96Total read plot
plot(r1,r2,ylab='HS959 reads',xlab='S96 reads',pch='*')
points(r1[id1],r2[id1],pch=1,col='pink')
points(r1[id2],r2[id2],pch=1,col='green')
title('Total reads S96 peaks')
plot(r1,r2,ylab='HS959 reads',xlab='S96 reads',pch='*',xlim=c(100,1100),ylim=c(0,600))
points(r1[id1],r2[id1],pch=1,col='pink')
points(r1[id2],r2[id2],pch=1,col='green')
title('Total reads S96 peaks (zoom)')






#############################################
# HS959 peaks
hs959bed=read.table('HS959/HS959_peaks.bed')
HS959overlap=read.table('HS959/HS959_overlap.bed')
hs959unique=read.table('HS959/HS959_unique.bed')
# Match indexes
id3=match(HS959overlap$V4,hs959bed$V4)
id4=match(hs959unique$V4,hs959bed$V4)
r3=hs959$getTotalReads(hs959bed)
r4=s96$getTotalReads(hs959bed)
###
plot(r3,r4,ylab='S96 reads',xlab='HS959 reads',pch='*')
points(r3[id3],r4[id3],pch=1,col='lightblue')
points(r3[id4],r4[id4],pch=1,col='green')
title('Total reads HS959 peaks')
plot(r7,r8,ylab='S96 reads',xlab='HS959 reads',pch='*',xlim=c(50,1000),ylim=c(0,1400))
points(r9,r10,pch=1,col='lightblue')
points(r11,r12,pch=1,col='green')
title('Total reads HS959 peaks (Zoom)')





#####
# Use average reads over peaks
ra1=s96$getAvgReads(s96bed)
ra2=hs959$getAvgReads(s96bed)
# S96 Avg read plot
plot(ra1,ra2,ylab='HS959 reads',xlab='S96 reads',pch='*')
points(ra1[id1],ra2[id1],pch=1,col='lightblue')
points(ra1[id2],ra2[id2],pch=1,col='orange')
title('Avg S96 peak reads vs HS959 synteny')
plot(ra1,ra2,ylab='HS959 reads',xlab='S96 reads',pch='*',xlim=c(9,25),ylim=c(1,15))
points(ra1[id1],ra2[id1],pch=1,col='lightblue')
points(ra1[id2],ra2[id2],pch=1,col='orange')
title('Avg S96 peak reads vs HS959 synteny (zoom)')

#######
ra3=hs959$getAvgReads(hs959bed)
ra4=s96$getAvgReads(hs959bed)
plot(ra3,ra4,ylab='S96 reads',xlab='HS959 reads',pch='*')
points(ra3[id3],ra4[id3],pch=1,col='lightgreen')
points(ra3[id4],ra4[id4],pch=1,col='orange')
title('Average HS959 peak reads vs S96 synteny')
##
plot(ra3,ra4,ylab='S96 reads',xlab='HS959 reads',pch='*',xlim=c(6,16),ylim=c(1,27))
points(ra3[id3],ra4[id3],pch=1,col='lightgreen')
points(ra3[id4],ra4[id4],pch=1,col='orange')
title('Average HS959 peak reads vs S96 synteny (zoom)')






###########
# Use Max avg reads over windows
rma1=s96$getMaxAvgReads(s96bed,100)
rma2=hs959$getMaxAvgReads(s96bed,100)
######################
plot(rma1,rma2,ylab='Max Avg HS959 reads',xlab='Max Avg S96 reads',pch='*')
points(rma1[id1],rma2[id1],pch=1,col='red')
points(rma1[id2],rma2[id2],pch=1,col='green')
title('Max Average reads S96 peaks  w=100')
legend('bottomright', legend=c('shared', 'unique'), fill=c('red', 'green'))
plot(rma1,rma2,ylab='Max Avg HS959 reads',xlab='Max Avg S96 reads',pch='*',xlim=c(10,30),ylim=c(1,20))
points(rma1[id1],rma2[id1],pch=1,col='red')
points(rma1[id2],rma2[id2],pch=1,col='green')
title('Max Avg reads S96 peaks w=100 (Zoom)')
legend('bottomright', legend=c('shared', 'unique'), fill=c('red', 'green'))




rma3=hs959$getMaxAvgReads(hs959bed,100)
rma4=s96$getMaxAvgReads(hs959bed,100)

plot(rma3,rma4,ylab='Max Avg S96 reads',xlab='Max Avg HS959 reads',pch='*')
points(rma3[id3],rma4[id3],pch=1,col='blue')
points(rma3[id4],rma4[id4],pch=1,col='red')
title('Max Average reads HS959 peaks w=100')
legend('bottomright', legend=c('shared', 'unique'), fill=c('blue', 'red'))

#
#Zoom
plot(rma3,rma4,ylab='Max Avg S96 reads',xlab='Max Avg HS959 reads',pch='*',xlim=c(6,30),ylim=c(1,31))
points(rma3[id3],rma4[id3],pch=1,col='blue')
points(rma3[id4],rma4[id4],pch=1,col='red')
title('Max Average reads HS959 peaks w=100 (Zoom)')
legend('bottomright', legend=c('shared', 'unique'), fill=c('blue', 'red'))






################
# Get NormDiff scaling factor, variance
s96scale=s96$estimateScalingFactor()
hs959scale=hs959$estimateScalingFactor()
s96var=estimate_variance_all(s96treat,s96control,s96scale)
hs959var=estimate_variance_all(hs959treat,hs959control,hs959scale)




##########
# Get S96 normdiff
#s96ndpeak<-getAvgNormDiff(s96bed, s96t,s96c, s96scale, s96var)
#s96ndoverlap<-getAvgNormDiff(s96overlap, s96t, s96c, s96scale, s96var)
#s96ndunique<-getAvgNormDiff(s96unique, s96t,s96c, s96scale, s96var)
######
# Get HS959 normdiff
#hs959nd<-getAvgNormDiff(s96bed, hs959t, hs959c, hs959scale, hs959var)
#hs959ndoverlap<-getAvgNormDiff(s96overlap, hs959t,hs959t, hs959scale, hs959var)
#hs959ndunique<-getAvgNormDiff(s96unique, hs959t, hs959c, hs959scale, hs959var)
################
# Plot s96 vs hs959 normdiff scores
#plot(s96ndpeak,hs959nd,pch='*',xlab='S96',ylab='HS959')
#points(s96ndoverlap,hs959ndoverlap,col='red')
#points(s96ndunique,hs959ndunique,col='green')
#title('S96 peaks vs HS959 syntenic NormDiff scores')










##########
# Get Z scores
Zs96=Z(s96treat,s96control,s96scale,s96var,'S96')
Zhs959=Z(hs959treat,hs959control,hs959scale,hs959var,'HS959')
s96z<-getAvgZ(s96bed,Zs96)
s96zoverlap<-getAvgZ(s96overlap,Zs96)
s96zunique<-getAvgZ(s96unique,Zs96)
hs959z<-getAvgZ(s96bed,Zhs959)
hs959zoverlap<-getAvgZ(s96overlap,Zhs959)
hs959zunique<-getAvgZ(s96unique,Zhs959)

plot(s96z,hs959z,pch='*',xlab='Avg S96 peak normdiff',ylab='hs959 syntenic')
points(s96zoverlap,hs959zoverlap,col='yellow')
points(s96zunique,hs959zunique,col='blue')
title('S96 average peak NormDiff')



hs959z<-getAvgZ(hs959bed,Zhs959)
hs959zoverlap<-getAvgZ(HS959overlap,Zhs959)
hs959zunique<-getAvgZ(hs959unique,Zhs959)
s96z<-getAvgZ(s96bed,Zs96)
s96zoverlap<-getAvgZ(s96overlap,Zs96)
s96zunique<-getAvgZ(s96unique,Zs96)
plot(hs959z,s96z,pch='*',xlab='Avg S96 peak normdiff',ylab='S96 syntenic')
points(hs959zoverlap,s96zoverlap,col='yellow')
points(hs959zunique,s96zunique,col='blue')
title('S96 average peak NormDiff')

#s96nd2<-getAvgNormDiff(hs959bed, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
#hs959nd2<-getAvgNormDiff(hs959bed, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
#s96ndoverlap2<-getAvgNormDiff(HS959overlap, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
#hs959ndoverlap2<-getAvgNormDiff(HS959overlap, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
#s96ndunique2<-getAvgNormDiff(hs959unique, 'S96_treat_afterfiting_', 'S96_control_afterfiting_', s96scale, s96var)
#hs959ndunique2<-getAvgNormDiff(hs959unique, 'HS959_treat_afterfiting_', 'HS959_control_afterfiting_', hs959scale, hs959var)
##########
#plot(s96nd2,hs959nd2,pch='*',xlab='S96',ylab='HS959')
#points(s96ndoverlap2,hs959ndoverlap2,col='yellow')
#points(s96ndunique2,hs959ndunique2,col='blue')
#title('HS959 peaks NormDiff')



hs959z<-getMaxAvgZ(hs959bed,Zhs959,100)
hs959zoverlap<-getMaxAvgZ(HS959overlap,Zhs959,100)
hs959zunique<-getMaxAvgZ(hs959unique,Zhs959,100)
s96z<-getMaxAvgZ(hs959bed,Zs96,100)
s96zoverlap<-getMaxAvgZ(HS959overlap,Zs96,100)
s96zunique<-getMaxAvgZ(hs959unique,Zs96,100)
plot(hs959z,s96z,pch='*',xlab='HS959 peak',ylab='S96 syntenic')
points(hs959zoverlap,s96zoverlap,col='blue')
points(hs959zunique,s96zunique,col='orange')
title('Max Avg HS959 peak normdiff w=100')
plot(hs959z,s96z,pch='*',xlab='HS959 peak',ylab='S96 syntenic',ylim=c(1,13),xlim=c(3.3,11))
points(hs959zoverlap,s96zoverlap,col='blue')
points(hs959zunique,s96zunique,col='orange')
title('Max Avg HS959 peak normdiff w=100 (Zoom)')




s96z2<-getMaxAvgZ(s96bed,Zs96,100)
s96zoverlap2<-getMaxAvgZ(s96overlap,Zs96,100)
s96zunique2<-getMaxAvgZ(s96unique,Zs96,100)
hs959z2<-getMaxAvgZ(s96bed,Zhs959,100)
hs959zoverlap2<-getMaxAvgZ(s96overlap,Zhs959,100)
hs959zunique2<-getMaxAvgZ(s96unique,Zhs959,100)
plot(s96z2,hs959z2,pch='*',xlab='S96 peak',ylab='HS959 syntenic')
points(s96zoverlap2,hs959zoverlap2,col='green')
points(s96zunique2,hs959zunique2,col='brown')
title('Max Avg S96 peak normdiff w=100')

plot(s96z2,hs959z2,pch='*',xlab='S96 peak',ylab='HS959 syntenic',ylim=c(0,10),xlim=c(3.3,13))
points(s96zoverlap2,hs959zoverlap2,col='green')
points(s96zunique2,hs959zunique2,col='brown')
title('Max Avg S96 peak normdiff w=100 (Zoom)')
##########3

################
# Sort and rank
# $ intersectBed -a S96/S96_peaks.bed -b HS959/HS959_peaks.bed -f 1.0 >> fullIntersect.bed

fullintersect=read.table('fullIntersect3.bed')
fullintersectZ<-getMaxAvgZ(fullintersect,Zs96,100)
sorts96Z=sort(s96z2)
plot(1:length(sorts96Z),sorts96Z,type='l',xlab='Rank', ylab='max avg NormDiff score')
matchpos=match(fullintersectZ,sorts96Z)
points(matchpos,fullintersectZ)
title('Sorted Max Avg. NormDiff (S96 peaks) w=100')

# $ intersectBed -a HS959/HS959_peaks.bed -b S96/S96_peaks.bed -f 1.0 >> fullIntersect2.bed
fullintersect2=read.table('fullIntersect2.bed')
fullintersectZ2<-getMaxAvgZ(fullintersect2,Zhs959,100)
sorths959Z=sort(hs959z)
plot(1:length(sorths959Z),sorths959Z,type='l',xlab='Rank', ylab='max avg NormDiff score')
matchpos=match(fullintersectZ2,sorths959Z)
points(matchpos,fullintersectZ2,pch='.',col='red')
title('Sorted Max Avg. NormDiff (HS959 peaks) w=100')
                   
                   
#####################
# Whole genome
Zs96genome=getAvgZ_WholeGenome(Zs96,100)
Zhs959genome=getAvgZ_WholeGenome(Zhs959,100)
Zs96genome=Zs96genome[-(length(Zs96genome):(length(Zhs959genome)+1))]
plot(Zs96genome,Zhs959genome,pch='.',xlim=c(-4,6),ylim=c(-3.5,6))


#########
# MANORM

M=log2(s96z2)-log2(hs959z2)
A=1/2*(log2(s96z2)+log2(hs959z2))
plot(A,M,ylim=c(-3,3),xlim=c(1,5.1),xlab='Normdiff Intensity', ylab='Normdiff ratio',pch='*')
M=log2(s96zoverlap2)-log2(hs959zoverlap2)
A=1/2*(log2(s96zoverlap2)+log2(hs959zoverlap2))
points(A,M,col='blue')
M=log2(s96zunique2)-log2(hs959zunique2)
A=1/2*(log2(s96zunique2)+log2(hs959zunique2))
points(A,M,col='red')
title('MA plot HS959 peaks (max avg NormDiff w=100)')
legend('bottomright', legend=c('shared', 'unique'), fill=c('blue', 'red'))



M=log2(hs959z)-log2(s96z)
A=1/2*(log2(s96z)+log2(hs959z))
plot(A,M,ylim=c(-3,3),xlim=c(1,5.1),xlab='Log product S96,HS959', ylab='Log ratio S96,HS959',pch='*')
M=log2(hs959zoverlap)-log2(s96zoverlap)
A=1/2*(log2(s96zoverlap)+log2(hs959zoverlap))
points(A,M,col='green')
M=log2(hs959zunique)-log2(s96zunique)
A=1/2*(log2(s96zunique)+log2(hs959zunique))
points(A,M,col='red')
title('MA plot S96 peaks (max avg NormDiff w=100)')
legend('bottomright', legend=c('shared', 'unique'), fill=c('green', 'red'))


########

M=log2(s96z)-log2(hs959z)
A=1/2*(log2(s96z)+log2(hs959z))
plot(A,M,ylim=c(-3,3),xlim=c(1,5.1),xlab='Log product S96,HS959', ylab='Log ratio S96,HS959',pch='*')
M=log2(s96zoverlap)-log2(hs959zoverlap)
A=1/2*(log2(s96zoverlap)+log2(hs959zoverlap))
points(A,M,col='blue')
M=log2(s96zunique)-log2(hs959zunique)
A=1/2*(log2(s96zunique)+log2(hs959zunique))
points(A,M,col='green')

M=log2(s96z)-log2(hs959z)
A=1/2*(log2(s96z)+log2(hs959z))
points(A,M,pch='*')
title('MA plot (max avg NormDiff w=100)')
legend('bottomright', legend=c('shared', 'S96', 'HS959'), fill=c('blue', 'red','green'))


s96sort=sort(s96nd)
plot(1:length(s96sort),s96sort,pch='*',xlab='Rank',ylab='Avg NormDiff')
title('S96 Normdiff Sort')


################
getPeakIndex(s96overlap)
