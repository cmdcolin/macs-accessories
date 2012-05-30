
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
title('Total S96 peak reads vs HS959 synteny')
plot(r1,r2,ylab='HS959 reads',xlab='S96 reads',pch='*',xlim=c(100,1100),ylim=c(0,600))
points(r1[id1],r2[id1],pch=1,col='pink')
points(r1[id2],r2[id2],pch=1,col='green')
title('Total S96 peak reads vs HS959 synteny (zoom)')






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
plot(r3,r4,ylab='S96 reads',xlab='HS959 reads',pch='*',xlim=c(50,1000),ylim=c(0,1400))
points(r3[id3],r4[id3],pch=1,col='lightblue')
points(r3[id4],r4[id4],pch=1,col='green')
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
s96$estimateScalingFactor()
s96$estimateVarianceAll()
hs959$estimateScalingFactor()
hs959$estimateVarianceAll()

s96=WiggleClass('S96')
Zs96=s96$Z(s96bed)




Rprof()
Zs96=s96$Z(s96bed[1:100,])
Rprof(NULL)
summaryRprof()



##########
# Get Z scores
s96=WiggleClass('S96')
hs959=WiggleClass('HS959')
Zs96=s96$Z(s96bed)
Zhs959=hs959$Z(s96bed)
s96z<-sapply(Zs96,function(x){mean(x[,3])})
hs959z<-sapply(Zhs959,function(x){mean(x[,3])})

plot(s96z,hs959z,pch='*',xlab='Avg S96 peak normdiff',ylab='hs959 syntenic')
points(s96z[id1],hs959z[id1],col='yellow')
points(s96z[id2],hs959z[id2],col='blue')
title('S96 peaks average NormDiff vs HS959 synteny')
plot(s96z,hs959z,pch='*',xlab='Avg S96 peak normdiff',ylab='Avg HS959 syntenic normdiff',xlim=c(1.5,9),ylim=c(-0.5,4))
points(s96z[id1],hs959z[id1],col='yellow')
points(s96z[id2],hs959z[id2],col='blue')
title('S96 peaks avg NormDiff vs HS959 synteny (zoom)')



###########
Zs96dos=s96$Z(hs959bed)
Zhs959dos=hs959$Z(hs959bed)
s96zdos<-sapply(Zs96dos,function(x){mean(x[,3])})
hs959zdos<-sapply(Zhs959dos,function(x){mean(x[,3])})

plot(hs959zdos,s96zdos,pch='*',xlab='Avg S96 peak normdiff',ylab='hs959 syntenic')
points(hs959zdos[id3],s96zdos[id3],col='lightblue')
points(hs959zdos[id4],s96zdos[id4],col='red')
title('S96 peaks average NormDiff vs HS959 synteny')
plot(hs959zdos,s96zdos,pch='*',xlab='Avg S96 peak normdiff',ylab='Avg HS959 syntenic normdiff',xlim=c(1,7),ylim=c(0,8))
points(hs959zdos[id3],s96zdos[id3],col='lightblue')
points(hs959zdos[id4],s96zdos[id4],col='red')
title('HS959 peaks avg NormDiff vs S96 synteny (zoom)')


###########
#Lol getmaxAvgZscore
s96=WiggleClass('S96')
hs959=WiggleClass('HS959')
maxs96<-s96$getMaxAvgZscore(Zs96)
maxhs959<-hs959$getMaxAvgZscore(Zhs959)


plot(maxs96,maxhs959,pch='*',xlab='S96 peak',ylab='HS959 syntenic')
points(maxs96[id1],maxhs959[id1],col='blue')
points(maxs96[id2],maxhs959[id2],col='orange')
title('Max Avg S96 peak normdiff vs HS959 synteny w=100')
plot(maxs96,maxhs959,pch='*',xlab='HS959 peak',ylab='S96 syntenic',xlim=c(1.6,10),ylim=c(-0.2,7))
points(maxs96[id1],maxhs959[id1],col='blue')
points(maxs96[id2],maxhs959[id2],col='orange')
title('Max Avg S96 peak normdiff vs HS959 synteny w=100 (Zoom)')


maxs96dos<-s96$getMaxAvgZscore(Zs96dos)
maxhs959dos<-hs959$getMaxAvgZscore(Zhs959dos)

plot(maxhs959dos,maxs96dos,pch='*',xlab='S96 peak',ylab='HS959 syntenic')
points(maxhs959dos[id3],maxs96dos[id3],col='darkblue')
points(maxhs959dos[id4],maxs96dos[id4],col='orange')
title('Max Avg HS959 peak normdiff vs S96 synteny w=100')
plot(maxhs959dos,maxs96dos,pch='*',xlab='HS959 peak',ylab='S96 syntenic',xlim=c(1.2,9),ylim=c(-0.2,7))
points(maxhs959dos[id3],maxs96dos[id3],col='darkblue')
points(maxhs959dos[id4],maxs96dos[id4],col='orange')
title('Max Avg HS959 peak normdiff vs S96 synteny w=100 (Zoom)')


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
