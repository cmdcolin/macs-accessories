ret1<-read.table('data/S96vsHS959-peaks-annotate.txt',sep='\t',header=TRUE)
ret1chrmod<-read.table('data/S96vsHS959-Diff_peaks.bed.bak',sep='\t',header=FALSE)

ret2<-read.table('data/HS959vsS96-Diff_peaks-annotate.txt',sep='\t',header=TRUE)
ret2chrmod<-read.table('data/HS959vsS96-Diff_peaks.bed.bak',sep='\t',header=FALSE)



fd<-featureData(gset)
pdiddy<-pData(fd)


genelist<-as.character(pData(fd)[,5])
geneIntersect<-intersect(genelist,as.character(ret1$Nearest.PromoterID))
geneIntersect2<-intersect(genelist,as.character(ret2$Nearest.PromoterID))
#subset gene table
pdiddyret<-pdiddy[pdiddy[,5]%in%geneIntersect,]
pdiddyret<-pdiddyret[pdiddyret[,5]!="",]

pdiddyret2<-pdiddy[pdiddy[,5]%in%geneIntersect2,]
pdiddyret2<-pdiddyret2[pdiddyret2[,5]!="",]


#merge peak table with gene table by the nearest orfs
gimme<-merge(pdiddyret,ret1,by.x="ORF",by.y="Nearest.PromoterID")
gimme<-gimme[!duplicated(gimme$ORF),]


gimme2<-merge(pdiddyret2,ret2,by.x="ORF",by.y="Nearest.PromoterID")
gimme2<-gimme2[!duplicated(gimme2$ORF),]



gsetexprs<-exprs(gset)[gimme$ID,]
gsetexprs2<-exprs(gset)[gimme2$ID,]


plot(log(gsetexprs[,1]/gsetexprs[,2]),gimme$Peak.Score,pch=16,cex=0.6,xlab="Log ratio (Gene expression)",ylab="Nearest ChIP-seq peak score -log(pval)",log="y")
points(log(gsetexprs2[,1]/gsetexprs2[,2]),gimme2$Peak.Score,pch=16,cex=0.6,col=2)

names(ret1chrmod)<-c('chr','start','end','V4','score')
names(ret2chrmod)<-c('chr','start','end','V4','score')
gimmePeakScore<-merge(gimme,ret1chrmod,by.x="PeakID..cmd.S96vsHS959.Diff_peaks.bed.sacCer3.",by.y="V4")
gimmePeakScore2<-merge(gimme2,ret2chrmod,by.x="PeakID..cmd.data.HS959vsS96.Diff_peaks.bed.sacCer3.",by.y="V4")

##SCALING
wiggleTableScale<-normalizeBetweenArrays(as.matrix(wiggleTable[,c(3,4,5,6,7,8,9,10)]),method="quantile")
wiggleTableScale<-cbind(wiggleTable[,c(1,2)],wiggleTableScale)

genepeakscore<-getPeakScores(gimmePeakScore,wiggleTableScale,do.rbind=FALSE)

genepeakscore<-do.call(rbind,lapply(genepeakscore,function(x)colSums(x[,c(4,6,8,10)])))



genepeakscore2<-getPeakScores(gimmePeakScore2,wiggleTableScale,do.rbind=FALSE)

genepeakscore2<-do.call(rbind,lapply(genepeakscore2,function(x)colSums(x[,c(4,6,8,10)])))


l1<-log2(genepeakscore[,3]/genepeakscore[,1])
l2<-log2(genepeakscore[,2]/genepeakscore[,4])

l3<-log2(genepeakscore2[,3]/genepeakscore2[,1])
l4<-log2(genepeakscore2[,2]/genepeakscore2[,4])



plot(log2(gsetexprs[1,]/gsetexprs[,2]),l1,pch=16,cex=0.9,xlab="Log ratio (Gene expression)",ylab="Log ratio (ChIP-seq)",xlim=c(-2,2),ylim=c(-3,3))
points(log2(gsetexprs2[,1]/gsetexprs2[,2]),l3,pch=16,cex=0.9,col=2)

legend("bottomright",fill=c(1,2),legend=c("S96 differential peaks","HS959 differential peaks"))
title("Differential peaks compared with differences in gene expression")

ret<-head(sort(log2(gsetexprs[,1]/gsetexprs[,2])),3)
gimme[gimme$ID%in%names(ret),]

ret2<-head(sort(log2(gsetexprs2[,1]/gsetexprs2[,2])),3)
gimme2[gimme2$ID%in%names(ret2),]



###BOXPLOT
mytitle<-"S96 replicates vs HS959 replicates"
par(mfrow=c(1,2))
palette(brewer.pal(4,"RdBu"))
fl<-as.factor(names(wiggleTable)[c(4,6,8,10)])
boxplot(wiggleTable[,c(4,6,8,10)], boxwex=0.5, notch=FALSE, main=mytitle, outline=FALSE, las=2, col=fl)
boxplot(wiggleTableScale[,c(4,6,8,10)], boxwex=0.6, notch=T, main=mytitle, outline=FALSE, las=2, col=fl)
par(mfrow=c(1,1))
palette("default")