ret1<-read.table('data/S96vsHS959-peaks-annotate.txt',sep='\t',header=TRUE)
ret1chrmod<-read.table('data/S96vsHS959-Diff_peaks.bed.bak',sep='\t',header=FALSE)
ret2<-read.table('data/s96-hs959-diff-genes.csv',sep=',')

fd<-featureData(gset)
pdiddy<-pData(fd)


genelist<-as.character(pData(fd)[,5])
geneIntersect<-intersect(genelist,as.character(ret2$V7))
geneIntersect<-intersect(genelist,as.character(ret1$Nearest.PromoterID))

#subset gene table
pdiddyret<-pdiddy[pdiddy[,5]%in%geneIntersect,]
pdiddyret2<-pdiddyret[pdiddyret[,5]!="",]
exprret<-exprs(gset)[pdiddy[,5]%in%geneIntersect,]


gimme<-merge(pdiddyret2,ret1,by.x="ORF",by.y="Nearest.PromoterID")

pdiddyret2ret2<-join(pdiddyret2,ret2,by="ORF")
plot(log(exprret[,1]/exprret[,2]))


gsetexprs<-exprs(gset)[gimme$ID,]
plot(log(gsetexprs[1,]/gsetexprs[,2]),gimme$Peak.Score,pch=16,cex=0.6,xlab="Log ratio (Gene expression)",ylab="Nearest ChIP-seq peak -log10(pvalue)")


