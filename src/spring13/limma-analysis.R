#LIMMA-ANALYSIS.R


## QQPlot example for log data
y=log(wiggleTable[,4])
qqnorm(y); qqline(y, col = 2)

##ddply example
ddply(wiggleTable,.(chr),summarize,)



ddply(wiggleTable,.(chr),summarize,outm=mean(V4))
ddply(wiggleTable,.(chr),head)


# plot(ret2$V1,ret2$V2)
# plot(ret2$V1,ret2$V2,type='l')
# lines(ret2$V1,ret2$V3)
# lines(ret2$V1,ret2$V3,col=2)
# plot(ret2$V1,ret2$V2,type='l',xlim=c(1,10000))
# plot(ret2$V1,ret2$V2,type='l',xlim=c(1,100000))
# plot(ret2$V1,ret2$V2,type='l',xlim=c(1,100000))
# lines(ret2$V1,ret2$V3,col=2)
# plot(ret2$V1,ret2$V2,type='l',xlim=c(1,100000),ylim=c(15,20))
# plot(ret2$V1,log(ret2$V2),type='l',xlim=c(1,100000),ylim=c(15,20))
# plot(ret2$V1,log(ret2$V2),type='l',xlim=c(1,100000))
# lines(ret2$V1,log(ret2$V3),col=2)
# plot(ret2$V1,log(ret2$V2),type='l',xlim=c(1,100000),ylim(0,6))
# plot(ret2$V1,log(ret2$V2),type='l',xlim=c(1,100000),ylim=c(0,6))
# lines(ret2$V1,log(ret2$V3),col=2)
# qqplot(log(ret2$V3))
# ?qqplot
# qqplot(log(ret2$V3),rnorm(100))
# qqplot(log(ret2$V3),rnorm(1000))
# qqplot(log(ret2$V3),rnorm(10000))
# qqplot(log(ret2$V3),rnorm(100000))
# qqplot(log(ret2$V3),rnorm(100000))
# qqplot(log(ret2$V2),rnorm(100000))
# qqplot(log(ret2$V3),rnorm(100000))
# qqline(rnorm(100000))
# qqline(log(rnorm(100000)))

getVennDiagram<-function(bed1,bed2) {
  r1<-intersectBed(bed2,bed1)
  r2<-uniqueBed(bed2,bed1)
  r3<-uniqueBed(bed1,bed2)
  s1<-c(r1$name, paste(r2$name,"unique"))
  s2<-c(r1$name,paste(r3$name,"unique2"))
  
  #r3<-unlist(apply(cbind(r1,r2),1,function(x) !(r1||r2)))
  #print(head(r3))
  #r3<-uniqueBedLimma(bed1,bed2)
  ret<-list(s1,s2)
  venn(ret)
  
}
bed1<-loadBed('s96rep1-high_peaks.bed')
bed2<-loadBed('hs959rep1-new_peaks.bed')
ret<-getVennDiagram(bed1,bed2)

rr<-read.table('GSE19635_HS_peaks.txt',header=TRUE)
rr2<-read.table('GSE19635_s96a_peaks.txt',header=TRUE)
rr$name<-paste0('MACS_PEAK_',1:nrow(rr))
rr2$name<-paste0('MACS_PEAK_',1:nrow(rr2))
ret<-getVennDiagram(rr,rr2)




getVennDiagram(rr,bed2)

scaler<-mean(apply(wiggleTable[,3:ncol(wiggleTable)],2,median))
for(i in 3:ncol(wiggleTable)) {
  wiggleScale[,i]<-wiggleTable[,i]/(mean(wiggleTable[,i])/scaler)
}

ret2<-getPeakScores(bed2,wiggleScale,FALSE)

mytests<-lapply(ret2,function(rtest) { 
  if(nrow(rtest)==0) {
    0
  } else {
    t.test(log(rtest[,4]/rtest[,6]),mu=0)
  }
})
pvals<-unlist(lapply(mytests,function(x)x['p.value']))
hist(pvals)




makeComparisonPlot(loadBed('GSE19635_s96a_peaks.txt',header=TRUE),loadBed('GSE19635_HS_peaks.txt',header=TRUE),wiggleTable,4,8,'Comparison of raw read scores for S96 vs HS959','S96', 'HS959', c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Dark2"),FALSE,TRUE)


out<-normalizeBetweenArrays(as.matrix(wiggleTable[,c(4,6,8,10,12,14,16,18)]),method="scale")
noNormWigTab<-wiggleTable
wiggleTable<-out
rm(out)
M=log2(wiggleTable[,c('V4','V6')]/wiggleTable[,c('V8','V10')])
#M=log2(wiggleTable[,c(4,6,8,10)])
fit<-lmFit(M)
fit<-eBayes(fit)
###############################################
set.seed(2004); invisible(runif(100))
M <- matrix(rnorm(100*6,sd=0.3),100,6)
M[1,] <- M[1,] + 1
fit <- lmFit(M)

#  Ordinary t-statistic
par(mfrow=c(1,2))
ordinary.t <- fit$coef / fit$stdev.unscaled / fit$sigma
qqt(ordinary.t,df=fit$df.residual,main="Ordinary t")
abline(0,1)

#  Moderated t-statistic
eb <- eBayes(fit)
qqt(eb$t,df=eb$df.prior+eb$df.residual,main="Moderated t")
abline(0,1,col=2)
#  Points off the line may be differentially expressed
par(mfrow=c(1,1))
#####################################################
#########################

x<-topTable(fit,number=5000)
volcanoplot(eb,coef=1)
points(x$logFC,x$B,col="#FF000077",pch=16,cex=0.35)
title('Log Fold change vs Log odds (S96 vs HS959)')
title('Log difference vs Log odds S96vsHS959')



x2<-wiggleTable[x$ID,]
write.table(cbind(wiggleTable[x$ID,1],wiggleTable[x$ID,2],wiggleTable[x$ID,2]+1),col.names=FALSE,row.names=FALSE,sep='\t')
plot(wiggleTable[,c(4,8)],pch=16,cex=.25,col=rgb(1,0,0,0.5),xlab)
points(x2[,c(4,8)],pch=16,cex=.25,col=rgb(0,0,0,0.5))









out<-normalizeBetweenArrays(as.matrix(wiggleTable[,c(4,6,8,10,12,14,16,18)]),method="scale")


# K-Means Clustering with 5 clusters
fit <- kmeans(mydata, 5)

# Cluster Plot against 1st 2 principal components

# vary parameters for most readable graph
library(cluster)
clusplot(mydata, fit$cluster, color=TRUE, shade=TRUE,
         labels=2, lines=0)

y<-predFC(D2,prior.count.total=2*ncol(D2))
heatmap.2(y,Rowv=NA,Colv=NA,trace="none")


et <- exactTest(D2)
topTags(et, n=20)



makeComparisonPlotHelp2(wiggleTable[x$ID,],loadBed("s96vshs959overlap.bed")wiggleTable,4,8,'Differential peaks highlighted with t-test','S96', 'HS959', c("Overlap","Differential"),brewer.pal(3,"Dark2"),FALSE,FALSE)


makeComparisonPlotHelp2(wiggleTable[x$ID,],loadBed("s96vshs959overlap.bed"),wiggleTable,4,8,'Differential positions highlighted with t-test','S96', 'HS959', c("Overlap","Differential"),brewer.pal(3,"Dark2"),FALSE,FALSE)


x<-wiggleTable[ret$best$class==1,]

makeComparisonPlotHelp2(x,loadBed("s96vshs959overlap.bed"),wiggleTable,4,8,'Differential positions highlighted with DIME classifier','S96', 'HS959', c("Overlap","Differential","Background"),c(brewer.pal(7,"Spectral")[1:2],"#77777777"),FALSE,FALSE)



wiggleTableNorm<-normalizeBetweenArrays(as.matrix(wiggleTable[,3:6]),method="scale")


wiggleTableNorm<-cbind(wiggleTable[,c(1,2)],wiggleTableNorm)

wigglePeaksNorm<-getPeakScores(loadBed('data/s96overlap-high-peaks.bed'),wiggleTableNorm)
doheatmap(wigglePeaksNorm[,c(-1,-2)],Colv=TRUE,Rowv=TRUE)



#write out positions?!
write.table(cbind(wiggleTable[x$ID,1],wiggleTable[x$ID,2],wiggleTable[x$ID,2]+1),col.names=FALSE,row.names=FALSE,sep='\t')







t1<-wiggleTable[,4]
t2<-wiggleTable[,6]
plot(t1,t2)
myma<-list(M=t1-t2,A=((t1+t2))/2)
plotMA(myma)

lib.sizes=c(2218566,1292603,1106212,1300837)
voom1<-voom(wiggleTable[,c(4,6,8,10)],lib.size=lib.sizes)
exp<-voom1[['E']]
t1<-exp[,1]
t2<-exp[,2]
MyMAList<-list(M=treat-control,A=(treat+control)/2)
plotMA(MyMAList)


