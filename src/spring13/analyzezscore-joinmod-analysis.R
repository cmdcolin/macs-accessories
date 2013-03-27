################
# Main file for analyzezscore-joinmod.R 
#






## See histogram of average read depth from all experiments
dist1<-apply(wiggleTable[,c(-1,-2)],2,mean)
hist(dist1)



#resize for heatmap
resize.win(10,30)



# Conduct heatmap (omit chr and pos columns)
doheatmap(wiggleTable[,c(-1,-2)],1000)




# Background subtraction and scaling
for(i in 1:nsamples) {
  r1<-paste0('V',i+1)
  r2<-paste0('V',i+2)
  ratio=median(table[[r1]]/table[[r2]])
  if(ratio>1) {
    printf("scale down control %f\n", ratio)
    table[[r1]]<-table[[r1]]/ratio
  } else {
    printf("scale down treat %f\n", 1/ratio)
    table[[r2]]<-table[[r2]]*ratio
  }
}

tablescale<-table

for(i in 1:length(dist1)) {
  r1<-paste0('V',i)
  printf("Scale %s by mean total read depth %f\n", r1, dist1[i])
  
  tablescale[[r1]]<-table[[r1]]/dist1[i]
}


doheatmap(table[,3:ncol(ret4)],1000)
doheatmap(tablescale[,3:ncol(ret4)],10000)




table<-wiggleTable
backSubList<-lapply(1:nsamples,function(i) {
  pos=i*2
  str1<-paste0('V',pos-1)
  str2<-paste0('V',pos)
  control<-table[[str1]]
  treat<-table[[str2]]
  
  m1=median(control/treat)
  treat-control/m1
})

# combine rows into table
backSubTable<-as.data.frame(do.call(cbind,backSubList))

#get pos and chr columns
backSubTable<-with(wiggleTable, cbind(pos,chr,backSubTable))



doheatmap(normdifftable,1000)


makeComparisonPlot<-function(bed1,bed2,table,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale=FALSE) {
  
  bedOverlap<-intersectBed(bed1,bed2)
  bedUnique1<-uniqueBed(bed1,bed2)
  bedUnique2<-uniqueBed(bed2,bed1)
  makeComparisonPlotHelp(bedOverlap,bedUnique1,bedUnique2,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale)
}

makeComparisonPlotHelp<-function(bedOverlap,bedUnique1,bedUnique2,table,c1,c2,titlex,xlab,ylab,legendx,fillx,logscale=FALSE){
  

  wtoverlap<-getPeakScores(bedOverlap,table)
  wtrep1<-getPeakScores(bedUnique1,table)
  wtrep2<-getPeakScores(bedUnique2,table)
  p<-table[sample(1:nrow(table),nrow(table)/10),c(c1,c2)]
  p1<-wtoverlap[,c(c1,c2)]
  p2<-wtrep1[,c(c1,c2)]
  p3<-wtrep2[,c(c1,c2)]
  if(logscale==TRUE) {
    p<-log2(p)
    p1<-log2(p1)
    p2<-log2(p2)
    p3<-log2(p3)
  }
  plot(p1,pch=19,cex=0.9,col="#aaaaaa22",xlab=xlab,ylab=ylab)
  #,log=if(log)"xy"else""
  points(p1,pch=19,cex=0.9,col=paste0(fillx[1],"77"))
  points(p3,pch=19,cex=0.9,col=paste0(fillx[3],"77"))
  points(p2,pch=19,cex=0.9,col=paste0(fillx[2],"77"))
  legend("bottomright",fill=fillx,legend=legendx)
  title(titlex)
  if(logscale==TRUE){
    
    lines(1:9,(-1:7),lwd=3,col=fillx[3])
    lines(-1:7,(1:9),lwd=3,col=fillx[3])
  }
}




makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('s96rep2-high_peaks.bed'),wiggleTable,4,6,'Comparison of raw read scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set1"))
mycor=cor(wiggleTable[,4],wiggleTable[,6])
text(50,220,substitute(paste(rho,"=",mycor),list(mycor=mycor)))


makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('hs959rep1-new_peaks.bed'),wiggleTable,4,8,'Comparison of raw read scores for S96 vs HS959','S96', 'HS959', c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Dark2"))
mycor=cor(wiggleTable[,4],wiggleTable[,8])
text(40,250,substitute(paste(rho,"=",mycor),list(mycor=mycor)))
mypal<-brewer.pal(3,"Accent")







makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('s96rep2-high_peaks.bed'),normDiffTable,3,4,'Comparison of NormDiff scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set2"))
mycor=cor(normDiffTable[,3],normDiffTable[,4])
text(2,30,substitute(paste(rho,"=",mycor),list(mycor=mycor)))


makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('hs959rep1-new_peaks.bed'),normDiffTable,3,5,'Comparison of NormDiff scores for S96 vs HS959','S96','HS959',c("Overlap","S96","HS959"),brewer.pal(3,"BrBG"))
mycor=cor(normDiffTable[,4],normDiffTable[,8])
text(2,30,substitute(paste(rho,"=",mycor),list(mycor=mycor)))




makeComparisonPlot(loadBed('s96rep1-high_peaks.bed'),loadBed('s96rep2-high_peaks.bed'),backSubTable,3,4,'Comparison of BackSub scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set2"))


makeComparisonPlotHelp(loadBed('s96overlap-high-peaks.bed'),loadBed('hs959rep1-new_peaks.bed'),backSubTable,3,5,'Comparison of BackSub scores for S96 vs HS959','S96','HS959',c("Overlap","S96","HS959"),brewer.pal(3,"Spectral"))



######################3
# MACS DIFF
######################3
makeComparisonPlotHelp(loadBed('s96overlap-high-peaks.bed'),loadBed('S96rep1-Diff_peaks.bed'),loadBed('S96rep2-Diff_peaks.bed'),wiggleTable,4,6,'Comparison of S96 replicates','S96rep1','S96rep2',c("Overlap","S96rep1","S96rep2"),brewer.pal(3,"Spectral"))


makeComparisonPlotHelp(loadBed('s96vshs959overlap.bed'),loadBed('S96vsHS959-Diff_peaks.bed'),loadBed('HS959vsS96-Diff_peaks.bed'),wiggleTable,4,8,'Comparison of S96 vs HS959','S96rep1','HS959rep1',c("Overlap","S96rep1","HS959rep1"),brewer.pal(3,"Spectral"))
#plot(retbin,main="Hexagonal binning")

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
