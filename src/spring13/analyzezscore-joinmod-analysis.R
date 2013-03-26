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




table<-ret4
retlist<-lapply(1:nsamples,function(i) {
  pos=i*2
  str1<-paste0('V',pos)
  str2<-paste0('V',pos+1)
  control<-table[[str1]]
  treat<-table[[str2]]
  
  m1=median(control/treat)
  treat-control/m1
})


normdifflist<-lapply(1:nsamples,function(i) { 
  pos=i*2
  str1<-paste0('V',pos)
  str2<-paste0('V',pos+1)
  control<-wiggleTable[[str1]]
  treat<-wiggleTable[[str2]]
  getnormdiff(treat,control)
})



# combine rows into table
normdifftable<-as.data.frame(do.call(cbind,normdifflist))
#get pos and chr columns
normdifftable<-with(wiggleTable, cbind(pos,normdifftable))
normdifftable<-with(wiggleTable, cbind(chr,normdifftable))


doheatmap(normdifftable,1000)





bed1<-loadBed('s96rep1-high_peaks.bed')
bed2<-loadBed('s96rep2-high_peaks.bed')

nd1<-getPeakScores(bed1,normdifftable)
nd2<-getPeakScores(bed2,normdifftable)
wt1<-getPeakScores(bed1,wiggleTable)
wt2<-getPeakScores(bed2,wiggleTable)



require(RColorBrewer)
plot(wiggleTable[,c(4,6)],pch=20,cex=0.7)
points(wt1[,c(4,6)],pch=20,cex=0.7,col=paste0(x[2],"55"))
points(wt2[,c(4,6)],pch=20,cex=0.7,col=paste0(x[3],"55"))

plot(log(wiggleTable[,c(4,6)]),pch=20,cex=0.7)
points(log(wt1[,c(4,6)]),pch=20,cex=0.7,col=paste0(x[2],"55"))
points(log(wt2[,c(4,6)]),pch=20,cex=0.7,col=paste0(x[3],"55"))

plot(jitter(wiggleTable[,4]),jitter(wiggleTable[,6]),pch=20,cex=0.7)
points(jitter(wt1[,4]),jitter(wt1[,6]),pch=20,cex=0.7,col=paste0(x[2],"55"))
points(jitter(wt2[,4]),jitter(wt2[,6]),pch=20,cex=0.7,col=paste0(x[3],"55"))



makeComparisonPlot<-function(bedOverlap,bedUnique1,bedUnique2,table,c1,c2,titlex,xlab,ylab,legendx,fillx,log=FALSE){
  
  p<-table[sample(1:nrow(table),nrow(table)/10),c(c1,c2)]
  wtoverlap<-getPeakScores(bedOverlap,table)
  wtrep1<-getPeakScores(bedUnique1,table)
  wtrep2<-getPeakScores(bedUnique2,table)
  p1<-wtoverlap[,c(c1,c2)]
  p2<-wtrep1[,c(c1,c2)]
  p3<-wtrep2[,c(c1,c2)]
  if(log==TRUE) {
    p<-log2(p)
    p1<-log2(p1)
    p2<-log2(p2)
    p3<-log2(p3)
  }
  plot(p,pch=19,cex=0.9,col="#aaaaaa22",xlab=xlab,ylab=ylab)
  #,log=if(log)"xy"else""
  points(p1,pch=19,cex=0.9,col=paste0(fillx[1],"77"))
  points(p3,pch=19,cex=0.9,col=paste0(fillx[3],"77"))
  points(p2,pch=19,cex=0.9,col=paste0(fillx[2],"77"))
  legend("bottomright",fill=fillx,legend=legendx)
  title(titlex)
  
}



makeComparisonPlot(loadBed('s96overlap-unique-high-peaks.bed'),loadBed('s96rep1-unique-high-peaks.bed'),loadBed('s96rep2-unique-high-peaks.bed'),wiggleTable,4,6,'Comparison of raw read scores for S96 replicates','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set1"))
mycor=cor(wiggleTable[,4],wiggleTable[,6])
text(50,220,substitute(paste(rho,"=",mycor),list(mycor=mycor)))


makeComparisonPlot(loadBed('s96overlap-unique-high-peaks.bed'),loadBed('s96rep1-unique-high-peaks.bed'),loadBed('s96rep2-unique-high-peaks.bed'),wiggleTable,4,6,'Comparison of raw read scores for S96 replicates (log2)','S96rep1','S96rep2',c("Overlap","Rep1 unique","Rep2 unique"),brewer.pal(3,"Set1"),log=TRUE)
mycor=cor(log2(wiggleTable[,4]),log2(wiggleTable[,6]))
text(1,7,substitute(paste(rho,"=",mycor),list(mycor=mycor)))

makeComparisonPlot(loadBed('s96vshs959overlap.bed'),loadBed('s96vshs959unique.bed'),loadBed('s96vshs959unique2.bed'),wiggleTable,4,8,'Comparison of raw read scores for S96 vs HS959','S96', 'HS959', c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Dark2"))
mycor=cor(wiggleTable[,4],wiggleTable[,8])
text(40,250,substitute(paste(rho,"=",mycor),list(mycor=mycor)))

#try log scale
makeComparisonPlot(loadBed('s96vshs959overlap.bed'),loadBed('s96vshs959unique.bed'),loadBed('s96vshs959unique2.bed'),wiggleTable,4,8,'Comparison of raw read scores for S96 vs HS959 (log2)','S96', 'HS959',c("Overlap","S96 unique","HS959 unique"),brewer.pal(3,"Dark2"),log=TRUE)
mycor=cor(log2(wiggleTable[,4]),log2(wiggleTable[,8]))
text(1,8,substitute(paste(rho,"=",mycor),list(mycor=mycor)))
mypal<-brewer.pal(3,"Accent")
#lines(0:8,(0:8)-0,lwd=3,col=)
lines(1:9,(-1:7),lwd=3,col=mypal[3])
lines(-1:7,(1:9),lwd=3,col=mypal[3])

#lm1<-lm(log2(wiggleTable[,8])~log2(wiggleTable[,4]))
#lines(log2(wiggleTable[,4]),lm1$fitted)
#plot(wiggleTable[sample(1:nrow(wiggleTable),nrow(wiggleTable)/10),c(4,6)])


ex1<-loadBed('s96rep1-new_peaks.bed')
ex2<-loadBed('hs959rep1-new_peaks.bed')

nd1<-getPeakScores(bed1,normdifftable)
nd2<-getPeakScores(bed2,normdifftable)
nd3<-getPeakScores(bed3,normdifftable)

#retbin<-hexbin(wiggleTable[,c(4,6)])
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