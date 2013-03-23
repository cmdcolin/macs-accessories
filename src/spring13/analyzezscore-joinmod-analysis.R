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