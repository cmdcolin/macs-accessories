################
# Main file for analyzezscore-joinmod.R 
#


# Set wiggle table names as abbreviated sample names
prettyNames<-function(wiggleTable) {
  
  fixNames<-sapply(names(macswiggle),function(x) strsplit(x,'-new')[[1]])
  prettyNames<-paste0(sort(rep(fixNames,2)),'-',rep(c('c','t'),nsamples))
  names(wiggleTable)<-c('chr','pos',prettyNames)
}

## Set wiggle table names as chr, pos, V1-Vn
plainNames<-function(wiggleTable) {
  
  fixNames<-c('chr','pos',paste0(rep('V',nsamples),1:(nsamples*2)))
  names(wiggleTable)<-fixNameszx
  
}


# Join wiggle files with matching positions into a table
wiggleTable<-joinWiggleFiles(chrnames, macswiggle)

head(wiggleTable)



## See histogram of average read depth from all experiments
dist1<-apply(wiggleTable[,c(-1,-2)],2,mean)
hist(dist1)



# Accessory function for plotting large heatmaps
resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
resize.win(10,30)



# Conduct heatmap (omit chr and pos columns)
doheatmap(wiggleTable[,c(-1,-2)],1000)








#########################
# Fix plots




# Background subtraction and scaling
table<-wiggleTable
table<-plainNames(table)


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


dist1<-sapply(names(table)[3:ncol(table)],function(i) mean(table[[i]]))
mean(dist1)
tablescale<-table

for(i in 1:length(dist1)) {
  r1<-paste0('V',i)
  printf("Scale %s by mean total read depth %f\n", r1, dist1[i])
  
  tablescale[[r1]]<-table[[r1]]/dist1[i]
}


doheatmap(table[,3:ncol(ret4)],1000)
doheatmap(tablescale[,3:ncol(ret4)],10000)




table<-ret4
retlist<-lapply(1:lenmacswiggle,function(i) {
  pos=i*2
  str1<-paste0('V',pos)
  str2<-paste0('V',pos+1)
  control<-table[[str1]]
  treat<-table[[str2]]
  
  m1=median(control/treat)
  treat-control/m1
})


normdifflist<-lapply(1:lenmacswiggle,function(i) { 
  pos=i*2
  str1<-paste0('V',pos)
  str2<-paste0('V',pos+1)
  control<-table[[str1]]
  treat<-table[[str2]]
  getnormdiff(treat,control)
})




caca<-as.data.frame(do.call(cbind,retlist))
names(caca)<-names(table)[seq(3,ncol(table),by=2)]
doheatmap(caca,10000)

caca2<-as.data.frame(do.call(cbind,normdifflist))
caca2<-cbind(ret4$pos,caca2)
caca2<-cbind(ret4$chr,caca2)
names(caca2)<-names(table)[seq(3,ncol(table),by=2)]


doheatmap(caca2,1000)





bed1<-loadBed('s96rep1-high_peaks.bed')
bed2<-loadBed('hs959rep1-new_peaks.bed')

nd1<-getPeakNormDiff(bed1,caca2)


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