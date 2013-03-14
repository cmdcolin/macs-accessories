library(plyr)
library(gplots)





# Get chromosomes list
chrnames<-names(macswiggle[[1]]$treat)


lenmacswiggle<-length(macswiggle)

### Get wig scores for treat only
getjoinscores<-function(chrnames,t1,t2,currpos) {
  
  ret<-lapply(chrnames, function(chr) {
    if(debug)
      printf("processing %s\tat V%d\n", chr,currpos);
    str<-paste0("V",currpos)
    colnames(t2[[chr]])[2]<-str
    ret<-join(t1[[chr]], t2[[chr]],by="V1")
    ret[complete.cases(ret),]
  })
  names(ret)<-chrnames
  #print(str(ret))
  # from R inferno, Burns (2011)
  #do.call('rbind', ret) 
  ret
}



flatten<-function(chrnames,ret) {
  RTE<-lapply(chrnames, function(chr) {
    if(debug)
      printf("processing %s\n", chr);
    data.frame(chr=chr,ret[[chr]])
  })
  
  do.call('rbind', RTE) 
}

joinWiggleFiles<-function(chrnames,macswig) {
  ret<-getjoinscores(chrnames,macswig[[1]]$control,macswig[[1]]$treat,3)
  currpos<-4
  for(i in 2:lenmacswiggle) {
    ret<-getjoinscores(chrnames,ret,macswig[[i]]$control,currpos)
    ret<-getjoinscores(chrnames,ret,macswig[[i]]$treat,currpos+1)
    currpos<-currpos+2
  }
  flatten(chrnames,ret)
}






doheatmap<-function(table,granularity=1) {
  
  if(granularity!=1) {
    sret<-apply(table[,1:ncol(table)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
  }
  else {
    sret<-table
  }
  
  heatmap.2(as.matrix(sret),col=redgreen(75), scale="none",key=TRUE, 
            density.info='none',trace='none',Rowv=FALSE)
}

doTreeView<-function(table,file=tempfile(),granularity=1) {
  
  if(granularity!=1) {
    lemm<-seq(1,nrow(table),by=granularity)
    coord<-paste0(table$chr[lemm],'P',table$pos[lemm])
    
    temp<-apply(table[,c(-1,-2)],2,function(x,g){
      slideMean(x,g,g)
    },granularity)
    sret<-cbind(coord,temp)
  }
  else {
    sret<-table
  }
  write.table(sret,file=file,row.names=FALSE,quote=FALSE,sep='\t')
}


getnormdiff<-function(treat,control) {
  
  treatSmooth<-slideMean(treat) #default params
  controlSmooth<-slideMean(control)
  
  
  scalingFactor=median(control/treat)
  globalVariance=sqrt(mean(treat)+mean(control)/scalingFactor^2)
  varlist<-sqrt(treatSmooth+controlSmooth/scalingFactor^2)
  
  
  #normdiff local
  normVar<-sapply(varlist,function(x){
    max(x,globalVariance)
  })
  
  (treat-control/scalingFactor)/normVar
}




getPeakNormDiff<-function(bed,scores) {
  chrmatch="NA"
  chrsub=data.frame()
  ret<-apply(bed,1,function(row) {
    printf("Processing peak %s (%d,%d)\n",row[4],row[2],row[3])
    chrselect<-strsplit(row[1],'.fsa')[[1]]
    if(chrmatch!=chrselect) {
      chrmatch<<-chrselect
      chrsub<<-scores[scores$chr==chrselect,]
    }
    chrsub[chrsub$pos>row[2]&chrsub$pos<row[3],]
  })
  
  # row bind dataframes. from R inferno, Burns (2011)
  do.call('rbind', ret) 
}





ret4<-joinWiggleFiles(chrnames, macswiggle)

x<-sapply(names(macswiggle),function(x) strsplit(x,'-new')[[1]])
n<-lenmacswiggle
strs<-paste0(sort(rep(x,2)),'-',rep(c('c','t'),n))
#strs2<-c('chr','pos',paste0(rep('V',n),1:(n*2)))
#names(ret4)<-strs2
names(ret4)<-c('chr','pos',strs)
head(ret4)

## See average read depth
dist1<-sapply(3:ncol(ret4),function(i) mean(ret4[,i]))
mean(dist1)


resize.win <- function(Width=6, Height=6)
{
  # works for windows
  dev.off(); # dev.new(width=6, height=6)
  windows(record=TRUE, width=Width, height=Height)
}
resize.win(10,30)




doheatmap(ret4[,3:ncol(ret4)],1000)








#########################
# Fix plots




# Background subtraction and scaling
table<-ret4
strs2<-c('chr','pos',paste0(rep('V',n),1:(n*2)))
names(table)<-strs2


for(i in 1:lenmacswiggle) {
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