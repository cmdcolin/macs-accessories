library(plyr)
library(gplots)





# Get chromosomes list
chrnames<-names(macswiggle[[1]]$treat)



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
  for(i in 2:length(macswiggle)) {
    ret<-getjoinscores(chrnames,ret,macswig[[i]]$control,currpos)
    ret<-getjoinscores(chrnames,ret,macswig[[i]]$treat,currpos+1)
    currpos<-currpos+2
  }
  flatten(chrnames,ret)
}



ret4<-joinWiggleFiles(macswiggle)

x<-sapply(names(macswiggle),function(x) strsplit(x,'-new')[[1]])
n<-length(macswiggle)
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
resize.win(10,20)




sret<-apply(ret4[,3:ncol(ret4)],2,function(x,g){slideMean(x,g,g)},1000)

heatmap.2(as.matrix(sret),col=redgreen(75), scale="none",key=TRUE, density.info='none',trace='none',Rowv=NULL)









#########################
# Fix plots




# Background subtraction and scaling
table<-ret4
strs2<-c('chr','pos',paste0(rep('V',n),1:(n*2)))
names(table)<-strs2


for(i in 1:length(macswiggle)) {
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


sret<-apply(table[,3:ncol(table)],2,function(x,g){slideMean(x,g,g)},1000)

heatmap.2(as.matrix(sret),col=redgreen(75), scale="none",key=TRUE, density.info='none',trace='none',Rowv=NULL)

sret<-apply(tablescale[,3:ncol(tablescale)],2,function(x,g){slideMean(x,g,g)},10000)

heatmap.2(as.matrix(sret),col=redgreen(75), scale="none",key=TRUE, density.info='none',trace='none',Rowv=NULL,Colv=NULL)



p1<-r1$V2-c1[cmatch,2]/m1
p2<-r2$V2-c2[cmatch2,2]/m2



#Background subtraction and normalization
v1<-mean(r1$V2)+mean(c1$V2[cmatch])/m1^2
v2<-mean(r2$V2)+mean(c2$V2[cmatch2])/m2^2
plot(p1/v1,p2[match]/v2,pch='.',cex=2)


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