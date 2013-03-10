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
    print(colnames(t2[[chr]]))
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

docol<-function(macswig) {
}
ret<-getjoinscores(chrnames,macswiggle[[1]]$control,macswiggle[[1]]$treat,3)
ret<-getjoinscores(chrnames,ret,macswiggle[[2]]$control,4)
ret<-getjoinscores(chrnames,ret,macswiggle[[2]]$treat,5)
ret<-getjoinscores(chrnames,ret,macswiggle[[3]]$control,6)
ret<-getjoinscores(chrnames,ret,macswiggle[[3]]$treat,7)
ret<-getjoinscores(chrnames,ret,macswiggle[[4]]$control,8)
ret<-getjoinscores(chrnames,ret,macswiggle[[4]]$treat,9)
ret<-getjoinscores(chrnames,ret,macswiggle[[5]]$control,10)
ret<-getjoinscores(chrnames,ret,macswiggle[[5]]$treat,11)
ret<-getjoinscores(chrnames,ret,macswiggle[[6]]$control,12)
ret<-getjoinscores(chrnames,ret,macswiggle[[6]]$treat,13)
ret<-getjoinscores(chrnames,ret,macswiggle[[7]]$control,14)
ret<-getjoinscores(chrnames,ret,macswiggle[[7]]$treat,15)
ret4<-flatten(chrnames,ret)


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