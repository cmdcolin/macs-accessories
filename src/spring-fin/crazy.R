set.seed(53535)
xValues = seq(0,2*pi,length=100)
yValues = rnorm(100) + sin(xValues)


lo<-sapply(1:10,function(i) {
  ns1<-ns(xValues,df=i)
  lm1<-lm(yValues~ns1)
  summary(lm1)$sigma
})
plot(lo)




set.seed(883833 )
quantile(Wind)


boot1<-one.boot(Wind,median,1000)


set.seed(7363 )




ret<-tail(d,650)
ret[,2]<-as.Date(ret[,2])
plot(ret[,2],ret[,8],type='l')


dimetab<-read.table('../mydimeout.txt')
head(dimetab)
pvalsequence<-seq(1,200,by=10)
pvals<-sort(dimetab$V5)
for(i in pvalsequence) {
  print(pvals[i])
  write.table(file=paste0('mydimeoutpval',i,'.tab'),dimetab[dimetab$V5<pvals[i],],col.names=F,row.names=F,quote=F,sep='\t')
}

mergetab<-read.table('myoutputs.txt')
myx<-seq(100,90+10*nrow(mergetab),by=10)
plot(myx,mergetab[,1])
plot(exp(-unlist(pvals[1:147])),mergetab[1:147,1],xlab="P-value threshold", ylab="Number of peaks")
title('Number of differential peaks identified using p-value threshold')




pvals<-lapply(list.files(pattern=".crazy"),function(filename)  {
  ret<-read.table(filename)
  #print(filename)
  print(max(as.numeric(ret[,4])))
  max(as.numeric(ret[,4]))
})

nsites<-lapply(list.files(pattern=".crazy"),function(filename)  {
  ret<-read.table(filename)
  nrow(ret)
})


points(pvals,nsites)
title('Significance of ')