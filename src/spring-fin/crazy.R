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


dimetab<-read.table('dimetable.txt')
head(dimetab)
pvals<-sort(unique(dimetab$V5))

pvalsequence<-seq(1,length(pvals))
for(i in pvalsequence) {
  print(pvals[i])
  write.table(file=paste0('mydimeoutpval',i,'.tab'),dimetab[dimetab$V5<pvals[i],],col.names=F,row.names=F,quote=F,sep='\t')
}

mergetab<-read.table('myoutputs.txt')
myx<-seq(100,90+10*nrow(mergetab),by=10)
plot(myx,mergetab[,1])
plot(exp(-unlist(pvals[1:147])),mergetab[1:147,1],xlab="P-value threshold", ylab="Number of peaks")
title('Number of differential peaks identified using p-value threshold')




pvals<-lapply(list.files(pattern=".rando"),function(filename)  {
  ret<-read.table(filename)
  #print(filename)
  print(max(as.numeric(ret[,4])))
  max(as.numeric(ret[,4]))
})

nsites<-lapply(list.files(pattern=".rando"),function(filename)  {
  ret<-read.table(filename)
  nrow(ret)
})

nsites2<-lapply(list.files(pattern="tab$"),function(filename)  {
  ret<-read.table(filename)
  print(nrow(ret))
  nrow(ret)
})

points(pvals,nsites)
title('Significance of ')






pvals<-lapply(list.files(pattern=".out"),function(filename)  {
  ret<-read.table(filename)
  #print(filename)
  print(min(as.numeric(ret[,4])))
  min(as.numeric(ret[,4]))
})

nsites<-lapply(list.files(pattern=".out"),function(filename)  {
  ret<-read.table(filename)
  nrow(ret)
})


nsites2<-lapply(list.files(pattern="tab$"),function(filename)  {
  ret<-read.table(filename)
  nrow(ret)
})


mergetab<-read.table('myoutputs.txt')
myx<-seq(100,90+10*nrow(mergetab),by=10)
plot(myx,mergetab[,1])
plot(exp(-unlist(pvals[1:147])),mergetab[1:147,1],xlab="P-value threshold", ylab="Number of peaks")
title('Number of differential peaks identified using p-value threshold')



par(mar=c(5,4,4,5)+.1)
plot(pvals[pvals<0.2],nsites[pvals<0.2],xlab='p-value cutoff',ylab='Number of peaks merged')
par(new=TRUE)
plot(pvals[pvals<0.2],nsites2[pvals<0.2],col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
axis(4)
mtext("Number of sites inputted",side=4,line=3)
title('Differential peaks identified using p-value threshold (iNUDGE)')
legend('bottomright',col=c('black','blue'),pch=1,legend=c('Merged peaks','Differential inputs'))


plot(exp(-unlist(pvals)),unlist(nsites2)/unlist(nsites),ylab="Merged peaks/Inputted sites",xlab="p-value")
title('Ratio of Merged peaks/Inputted sites for each p-value threshold')




axis(4)
mtext("Number of sites inputted",side=4,line=3)
title('Differential peaks identified using p-value threshold (iNUDGE)')
legend('bottomright',col=c('black','blue'),pch=1,legend=c('Merged peaks','Differential inputs'))
