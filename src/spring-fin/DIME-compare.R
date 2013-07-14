require(DIME)
require(RColorBrewer)
paste0<-function(...) paste(...,sep='')
load('wiggleTable.RData')
slideMean<-function(x,windowsize=100,slide=1){
  idx1<-seq(1,length(x),by=slide);
  idx1+windowsize->idx2;
  idx2[idx2>(length(x)+1)]<-length(x)+1;
  c(0,cumsum(x))->cx;
  return((cx[idx2]-cx[idx1])/windowsize);
}
sel<-my.classify$pdiff[my.classify$class==1]


ret<-log2(wiggleTableScale[,3]/wiggleTableScale[,7])
my.obj<-inudge.fit(ret,max.iter=50)
system.time(dime1<-inudge.classify(ret))
system.time(my.obj<-inudge.fit(ret,max.iter=10))
my.classify<-inudge.classify(ret,my.obj)
pal<-brewer.pal(4,"RdYlGn")

cache('my.classify')
cache('my.obj')
cache(dime1)

cex=0.6
png(filename="plot1.png",width=640,height=480)
plot(slideMean(wiggleTable[,4],10,10),slideMean(wiggleTable[,8],10,10),pch=19,cex=cex,col=pal[1],xlab="S96",ylab="HS959")
points(slideMean(wiggleTable[,4],10,10)[dime1$nudge$class==1],slideMean(wiggleTable[,8],10,10)[dime1$nudge$class==1],pch=19,cex=cex,col=paste0(pal[2],'77'))
title("S96 vs HS959 differential components (NUDGE)")
legend('bottomright',legend=c('Conserved','Differential'),fill=c(pal[1],pal[2]))


pal<-brewer.pal(4,"RdYlGn")
cex=0.7
plot(wiggleTableScale[,3],wiggleTableScale[,7],pch=19,cex=cex,col=pal[1],xlab="S96",ylab="HS959")
points(wiggleTableScale[,3][my.classify$class==1],wiggleTableScale[,7][my.classify$class==1],pch=19,cex=cex,col=pal[2])
title("S96 vs HS959 differential components (iNUDGE)")
legend('bottomright',legend=c('Conserved','Differential'),fill=c(pal[1],pal[2]))



my.classify$fdr[my.classify$class==1]

sel<-my.classify$fdr[my.classify$class==1]
sel2<-my.classify$fdr
order(sel)

hist((my.classify$pdiff[my.classify$class==1]))

dimetable<-data.frame(chr=wiggleTableScale[my.classify$class==1,1],
           start=wiggleTableScale[my.classify$class==1,2],
           end=wiggleTableScale[my.classify$class==1,2]+10,
           name=paste0('S',wiggleTableScale[my.classify$class==1,1],'B',wiggleTableScale[my.classify$class==1,2],'E',wiggleTableScale[my.classify$class==1,2]+10),
           pdiff=my.classify$pdiff[my.classify$class==1])


write.table(dimetable,col.names=F,row.names=F,quote=F,file='dimetable.txt',sep='\t')

sel<-sel[order(sel)]
dev.off()
png(filename="plot2.png",width=640,height=480)
plot(slideMean(wiggleTable[,4],10,10),slideMean(wiggleTable[,8],10,10),pch=19,cex=cex,col=pal[1],xlab="S96",ylab="HS959")
points(slideMean(wiggleTable[,4],10,10)[dime1$inudge$class==1],slideMean(wiggleTable[,8],10,10)[dime1$inudge$class==1],pch=19,cex=cex,col=paste0(pal[3],'77'))
title("S96 vs HS959 differential components (iNUDGE)")
legend('bottomright',legend=c('Conserved','Differential'),fill=c(pal[1],pal[3]))

dev.off()

png(filename="plot3.png",width=640,height=480)
plot(slideMean(wiggleTable[,4],10,10),slideMean(wiggleTable[,8],10,10),pch=19,cex=cex,col=pal[1],xlab="S96",ylab="HS959")
points(slideMean(wiggleTable[,4],10,10)[dime1$gng$class==1],slideMean(wiggleTable[,8],10,10)[dime1$gng$class==1],pch=19,cex=cex,col=paste0(pal[4],'77'))
title("S96 vs HS959 differential components (GNG)")
legend('bottomright',legend=c('Conserved','Differential'),fill=c(pal[1],pal[4]))
dev.off()
#points(slideMean(wiggleTable[,4],10,10)[dime1$inudge$class==1],slideMean(wiggleTable[,8],10,10)[dime1$inudge$class==1],pch=16,cex=cex,col=paste0(pal[3],'77'))
#points(slideMean(wiggleTable[,4],10,10)[dime1$gng$class==1],slideMean(wiggleTable[,8],10,10)[dime1$gng$class==1],pch=16,cex=cex,col=paste0(pal[4],'77'))




dev.off()



mytable<-cbind(as.character(wiggleTable[my.classify$class==1,1]),wiggleTable[my.classify$class==1,2],wiggleTable[my.classify$class==1,2]+10,paste0('mydiff',wiggleTable[my.classify$class==1,1],wiggleTable[my.classify$class==1,2]),my.classify$pdiff[my.classify$class==1])
write.table(mytable,quote=F,row.names=FALSE,col.names=FALSE,file='mydimeout.txt')
