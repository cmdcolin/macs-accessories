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



ret<-log2(slideMean(wiggleTable[,4],10,10)/slideMean(wiggleTable[,8],10,10))
system.time(dime1<-DIME(ret,gng.max.iter=500,inudge.max.iter=500,nudge.max.iter=500))

pal<-brewer.pal(4,"RdYlGn")

cex=0.6
png(filename="plot1.png",width=640,height=480)
plot(slideMean(wiggleTable[,4],10,10),slideMean(wiggleTable[,8],10,10),pch=19,cex=cex,col=pal[1],xlab="S96",ylab="HS959")
points(slideMean(wiggleTable[,4],10,10)[dime1$nudge$class==1],slideMean(wiggleTable[,8],10,10)[dime1$nudge$class==1],pch=19,cex=cex,col=paste0(pal[2],'77'))
title("S96 vs HS959 differential components (NUDGE)")
legend('bottomright',legend=c('Conserved','Differential'),fill=c(pal[1],pal[2]))

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
